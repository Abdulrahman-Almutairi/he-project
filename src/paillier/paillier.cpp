#include "paillier/paillier.h"
#include "bigint/mod_arith.h"
#include "bigint/prime.h"
#include "bigint/inv.h"
#include "bigint/div_exact.h"
#include <iostream>
#include <random>
#include <stdexcept>
#include <algorithm>

namespace he
{
  using bi::BigInt;
  using bi::gcd_u64;
  using bi::lcm_u64;
  using bi::modinv_u64;
  using bi::mul_mod;
  using bi::mul_mod_u64;
  using bi::pow_mod;
  using bi::pow_mod_bigexp;
  using bi::pow_mod_u64;
  using bi::to_u64;

  static BigInt to_big(std::uint64_t x) { return BigInt::from_u64(x); }

  KeyPair keygen(unsigned bits)
  {
    using namespace bi;

    if (bits < 128)
      bits = 128; // safety for demo

    KeyPair kp;

    std::random_device rd;
    std::mt19937_64 rng(rd());

    // 1) generate big primes p,q
    BigInt p = generate_prime(bits / 2, rng, 20);
    BigInt q = generate_prime(bits / 2, rng, 20);
    while (cmp(p, q) == 0)
    {
      q = generate_prime(bits / 2, rng, 20);
    }

    // 2) n = p*q, n^2 = n*n
    BigInt n = mul_schoolbook(p, q);
    BigInt n2 = mul_schoolbook(n, n);

    // 3) g = n + 1
    BigInt g = n;
    // g = bi::add_u64(g, 1); // if you don't have add_u64 globally, see note below
    if (g.limbs.empty())
      g.limbs.push_back(1);
    else
    {
      std::uint64_t carry = 1;
      for (std::size_t i = 0; i < g.limbs.size() && carry; ++i)
      {
        unsigned __int128 s = (unsigned __int128)g.limbs[i] + carry;
        g.limbs[i] = (std::uint64_t)s;
        carry = (std::uint64_t)(s >> 64);
      }
      if (carry)
        g.limbs.push_back(carry);
    }
    g.normalize();

    // 4) lambda = (p-1)(q-1)  (works with g=n+1 trick)
    BigInt pm1 = p;
    sub_inplace(pm1, BigInt::from_u64(1));
    BigInt qm1 = q;
    sub_inplace(qm1, BigInt::from_u64(1));
    BigInt lambda = mul_schoolbook(pm1, qm1);

    // 5) mu = lambda^{-1} mod n
    BigInt mu = inv_mod_odd(lambda, n);

    kp.pk.n = n;
    kp.pk.n2 = n2;
    kp.pk.g = g;

    kp.sk.lambda = lambda;
    kp.sk.mu = mu;

    return kp;
  }

  bi::BigInt encrypt(const PublicKey &pk, std::uint64_t m)
  {
    // cipher = g^m * r^n mod n^2
    // 1) g^m mod n^2
    BigInt gm = pow_mod(pk.g, m, pk.n2);
    // 2) choose random r in [1, n-1] (dummy: 64-bit, no gcd check yet)
    // BigInt r: random in [1, n-1], must be coprime with n.
    // For now we generate random odd r < n and skip gcd check (acceptable for development);
    // later we add BigInt gcd and enforce gcd(r,n)=1.

    std::random_device rd;
    std::mt19937_64 gen(rd());

    // random r < n (rejection sampling by limbs)
    auto rand_below = [&](const BigInt &limit) -> BigInt
    {
      BigInt L = limit;
      L.normalize();
      BigInt x;
      x.limbs.resize(L.limbs.size());
      while (true)
      {
        for (std::size_t i = 0; i < x.limbs.size(); ++i)
          x.limbs[i] = gen();
        x.normalize();
        if (bi::cmp(x, L) < 0 && !x.limbs.empty())
          return x;
      }
    };

    BigInt one = BigInt::from_u64(1);
    BigInt n_minus_1 = pk.n;
    bi::sub_inplace(n_minus_1, one);

    BigInt r = rand_below(n_minus_1);
    r.limbs[0] |= 1ULL; // make it odd
    r.normalize();

    // 3) r^n mod n^2
    // std::uint64_t n_u64 = to_u64(pk.n);
    BigInt rn = pow_mod_bigexp(r, pk.n, pk.n2);
    // 4) multiply mod n^2
    BigInt c = bi::mul_mod(gm, rn, pk.n2);
    return c;
  }

  bi::BigInt decrypt(const PublicKey &pk, const SecretKey &sk, const BigInt &c)
  {
    // u = c^lambda mod n^2
    BigInt u = bi::pow_mod_bigexp(c, sk.lambda, pk.n2);

    // L(u) = (u - 1) / n
    BigInt u_minus_1 = u;
    bi::sub_inplace(u_minus_1, BigInt::from_u64(1));
    BigInt L = bi::div_exact(u_minus_1, pk.n);

    // m = L * mu mod n
    BigInt m = bi::mul_mod(L, sk.mu, pk.n);
    return m;
  }

  std::uint64_t decrypt_u64(const PublicKey &pk, const SecretKey &sk, const BigInt &c)
  {
    BigInt m = decrypt(pk, sk, c);
    return bi::to_u64(m); // safe if your plaintext fits in u64 (your tests do)
  }

  bi::BigInt add(const PublicKey &pk, const bi::BigInt &c1, const bi::BigInt &c2)
  {
    return bi::mul_mod(c1, c2, pk.n2);
  }

  bi::BigInt scale(const PublicKey &pk, const bi::BigInt &c, std::uint64_t k)
  {
    return pow_mod(c, k, pk.n2);
  }

  bi::BigInt enc_dot_cp(const PublicKey &pk, const std::vector<bi::BigInt> &enc_a, const std::vector<std::uint64_t> &b)
  {
    if (enc_a.size() != b.size())
    {
      return bi::BigInt::from_u64(0);
    }

    bi::BigInt acc = bi::BigInt::from_u64(1);

    for (std::size_t i = 0; i < enc_a.size(); ++i)
    {
      // Enc(a_i * b_i) = Enc(a_i) ^ b_i
      auto term = he::scale(pk, enc_a[i], b[i]);
      // Accumulate sum via multiplication in ciphertext space
      acc = he::add(pk, acc, term);
    }

    return acc;
  }

  std::vector<bi::BigInt> enc_gemv_cp(const PublicKey &pk,
                                      const std::vector<bi::BigInt> &enc_A,
                                      std::size_t rows,
                                      std::size_t cols,
                                      const std::vector<std::uint64_t> &x)
  {
    std::vector<bi::BigInt> y;

    // shape checks
    if (x.size() != cols || enc_A.size() != rows * cols)
    {
      return y; // empty indicates mismatch
    }

    // Pre-size output so each thread writes to y[i] safely
    y.resize(rows);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < rows; ++i)
    {
      // Enc(0) identity under ciphertext multiplication is 1 mod n^2
      bi::BigInt acc = bi::BigInt::from_u64(1);

      const std::size_t row_off = i * cols;
      for (std::size_t j = 0; j < cols; ++j)
      {
        // Enc(Aij * xj) = Enc(Aij) ^ xj
        bi::BigInt term = he::scale(pk, enc_A[row_off + j], x[j]);

        // Accumulate: Enc(sum) = Enc(sum) * Enc(term)
        acc = he::add(pk, acc, term);
      }

      y[i] = acc;
    }

    return y;
  }

  std::vector<bi::BigInt> enc_gemm_cp(const PublicKey &pk,
                                      const std::vector<bi::BigInt> &enc_A,
                                      std::size_t A_rows,
                                      std::size_t A_cols,
                                      const std::vector<std::uint64_t> &B,
                                      std::size_t B_cols)
  {
    // shape checks
    if (enc_A.size() != A_rows * A_cols)
      return {};
    if (B.size() != A_cols * B_cols)
      return {};

    // Output C (row-major)
    std::vector<bi::BigInt> enc_C(A_rows * B_cols);

    // --- block sizes (tune later) ---
    const std::size_t Bi = 16; // rows tile
    const std::size_t Bk = 8;  // cols tile
    const std::size_t Bj = 16; // inner (shared) dimension tile

// Parallelize across (i0,k0) tiles (independent)
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (std::size_t i0 = 0; i0 < A_rows; i0 += Bi)
    {
      for (std::size_t k0 = 0; k0 < B_cols; k0 += Bk)
      {

        const std::size_t i_max = std::min(A_rows, i0 + Bi);
        const std::size_t k_max = std::min(B_cols, k0 + Bk);

        const std::size_t ti = i_max - i0;
        const std::size_t tk = k_max - k0;

        // Accumulators for this output tile.
        // Enc(0) identity under ciphertext multiplication is 1 mod n^2.
        std::vector<bi::BigInt> acc(ti * tk, bi::BigInt::from_u64(1));

        // Iterate over the shared dimension in blocks
        for (std::size_t j0 = 0; j0 < A_cols; j0 += Bj)
        {
          const std::size_t j_max = std::min(A_cols, j0 + Bj);

          // Multiply-accumulate into this tile
          for (std::size_t i = i0; i < i_max; ++i)
          {
            const std::size_t a_row = i * A_cols;

            for (std::size_t j = j0; j < j_max; ++j)
            {
              const bi::BigInt &enc_aij = enc_A[a_row + j];
              const std::size_t b_row = j * B_cols;

              for (std::size_t k = k0; k < k_max; ++k)
              {
                const std::uint64_t b_jk = B[b_row + k];

                // Enc(Aij * Bjk) = Enc(Aij) ^ Bjk
                bi::BigInt term = he::scale(pk, enc_aij, b_jk);

                // Enc(sum + term) = Enc(sum) * Enc(term)
                const std::size_t ai = i - i0;
                const std::size_t ak = k - k0;
                acc[ai * tk + ak] = he::add(pk, acc[ai * tk + ak], term);
              }
            }
          }
        }

        // Store the tile to output
        for (std::size_t i = i0; i < i_max; ++i)
        {
          for (std::size_t k = k0; k < k_max; ++k)
          {
            const std::size_t ai = i - i0;
            const std::size_t ak = k - k0;
            enc_C[i * B_cols + k] = acc[ai * tk + ak];
          }
        }
      }
    }

    return enc_C;
  }

}
