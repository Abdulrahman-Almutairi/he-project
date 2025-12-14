#include "paillier/paillier.h"
#include "bigint/mod_arith.h"
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
  using bi::pow_mod_u64;
  using bi::to_u64;

  static BigInt to_big(std::uint64_t x) { return BigInt::from_u64(x); }

  KeyPair keygen(unsigned bits)
  {
    (void)bits;
    KeyPair kp;

    std::uint64_t p = 50021; // prime
    std::uint64_t q = 50023; // prime
    std::uint64_t n_u64 = p * q;
    std::uint64_t n2_u64 = (std::uint64_t)((unsigned __int128)n_u64 * (unsigned __int128)n_u64);

    // Public key
    kp.pk.n = to_big(n_u64);
    kp.pk.n2 = to_big(n2_u64);
    kp.pk.g = to_big(n_u64 + 1);

    // Private key: lambda and mu
    std::uint64_t lambda_u64 = lcm_u64(p - 1, q - 1);

    // compute L (g ^ lambda mod n^2) = (u - 1) / n
    std::uint64_t u = pow_mod_u64(n_u64 + 1, lambda_u64, n2_u64);
    std::uint64_t L = (u - 1) / n_u64;

    std::uint64_t mu_u64 = 0;
    bool ok = modinv_u64(L % n_u64, n_u64, mu_u64);
    if (!ok)
      throw std::runtime_error("modinv failed in keygen (toy)");

    kp.sk.lambda = to_big(lambda_u64);
    kp.sk.mu = to_big(mu_u64);

    return kp;
  }

  bi::BigInt encrypt(const PublicKey &pk, std::uint64_t m)
  {
    // cipher = g^m * r^n mod n^2
    // 1) g^m mod n^2
    BigInt gm = pow_mod(pk.g, m, pk.n2);
    // 2) choose random r in [1, n-1] (dummy: 64-bit, no gcd check yet)
    std::uint64_t n_u64 = to_u64(pk.n);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<std::uint64_t> dist(1, n_u64 - 1);
    std::uint64_t r_u64 = 0;
    do
    {
      r_u64 = dist(gen);
    } while (gcd_u64(r_u64, n_u64) != 1);
    BigInt r = BigInt::from_u64(r_u64);

    // 3) r^n mod n^2
    // std::uint64_t n_u64 = to_u64(pk.n);
    BigInt rn = pow_mod(r, n_u64, pk.n2);
    // 4) multiply mod n^2
    BigInt c = bi::mul_mod(gm, rn, pk.n2);
    return c;
  }

  std::uint64_t decrypt_u64(const PublicKey &pk, const SecretKey &sk, const BigInt &c)
  {
    std::uint64_t n_u64 = to_u64(pk.n);
    std::uint64_t n2_u64 = to_u64(pk.n2);
    std::uint64_t lambda = to_u64(sk.lambda);
    std::uint64_t mu = to_u64(sk.mu);

    // u = c^lambda mod n^2
    // (since n^2 fits in u64, convert to u64 fast path cleanly)
    std::uint64_t c_u64 = to_u64(c);
    std::uint64_t u = pow_mod_u64(c_u64 % n2_u64, lambda, n2_u64);

    // L(u) = (u - 1) / n
    std::uint64_t L = (u + n2_u64 - 1) % n2_u64; // safe (avoid underflow), but u>=1 anyway
    L = L / n_u64;

    // m = (L * mu) mod n
    std::uint64_t m = mul_mod_u64(L % n_u64, mu, n_u64);
    return m;
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
