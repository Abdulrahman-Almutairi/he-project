#include "paillier/paillier.h"
#include "bigint/mod_arith.h"
#include <iostream>
#include <random>
#include <stdexcept>

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
    if (x.size() != cols || enc_A.size() != rows * cols)
    {
      return y; // empty indicates mismatch
    }

    y.reserve(rows);

    for (std::size_t i = 0; i < rows; ++i)
    {
      // Enc(0) identity under ciphertext multiplication is 1 mod n^2
      bi::BigInt acc = bi::BigInt::from_u64(1);

      const std::size_t row_off = i * cols;
      for (std::size_t j = 0; j < cols; ++j)
      {
        // term = Enc(A[i,j] * x[j]) = Enc(A[i,j]) ^ x[j]
        bi::BigInt term = he::scale(pk, enc_A[row_off + j], x[j]);
        // acc = Enc(sum + Aij*xj)  (ciphertext multiply)
        acc = he::add(pk, acc, term);
      }

      y.push_back(acc);
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
    std::vector<bi::BigInt> enc_C;
    enc_C.reserve(A_rows * B_cols);

    // shape checks
    if (enc_A.size() != A_rows * A_cols)
      return {};
    if (B.size() != A_cols * B_cols)
      return {};

    // For each output cell C[i,k] = sum_j A[i,j] * B[j,k]
    for (std::size_t i = 0; i < A_rows; ++i)
    {
      for (std::size_t k = 0; k < B_cols; ++k)
      {
        // Accumulate encrypted sum using ciphertext multiplication identity Enc(0)=1
        bi::BigInt acc = bi::BigInt::from_u64(1);

        for (std::size_t j = 0; j < A_cols; ++j)
        {
          const auto &enc_aij = enc_A[i * A_cols + j];
          std::uint64_t b_jk = B[j * B_cols + k];

          auto term = he::scale(pk, enc_aij, b_jk); // Enc(Aij * Bjk)
          acc = he::add(pk, acc, term);             // Enc(sum + term)
        }

        enc_C.push_back(acc);
      }
    }
    return enc_C;
  }
}
