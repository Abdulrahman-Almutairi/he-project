#include "paillier/paillier.h"
#include "bigint/mod_arith.h"
#include <iostream>
#include <random>
#include <stdexcept>

namespace he
{
  using bi::BigInt;
  using bi::to_u64;
  using bi::pow_mod_u64;
  using bi::mul_mod_u64;
  using bi::gcd_u64;
  using bi::lcm_u64;
  using bi::modinv_u64;
  using bi::pow_mod;
  using bi::mul_mod;

  static BigInt to_big(std::uint64_t x) { return BigInt::from_u64(x); }

  KeyPair keygen(unsigned bits)
  {
    (void)bits;
    KeyPair kp;

    std::uint64_t p = 50021;  // prime
    std::uint64_t q = 50023;  // prime
    std::uint64_t n_u64 = p * q;
    std::uint64_t n2_u64 = (std::uint64_t)((unsigned __int128)n_u64 * (unsigned __int128)n_u64);

    // Public key
    kp.pk.n  = to_big(n_u64);
    kp.pk.n2 = to_big(n2_u64);
    kp.pk.g  = to_big(n_u64 + 1);

    // Private key: lambda and mu
    std::uint64_t lambda_u64 = lcm_u64(p - 1, q - 1);

    // compute L (g ^ lambda mod n^2) = (u - 1) / n
    std::uint64_t u = pow_mod_u64(n_u64 + 1, lambda_u64, n2_u64);
    std::uint64_t L = (u - 1) / n_u64;

    std::uint64_t mu_u64 = 0;
    bool ok = modinv_u64(L % n_u64, n_u64, mu_u64);
    if (!ok) throw std::runtime_error("modinv failed in keygen (toy)");

    kp.sk.lambda = to_big(lambda_u64);
    kp.sk.mu     = to_big(mu_u64);

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
    do { r_u64 = dist(gen); } while (gcd_u64(r_u64, n_u64) != 1);
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
    std::uint64_t n_u64  = to_u64(pk.n);
    std::uint64_t n2_u64 = to_u64(pk.n2);
    std::uint64_t lambda = to_u64(sk.lambda);
    std::uint64_t mu     = to_u64(sk.mu);
    
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

  bi:: BigInt add(const PublicKey &pk, const bi::BigInt &c1, const bi::BigInt &c2) {
    return bi::mul_mod(c1, c2, pk.n2);
  }

  bi::BigInt scale(const PublicKey &pk, const bi::BigInt &c, std::uint64_t k) {
    return pow_mod(c, k, pk.n2);
  }
}
