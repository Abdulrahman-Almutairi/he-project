#pragma once
#include "bfv/params.h"
#include "bfv/keys.h"
#include "bfv/ciphertext.h"
#include "bfv/ciphertext3.h"
#include <random>

namespace bfv
{
  // Key generation
  KeyPair keygen(const Params &p, std::mt19937_64 &rng);

  // Helpers: sampling
  Poly sample_uniform_q(const Params &p, std::mt19937_64 &rng);  // coeffs in [0,q)
  Poly sample_ternary(const Params &p, std::mt19937_64 &rng);    // coeffs in {-1,0,1} mod q
  Poly sample_small_error(const Params &p, std::mt19937_64 &rng); // small noise (for now ternary)

  // Encode integer m (0..t-1) as constant polynomial in R_t
  Poly encode_const(const Params &p, std::uint64_t m);

  // Decode constant coefficient back to integer in [0,t)
  std::uint64_t decode_const(const Params &p, const Poly &m);

  Ciphertext encrypt_const(const Params &p, const PublicKey &pk, std::uint64_t m, std::mt19937_64 &rng);
  Ciphertext encrypt_const_noiseless(const Params &p, const PublicKey &pk, std::uint64_t m, std::mt19937_64 &rng);
  std::uint64_t decrypt_const(const Params &p, const SecretKey &sk, const Ciphertext &ct);

  Ciphertext add(const Params &p, const Ciphertext &x, const Ciphertext &y);

  Ciphertext3 mul_raw(const Params &p, const Ciphertext &x, const Ciphertext &y);

  // Divide all coefficients by Delta= floor(q/t) with nearest rounding mod q.
  // This is the toy BFV "scale-down" step after multiplication.
  Ciphertext3 rescale_delta(const Params &p, const Ciphertext3 &ct3);

  // Convenience: full multiply = mul_raw + rescale
  Ciphertext3 mul(const Params &p, const Ciphertext &x, const Ciphertext &y);

  KeyPair keygen_noiseless(const Params &p, std::mt19937_64 &rng);
}
