#pragma once
#include "bigint/bigint.h"

namespace he
{
  struct PublicKey
  {
    bi::BigInt n;
    bi::BigInt n2;
    bi::BigInt g;
  };
  struct SecretKey
  {
    bi::BigInt lambda;
    bi::BigInt mu;
  };
  struct KeyPair
  {
    PublicKey pk;
    SecretKey sk;
  };

  KeyPair keygen(unsigned bits);                                                            // placeholder
  bi::BigInt encrypt(const PublicKey &pk, std::uint64_t m);                                 // placeholder
  std::uint64_t decrypt_u64(const PublicKey &pk, const SecretKey &sk, const bi::BigInt &c); // placeholder

  bi::BigInt add(const PublicKey &pk, const bi::BigInt &c1, const bi::BigInt &c2); // homomorphic addition
  bi::BigInt scale(const PublicKey &pk, const bi::BigInt &c, std::uint64_t k);
}
