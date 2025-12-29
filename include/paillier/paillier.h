#pragma once
#include <vector>
#include <cstdint>
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
  bi::BigInt decrypt(const PublicKey &pk, const SecretKey &sk, const bi::BigInt &c);
  std::uint64_t decrypt_u64(const PublicKey &pk, const SecretKey &sk, const bi::BigInt &c); // placeholder

  bi::BigInt add(const PublicKey &pk, const bi::BigInt &c1, const bi::BigInt &c2); // homomorphic addition
  bi::BigInt scale(const PublicKey &pk, const bi::BigInt &c, std::uint64_t k);

  bi::BigInt enc_dot_cp(const PublicKey &pk, const std::vector<bi::BigInt> &enc_a, const std::vector<std::uint64_t> &b);
  std::vector<bi::BigInt> enc_gemv_cp(const PublicKey &pk, 
                                      const std::vector<bi::BigInt> &enc_A, 
                                      std::size_t rows, 
                                      std::size_t cols, 
                                      const std::vector<std::uint64_t> &x);

  std::vector<bi::BigInt> enc_gemm_cp(const PublicKey& pk,
                                      const std::vector<bi::BigInt>& enc_A,
                                      std::size_t A_rows,
                                      std::size_t A_cols,
                                      const std::vector<std::uint64_t>& B,
                                      std::size_t B_cols);
}
