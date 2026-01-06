#include "bfv/bfv.h"
#include "bfv/eval.h"
#include "bfv/relin.h"
#include <cassert>
#include <iostream>
#include <random>

int main()
{
  // Use your "real-ish" NTT params
  bfv::Params P{1024, 12289, 2};
  std::mt19937_64 rng(123);

  auto kp = bfv::keygen(P, rng);
  auto rlk = bfv::relin_keygen_noiseless(P, kp.sk, rng, 4);

  // 2x3 * 3x2 example over Z_t (t=2)
  const std::size_t A_rows=2, A_cols=3, B_cols=2;

  // Use small values mod t
  std::uint64_t A_plain[] = {1,0,1,
                             1,1,0}; // mod 2
  std::uint64_t B_plain[] = {1,1,
                             0,1,
                             1,0};  // mod 2

  // Encrypt A and B entrywise
  std::vector<bfv::Ciphertext> enc_A;
  std::vector<bfv::Ciphertext> enc_B;
  enc_A.reserve(A_rows*A_cols);
  enc_B.reserve(A_cols*B_cols);

  for (std::size_t i=0;i<A_rows*A_cols;++i)
    enc_A.push_back(bfv::encrypt_const(P, kp.pk, A_plain[i] % P.t, rng));

  for (std::size_t i=0;i<A_cols*B_cols;++i)
    enc_B.push_back(bfv::encrypt_const(P, kp.pk, B_plain[i] % P.t, rng));

  // Homomorphic GEMM (ct×ct)
  auto enc_C = bfv::enc_gemm_cc(P, enc_A, A_rows, A_cols, enc_B, B_cols, rlk);
  assert(enc_C.size() == A_rows*B_cols);

  // Decrypt and compare with plaintext GEMM mod t
  auto idx = [&](std::size_t r, std::size_t c){ return r*B_cols + c; };

  for (std::size_t i=0;i<A_rows;++i)
  {
    for (std::size_t k=0;k<B_cols;++k)
    {
      std::uint64_t sum = 0;
      for (std::size_t j=0;j<A_cols;++j)
      {
        sum += (A_plain[i*A_cols+j] * B_plain[j*B_cols+k]);
      }
      sum %= P.t;

      std::uint64_t dec = bfv::decrypt_const(P, kp.sk, enc_C[idx(i,k)]);
      assert(dec == sum);
    }
  }

  std::cout << "BFV encrypted GEMM (ct×ct) passed.\n";
  return 0;
}
