#include "bfv/eval.h"
#include <stdexcept>

namespace bfv
{
  std::vector<Ciphertext> enc_gemm_cc(const Params &p,
                                      const std::vector<Ciphertext> &enc_A,
                                      std::size_t A_rows,
                                      std::size_t A_cols,
                                      const std::vector<Ciphertext> &enc_B,
                                      std::size_t B_cols,
                                      const RelinKey &rlk)
  {
    if (enc_A.size() != A_rows * A_cols) return {};
    if (enc_B.size() != A_cols * B_cols) return {};

    std::vector<Ciphertext> enc_C;
    enc_C.reserve(A_rows * B_cols);

    for (std::size_t i = 0; i < A_rows; ++i)
    {
      for (std::size_t k = 0; k < B_cols; ++k)
      {
        // acc = Enc(0)  (in our encoding Enc(0) is just encrypt_const(...,0))
        // but we don't want randomness in accumulation here; simplest is:
        // start with first term, then add remaining.
        bool first = true;
        Ciphertext acc(p);

        for (std::size_t j = 0; j < A_cols; ++j)
        {
          const auto &a = enc_A[i * A_cols + j];
          const auto &b = enc_B[j * B_cols + k];

          Ciphertext prod = mul(p, a, b, rlk); // ct×ct + relin

          if (first)
          {
            acc = prod;
            first = false;
          }
          else
          {
            acc = add(p, acc, prod);
          }
        }

        enc_C.push_back(acc);
      }
    }

    return enc_C;
  }
}
