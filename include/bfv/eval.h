#pragma once
#include "bfv/bfv.h"
#include "bfv/relin.h"
#include <cstddef>
#include <vector>

namespace bfv
{
  Ciphertext add(const Params &p, const Ciphertext &x, const Ciphertext &y);

  // Multiply + relinearize (returns 2-part ciphertext)
  Ciphertext mul(const Params &p, const Ciphertext &x, const Ciphertext &y, const RelinKey &rlk);

  // Helpers for matrix multiply (each entry is its own ciphertext)
  std::vector<Ciphertext> enc_gemm_cc(const Params &p,
                                      const std::vector<Ciphertext> &enc_A,
                                      std::size_t A_rows,
                                      std::size_t A_cols,
                                      const std::vector<Ciphertext> &enc_B,
                                      std::size_t B_cols,
                                      const RelinKey &rlk);
}
