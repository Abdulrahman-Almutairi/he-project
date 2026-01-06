#pragma once
#include "bfv/bfv.h"
#include <vector>
#include <random>


namespace bfv
{
  // Relinearization key: rk[i] encrypts (B^i * s^2) under secret key s
  struct RelinKey
  {
    std::uint64_t base = 4;          // gadget base B
    std::size_t L = 0;              // number of digits
    std::vector<Ciphertext> rk;     // rk[i] is a ciphertext (c0,c1)
  };

  // Build noiseless relin key (for correctness testing)
  RelinKey relin_keygen_noiseless(const Params &p, const SecretKey &sk, std::mt19937_64 &rng,
                                  std::uint64_t base = 4);

  // Decompose a polynomial c2 into digits polys d[i] in base B:
  // c2 = sum_i d[i] * B^i  (coeff-wise in Z_q)
  std::vector<Poly> decompose_base_B(const Params &p, const Poly &c2, std::uint64_t base, std::size_t L);

  // Relinearize: (c0,c1,c2) -> (c0',c1')
  Ciphertext relinearize(const Params &p, const Ciphertext3 &ct3, const RelinKey &rk);
}
