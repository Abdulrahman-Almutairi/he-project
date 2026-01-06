#pragma once
#include "bfv/bfv.h"
#include <cstdint>

namespace bfv
{
  struct NoiseStats
  {
    std::int64_t inf = 0;          // max |e[i]|
    std::uint64_t l2_sq = 0;       // sum e[i]^2
  };

  // Noise of a 2-part ciphertext assuming it encrypts constant m (in Z_t).
  NoiseStats noise_const(const Params &p, const SecretKey &sk, const Ciphertext &ct, std::uint64_t m);

  // Noise of a 3-part ciphertext (after raw multiply), assuming constant message.
  NoiseStats noise3_const(const Params &p, const SecretKey &sk, const Ciphertext3 &ct3, std::uint64_t m);
}
