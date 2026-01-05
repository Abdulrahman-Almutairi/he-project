#pragma once
#include "bfv/poly.h"

namespace bfv
{
  struct Ciphertext3
  {
    Poly c0;
    Poly c1;
    Poly c2;

    Ciphertext3() = default;
    explicit Ciphertext3(const Params &p) : c0(p), c1(p), c2(p) {}
  };
}
