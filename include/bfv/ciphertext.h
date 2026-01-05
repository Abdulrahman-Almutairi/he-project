#pragma once
#include "bfv/poly.h"

namespace bfv
{
  struct Ciphertext
  {
    Poly c0;
    Poly c1;

    Ciphertext() = default;
    explicit Ciphertext(const Params &p) : c0(p), c1(p) {}
  };
}
