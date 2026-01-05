#pragma once
#include "bfv/poly.h"

namespace bfv
{
  // Forward negacyclic NTT (in-place)
  void ntt(Poly &a);

  // Inverse negacyclic NTT (in-place)
  void intt(Poly &a);

  // Multiply using NTT (drop-in replacement for schoolbook)
  Poly mul_negacyclic_ntt(const Poly &a, const Poly &b);
}
