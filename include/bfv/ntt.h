#pragma once
#include "bfv/bfv.h"

namespace bfv
{
  // Negacyclic NTT in R_q = Z_q[X]/(X^N + 1)
  void ntt(const Params &p, Poly &a);
  void intt(const Params &p, Poly &a);

  Poly mul_negacyclic_ntt(const Params &p, const Poly &a, const Poly &b);
}
