#include "bfv/bfv.h"

namespace bfv
{
  Ciphertext add(const Params &p, const Ciphertext &x, const Ciphertext &y)
  {
    Ciphertext r(p);
    for (std::size_t i = 0; i < p.N; ++i)
    {
      r.c0[i] = add_mod(x.c0[i], y.c0[i], p.q);
      r.c1[i] = add_mod(x.c1[i], y.c1[i], p.q);
    }
    return r;
  }
}
