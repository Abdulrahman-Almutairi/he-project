#include "bfv/bfv.h"
#include "bfv/ntt.h"
#include <cassert>
#include <iostream>
#include <random>

int main()
{
  bfv::Params P{1024, 12289, 2};
  std::mt19937_64 rng(123);

  // random poly
  bfv::Poly a(P);
  for (std::size_t i = 0; i < P.N; ++i)
    a[i] = rng() % P.q;

  bfv::Poly orig = a;

  bfv::ntt(P, a);
  bfv::intt(P, a);

  for (std::size_t i = 0; i < P.N; ++i)
    assert(a[i] == orig[i]);

  std::cout << "BFV NTT roundtrip passed for N=1024, q=12289.\n";
  return 0;
}
