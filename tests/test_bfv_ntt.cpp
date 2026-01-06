#include "bfv/params.h"
#include "bfv/poly.h"
#include "bfv/ntt.h"
#include <cassert>
#include <iostream>

int main()
{
  bfv::Params P{1024, 12289, 2};
  bfv::Poly a(P), b(P);

  a[0] = 1; a[1] = 2; a[2] = 3;
  b[0] = 4; b[1] = 5; b[2] = 6;

  auto c1 = bfv::mul_negacyclic_schoolbook(a, b);
  auto c2 = bfv::mul_negacyclic_ntt(P, a, b);

  for (std::size_t i = 0; i < P.N; ++i)
    assert(c1[i] == c2[i]);

  std::cout << "BFV NTT multiplication passed.\n";
  return 0;
}
