#include "bfv/bfv.h"
#include <iostream>
#include <cassert>

int main()
{
  bfv::Params P{8, 97, 17};
  std::mt19937_64 rng(123);

  auto kp = bfv::keygen(P, rng);

  // Basic sanity: polynomials have correct size and mod range
  assert(kp.sk.s.a.size() == P.N);
  assert(kp.pk.a.a.size() == P.N);
  assert(kp.pk.b.a.size() == P.N);

  for (std::size_t i = 0; i < P.N; ++i)
  {
    assert(kp.pk.a[i] < P.q);
    assert(kp.pk.b[i] < P.q);
    assert(kp.sk.s[i] < P.q);
  }

  std::cout << "BFV keygen sanity passed.\n";
  return 0;
}
