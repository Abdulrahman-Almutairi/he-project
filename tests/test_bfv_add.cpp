#include "bfv/bfv.h"
#include <cassert>
#include <iostream>

int main()
{
  bfv::Params P{8, 97, 2};
  std::mt19937_64 rng(123);

  auto kp = bfv::keygen(P, rng);

  // test all combos in Z_t (t=2)
  for (std::uint64_t a = 0; a < P.t; ++a)
  {
    for (std::uint64_t b = 0; b < P.t; ++b)
    {
      auto ca = bfv::encrypt_const(P, kp.pk, a, rng);
      auto cb = bfv::encrypt_const(P, kp.pk, b, rng);

      auto cs = bfv::add(P, ca, cb);
      auto dec = bfv::decrypt_const(P, kp.sk, cs);

      assert(dec == ((a + b) % P.t));
    }
  }

  std::cout << "BFV homomorphic add passed.\n";
  return 0;
}
