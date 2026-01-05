#include "bfv/bfv.h"
#include <cassert>
#include <iostream>

int main()
{
  bfv::Params P{8, 97, 2};
  std::mt19937_64 rng(123);

  auto kp = bfv::keygen(P, rng);

  for (std::uint64_t m = 0; m < P.t; ++m)
  {
    auto ct = bfv::encrypt_const(P, kp.pk, m, rng);
    auto dec = bfv::decrypt_const(P, kp.sk, ct);
    assert(dec == m);
  }

  // randomness check: same message twice => ciphertext differs (likely)
  auto ct1 = bfv::encrypt_const(P, kp.pk, 5, rng);
  auto ct2 = bfv::encrypt_const(P, kp.pk, 5, rng);
  // Compare a coefficient or two (full compare is fine but verbose)
  bool different = (ct1.c0[0] != ct2.c0[0]) || (ct1.c1[0] != ct2.c1[0]);
  assert(different);

  std::cout << "BFV encrypt/decrypt (toy const) passed.\n";
  return 0;
}
