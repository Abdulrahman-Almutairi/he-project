#include "bfv/bfv.h"
#include "bfv/relin.h"
#include <cassert>
#include <iostream>

int main()
{
  bfv::Params P{1024, 12289, 2};
  std::mt19937_64 rng(123);

  auto kp = bfv::keygen_noiseless(P, rng); // includes kp.rlk
  auto rlk = bfv::relin_keygen_noiseless(P, kp.sk, rng, 4);

  for (std::uint64_t a = 0; a < P.t; ++a)
  {
    for (std::uint64_t b = 0; b < P.t; ++b)
    {
      auto ca = bfv::encrypt_const_noiseless(P, kp.pk, a, rng);
      auto cb = bfv::encrypt_const_noiseless(P, kp.pk, b, rng);

      auto cc = bfv::mul(P, ca, cb, rlk);

      // Check ring-equality against ct3 decrypt (same approach as relin test)
      // v2 = c0 + c1*s  (mod q)
      auto dec = bfv::decrypt_const(P, kp.sk, cc);

      // With toy params, decrypt_const is reliable for constant messages even after one mul_relin.
      // Expected: a*b mod t
      assert(dec == ((a * b) % P.t));
    }
  }

  std::cout << "BFV mul (with relinearization) passed.\n";
  return 0;
}
