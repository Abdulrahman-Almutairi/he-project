#include "bfv/bfv.h"
#include "bfv/ntt.h"
#include "bfv/relin.h"
#include <cassert>
#include <iostream>

int main()
{
  bfv::Params P{1024, 12289, 2};
  std::mt19937_64 rng(123);

  auto kp = bfv::keygen_noiseless(P, rng);
  auto rlk = bfv::relin_keygen_noiseless(P, kp.sk, rng, /*base=*/4);

  for (std::uint64_t a = 0; a < P.t; ++a)
  {
    for (std::uint64_t b = 0; b < P.t; ++b)
    {
      auto ca = bfv::encrypt_const_noiseless(P, kp.pk, a, rng);
      auto cb = bfv::encrypt_const_noiseless(P, kp.pk, b, rng);

      // 3-part ciphertext from raw multiply
      auto ct3 = bfv::mul_raw(P, ca, cb);

      // Relinearize -> 2-part ciphertext
      auto ct2 = bfv::relinearize(P, ct3, rlk);

      // Compute v3 = c0 + c1*s + c2*s^2  (mod q)
      auto c1s = bfv::mul_negacyclic_ntt(P, ct3.c1, kp.sk.s);
      auto s2  = bfv::mul_negacyclic_ntt(P, kp.sk.s, kp.sk.s);
      auto c2s2 = bfv::mul_negacyclic_ntt(P, ct3.c2, s2);

      bfv::Poly v3(P);
      for (std::size_t i = 0; i < P.N; ++i)
      {
        auto t0 = bfv::add_mod(ct3.c0[i], c1s[i], P.q);
        v3[i] = bfv::add_mod(t0, c2s2[i], P.q);
      }

      // Compute v2 = c0' + c1'*s (mod q)
      auto c1ps = bfv::mul_negacyclic_ntt(P, ct2.c1, kp.sk.s);

      bfv::Poly v2(P);
      for (std::size_t i = 0; i < P.N; ++i)
      {
        v2[i] = bfv::add_mod(ct2.c0[i], c1ps[i], P.q);
      }

      // Relinearization correctness: v2 == v3 coefficient-wise
      for (std::size_t i = 0; i < P.N; ++i)
      {
        assert(v2[i] == v3[i]);
      }
    }
  }

  std::cout << "BFV relinearization (toy ring-equality) passed.\n";
  return 0;
}
