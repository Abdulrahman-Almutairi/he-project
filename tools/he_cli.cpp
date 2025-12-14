#include <iostream>
#include <cstdint>
#include "paillier/paillier.h"
#include "bigint/bigint.h"

int main()
{
  auto kp = he::keygen(2048); // toy key

  std::uint64_t m1, m2, k;
  std::cout << "Enter m1: ";
  std::cin >> m1;
  std::cout << "Enter m2: ";
  std::cin >> m2;
  std::cout << "Enter k (for scaling m1): ";
  std::cin >> k;

  std::uint64_t n = bi::to_u64(kp.pk.n);
  m1 %= n;
  m2 %= n;

  auto c1 = he::encrypt(kp.pk, m1);
  auto c2 = he::encrypt(kp.pk, m2);

  // homomorphic add
  auto c_add = he::add(kp.pk, c1, c2);
  auto d_add = he::decrypt_u64(kp.pk, kp.sk, c_add);

  // homomorphic scale
  auto c_scale = he::scale(kp.pk, c1, k);
  auto d_scale = he::decrypt_u64(kp.pk, kp.sk, c_scale);

  std::cout << "\n=== Parameters ===\n";
  std::cout << "n                = " << n << "\n";

  std::cout << "\n=== Inputs ===\n";
  std::cout << "m1               = " << m1 << "\n";
  std::cout << "m2               = " << m2 << "\n";
  std::cout << "k                = " << k << "\n";

  std::cout << "\n=== Ciphertexts (low limb only) ===\n";
  std::cout << "Enc(m1).lo       = " << bi::to_u64(c1) << "\n";
  std::cout << "Enc(m2).lo       = " << bi::to_u64(c2) << "\n";
  std::cout << "Enc(m1+m2).lo    = " << bi::to_u64(c_add) << "\n";
  std::cout << "Enc(k*m1).lo     = " << bi::to_u64(c_scale) << "\n";

  std::cout << "\n=== Decrypted results ===\n";
  std::cout << "Dec(Enc(m1+m2))  = " << d_add
            << "  (expected " << ((m1 + m2) % n) << ")\n";
  std::cout << "Dec(Enc(k*m1))   = " << d_scale
            << "  (expected " << ((k * m1) % n) << ")\n";

  return 0;
}