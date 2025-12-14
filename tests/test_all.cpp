#include <cassert>
#include <iostream>
#include "bigint/bigint.h"
#include "bigint/mod_arith.h"
#include "paillier/paillier.h"
#include <vector>

void test_homomorphic_add_scale()
{
  auto kp = he::keygen(2048);

  std::uint64_t n = bi::to_u64(kp.pk.n);

  std::uint64_t m1 = 5;
  std::uint64_t m2 = 7;
  std::uint64_t k = 3;

  m1 %= n;
  m2 %= n;

  auto c1 = he::encrypt(kp.pk, m1);
  auto c2 = he::encrypt(kp.pk, m2);

  // Enc(m1 + m2)
  auto c_add = he::add(kp.pk, c1, c2);
  auto d_add = he::decrypt_u64(kp.pk, kp.sk, c_add);
  std::uint64_t expected_add = (m1 + m2) % n;
  assert(d_add == expected_add);

  // Enc(k * m1)
  auto c_scale = he::scale(kp.pk, c1, k);
  auto d_scale = he::decrypt_u64(kp.pk, kp.sk, c_scale);
  std::uint64_t expected_scale = (k * m1) % n;
  assert(d_scale == expected_scale);

  std::cout << "Homomorphic add/scale tests passed.\n";
}

void test_paillier_roundtrip()
{
  auto kp = he::keygen(2048);

  std::vector<std::uint64_t> msgs = {0, 1, 2, 7, 42, 123, 123456789 % bi::to_u64(kp.pk.n)};
  for (auto m : msgs)
  {
    auto c = he::encrypt(kp.pk, m);
    auto d = he::decrypt_u64(kp.pk, kp.sk, c);
    assert(d == (m % bi::to_u64(kp.pk.n)));
  }
  std::cout << "Paillier toy round-trip passed.\n";
}

void test_pow_mod()
{
  using namespace bi;
  BigInt base = BigInt::from_u64(8);
  BigInt mod = BigInt::from_u64(22);

  auto r = pow_mod(base, 3, mod); // 8^3 mod 22 = 6
  assert(r.limbs[0] == 6);
  std::cout << "pow_mod test passed " << r.limbs[0] << std::endl; // expect 6
}

int main()
{
  using namespace bi;
  BigInt a = BigInt::from_u64(0xffffffffffffffffULL);
  BigInt b = BigInt::from_u64(2);
  auto s = add(a, b);
  assert(s.limbs.size() == 2); // carry out
  auto m = mul_schoolbook(a, b);
  assert(m.limbs.size() >= 2);

  test_pow_mod();
  test_paillier_roundtrip();
  test_homomorphic_add_scale();
  std::cout << "All tests passed.\n";

  // std::cout << "lol.\n";

  return 0;
}
