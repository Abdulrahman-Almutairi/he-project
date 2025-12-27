#include <cassert>
#include <iostream>
#include "bigint/bigint.h"
#include "bigint/mod_arith.h"
#include "paillier/paillier.h"
#include <vector>

void test_pow_mod_bigexp_matches_u64() {
  using namespace bi;

  BigInt mod  = BigInt::from_u64(65537);     // small odd modulus
  BigInt base = BigInt::from_u64(12345);

  // pick a small exponent we can represent as both u64 and BigInt
  std::uint64_t e_u64 = 1234;
  BigInt e_big = BigInt::from_u64(e_u64);

  auto a = pow_mod(base, e_u64, mod);
  auto b = pow_mod_bigexp(base, e_big, mod);

  assert(to_u64(a) == to_u64(b));
  std::cout << "pow_mod_bigexp matches pow_mod for small exponent.\n";
}

void test_enc_gemm_cp() {
  auto kp = he::keygen(2048);
  std::uint64_t n = bi::to_u64(kp.pk.n);

  // A (2x3):
  // [1 2 3
  //  4 5 6]
  std::size_t A_rows = 2, A_cols = 3;
  std::vector<std::uint64_t> A_plain = {1,2,3, 4,5,6};

  std::vector<bi::BigInt> enc_A;
  enc_A.reserve(A_plain.size());
  for (auto v : A_plain) enc_A.push_back(he::encrypt(kp.pk, v % n));

  // B (3x2):
  // [7  8
  //  9 10
  // 11 12]
  std::size_t B_cols = 2;
  std::vector<std::uint64_t> B = {7,8, 9,10, 11,12}; // row-major (3x2)

  // Expected C = A*B (2x2):
  // C00 = 1*7 + 2*9 + 3*11 = 58
  // C01 = 1*8 + 2*10+ 3*12 = 64
  // C10 = 4*7 + 5*9 + 6*11 = 139
  // C11 = 4*8 + 5*10+ 6*12 = 154
  auto enc_C = he::enc_gemm_cp(kp.pk, enc_A, A_rows, A_cols, B, B_cols);
  assert(enc_C.size() == A_rows * B_cols);

  auto c00 = he::decrypt_u64(kp.pk, kp.sk, enc_C[0]);
  auto c01 = he::decrypt_u64(kp.pk, kp.sk, enc_C[1]);
  auto c10 = he::decrypt_u64(kp.pk, kp.sk, enc_C[2]);
  auto c11 = he::decrypt_u64(kp.pk, kp.sk, enc_C[3]);

  assert(c00 == (58  % n));
  assert(c01 == (64  % n));
  assert(c10 == (139 % n));
  assert(c11 == (154 % n));

  std::cout << "Encrypted GEMM (c×p) passed.\n";
}


void test_enc_gemv_cp() {
  auto kp = he::keygen(2048);
  std::uint64_t n = bi::to_u64(kp.pk.n);

  // A = [[1,2,3],
  //      [4,5,6]]
  std::size_t rows = 2, cols = 3;
  std::vector<std::uint64_t> A_plain = {1,2,3, 4,5,6};
  std::vector<bi::BigInt> enc_A;
  enc_A.reserve(A_plain.size());
  for (auto v : A_plain) enc_A.push_back(he::encrypt(kp.pk, v % n));

  // x = [7,8,9]
  std::vector<std::uint64_t> x = {7,8,9};

  // y = A*x = [1*7+2*8+3*9, 4*7+5*8+6*9] = [50, 122]
  auto enc_y = he::enc_gemv_cp(kp.pk, enc_A, rows, cols, x);
  assert(enc_y.size() == rows);

  auto y0 = he::decrypt_u64(kp.pk, kp.sk, enc_y[0]);
  auto y1 = he::decrypt_u64(kp.pk, kp.sk, enc_y[1]);

  assert(y0 == (50 % n));
  assert(y1 == (122 % n));

  std::cout << "Encrypted GEMV (c×p) passed.\n";
}


void test_enc_dot_cp()
{
  auto kp = he::keygen(2048);
  std::uint64_t n = bi::to_u64(kp.pk.n);

  // plaintext vectors
  std::vector<std::uint64_t> a = {2, 3, 4};
  std::vector<std::uint64_t> b = {5, 6, 7};

  // encrypt a
  std::vector<bi::BigInt> enc_a;
  for (auto x : a)
    enc_a.push_back(he::encrypt(kp.pk, x % n));

  // compute Enc(dot)
  auto c = he::enc_dot_cp(kp.pk, enc_a, b);

  // decrypt result
  auto d = he::decrypt_u64(kp.pk, kp.sk, c);

  // expected dot
  std::uint64_t expected = (2 * 5 + 3 * 6 + 4 * 7) % n; // 10+18+28=56
  assert(d == expected);

  std::cout << "Encrypted dot-product (c×p) passed.\n";
}

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
  test_enc_dot_cp();
  test_enc_gemv_cp();
  test_enc_gemm_cp();
  test_pow_mod_bigexp_matches_u64();
  std::cout << "All tests passed.\n";

  // std::cout << "lol.\n";

  return 0;
}
