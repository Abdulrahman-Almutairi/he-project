#include "bigint/gcd.h"
#include <cstdint>

namespace bi {

static inline bool is_zero(const BigInt& x) { return x.limbs.empty(); }
static inline bool is_even(const BigInt& x) { return x.limbs.empty() || ((x.limbs[0] & 1ULL) == 0); }

static void shr1_inplace(BigInt& x) {
  std::uint64_t carry = 0;
  for (std::size_t i = x.limbs.size(); i-- > 0;) {
    std::uint64_t new_carry = x.limbs[i] & 1ULL;
    x.limbs[i] = (x.limbs[i] >> 1) | (carry << 63);
    carry = new_carry;
  }
  x.normalize();
}

static void shl1_inplace(BigInt& x) {
  std::uint64_t carry = 0;
  for (std::size_t i = 0; i < x.limbs.size(); ++i) {
    unsigned __int128 s = ((unsigned __int128)x.limbs[i] << 1) | carry;
    x.limbs[i] = (std::uint64_t)s;
    carry = (std::uint64_t)(s >> 64);
  }
  if (carry) x.limbs.push_back(carry);
  x.normalize();
}

// gcd via Stein's algorithm
BigInt gcd(BigInt a, BigInt b) {
  a.normalize();
  b.normalize();

  if (is_zero(a)) return b;
  if (is_zero(b)) return a;

  // count common factors of 2
  std::size_t shift = 0;
  while (is_even(a) && is_even(b)) {
    shr1_inplace(a);
    shr1_inplace(b);
    ++shift;
  }

  // make a odd
  while (is_even(a)) shr1_inplace(a);

  do {
    while (is_even(b)) shr1_inplace(b);

    // now a,b odd; ensure a <= b
    if (cmp(a, b) > 0) {
      BigInt tmp = a;
      a = b;
      b = tmp;
    }

    // b = b - a (non-negative)
    sub_inplace(b, a);
    b.normalize();

  } while (!is_zero(b));

  // restore factor 2^shift
  while (shift--) shl1_inplace(a);
  a.normalize();
  return a;
}

} // namespace bi
