#include "bigint/inv.h"
#include <algorithm>

namespace bi {

// --- local helpers (so we don't depend on other files) ---

static bool is_zero(const BigInt& x) { return x.limbs.empty(); }
static bool is_one (const BigInt& x) { return x.limbs.size()==1 && x.limbs[0]==1; }
static bool is_even(const BigInt& x) { return x.limbs.empty() ? true : ((x.limbs[0] & 1ULL) == 0); }

static BigInt add_big(BigInt a, const BigInt& b) {
  std::size_t n = std::max(a.limbs.size(), b.limbs.size());
  a.limbs.resize(n, 0);
  std::uint64_t carry = 0;
  for (std::size_t i = 0; i < n; ++i) {
    unsigned __int128 sum = (unsigned __int128)a.limbs[i]
                          + (unsigned __int128)(i < b.limbs.size() ? b.limbs[i] : 0)
                          + carry;
    a.limbs[i] = (std::uint64_t)sum;
    carry = (std::uint64_t)(sum >> 64);
  }
  if (carry) a.limbs.push_back(carry);
  a.normalize();
  return a;
}

static BigInt sub_big(BigInt a, const BigInt& b) {
  // assumes a >= b
  sub_inplace(a, b);
  a.normalize();
  return a;
}

static void shr1_inplace(BigInt& x) {
  std::uint64_t carry = 0;
  for (std::size_t i = x.limbs.size(); i-- > 0;) {
    std::uint64_t new_carry = x.limbs[i] & 1ULL;
    x.limbs[i] = (x.limbs[i] >> 1) | (carry << 63);
    carry = new_carry;
  }
  x.normalize();
}

static BigInt div2_mod(BigInt x, const BigInt& mod) {
  // compute x/2 mod mod, assuming mod is odd.
  if (is_even(x)) {
    shr1_inplace(x);
    return x;
  }
  // if x is odd, (x+mod) is even, then divide by 2.
  x = add_big(x, mod);
  shr1_inplace(x);
  // keep it reduced
  if (cmp(x, mod) >= 0) sub_inplace(x, mod);
  x.normalize();
  return x;
}

// Keep value in [0, mod)
static void reduce_once(BigInt& x, const BigInt& mod) {
  x.normalize();
  if (cmp(x, mod) >= 0) sub_inplace(x, mod);
  x.normalize();
}

BigInt inv_mod_odd(const BigInt& a_in, const BigInt& mod_in) {
  BigInt mod = mod_in; mod.normalize();
  BigInt a   = a_in;   a.normalize();

  // mod must be odd
  // (We assume gcd(a,mod)=1 for Paillier)
  if (mod.limbs.empty() || is_even(mod)) return BigInt::from_u64(0);

  // Reduce a modulo mod (so 0 <= a < mod)
  while (cmp(a, mod) >= 0) sub_inplace(a, mod);
  a.normalize();

  // Binary extended gcd inverse algorithm (HAC style):
  // u=a, v=mod, x1=1, x2=0
  BigInt u = a;
  BigInt v = mod;
  BigInt x1 = BigInt::from_u64(1);
  BigInt x2 = BigInt::from_u64(0);

  while (!is_one(u) && !is_one(v)) {
    while (is_even(u)) {
      shr1_inplace(u);
      x1 = div2_mod(x1, mod);
    }
    while (is_even(v)) {
      shr1_inplace(v);
      x2 = div2_mod(x2, mod);
    }

    if (cmp(u, v) >= 0) {
      // u = u - v
      sub_inplace(u, v);
      u.normalize();

      // x1 = x1 - x2 (mod mod)
      if (cmp(x1, x2) >= 0) {
        x1 = sub_big(x1, x2);
      } else {
        // x1 = x1 + mod - x2
        x1 = add_big(x1, mod);
        x1 = sub_big(x1, x2);
      }
      reduce_once(x1, mod);
    } else {
      // v = v - u
      sub_inplace(v, u);
      v.normalize();

      // x2 = x2 - x1 (mod mod)
      if (cmp(x2, x1) >= 0) {
        x2 = sub_big(x2, x1);
      } else {
        x2 = add_big(x2, mod);
        x2 = sub_big(x2, x1);
      }
      reduce_once(x2, mod);
    }
  }

  BigInt res = is_one(u) ? x1 : x2;
  // ensure reduced
  while (cmp(res, mod) >= 0) sub_inplace(res, mod);
  res.normalize();
  return res;
}

} // namespace bi
