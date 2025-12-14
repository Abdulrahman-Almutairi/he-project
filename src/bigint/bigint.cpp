#include "bigint/bigint.h"
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace bi
{
  // remove leading zero limbs
  void BigInt::normalize()
  {
    while (!limbs.empty() && limbs.back() == 0)
      limbs.pop_back();
  }

  // convert BigInt to hexa representation
  std::string BigInt::to_hex() const
  {
    if (limbs.empty())
      return "0x0";
    std::ostringstream oss;
    oss << "0x" << std::hex << std::setfill('0');
    for (std::size_t i = limbs.size(); i-- > 0;)
    {
      if (i == limbs.size() - 1)
        oss << std::hex << limbs[i];
      else
        oss << std::setw(16) << limbs[i];
    }
    return oss.str();
  }

  // compute a = a + b
  void add_inplace(BigInt &a, const BigInt &b)
  {
    const std::size_t n = std::max(a.limbs.size(), b.limbs.size());
    if (a.limbs.size() < n)
      a.limbs.resize(n, 0);
    limb_t carry = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      limb_t ai = a.limbs[i];
      limb_t bi = i < b.limbs.size() ? b.limbs[i] : 0;
      dlimb_t s = (dlimb_t)ai + bi + carry;
      a.limbs[i] = (limb_t)s;
      carry = (limb_t)(s >> 64);
    }
    if (carry)
      a.limbs.push_back(carry);
  }

  // compute c = a + b
  BigInt add(const BigInt &a, const BigInt &b)
  {
    BigInt r = a;
    add_inplace(r, b);
    return r;
  }

  // compare a and b (limb by limb)
  int cmp(const BigInt &a, const BigInt &b)
  {
    if (a.limbs.size() != b.limbs.size())
      return (a.limbs.size() < b.limbs.size()) ? -1 : 1;
    for (std::size_t i = a.limbs.size(); i-- > 0;)
    {
      if (a.limbs[i] != b.limbs[i])
        return (a.limbs[i] < b.limbs[i]) ? -1 : 1;
    }
    return 0;
  }

  // compute a = a - b
  bool sub_inplace(BigInt &a, const BigInt &b)
  {
    // returns true if borrow occurred (i.e., a<b)
    bi::limb_t borrow = 0;
    const std::size_t n = std::max(a.limbs.size(), b.limbs.size());
    if (a.limbs.size() < n)
      a.limbs.resize(n, 0);
    for (std::size_t i = 0; i < n; ++i)
    {
      bi::limb_t ai = a.limbs[i];
      bi::limb_t bi_ = (i < b.limbs.size() ? b.limbs[i] : 0);
      bi::dlimb_t diff = (bi::dlimb_t)ai - bi_ - borrow;
      a.limbs[i] = (bi::limb_t)diff;
      // top bit of 128-bit diff is 1 if negative; convert to 0/1
      borrow = (diff >> 127) & 1;
    }
    a.normalize();
    return borrow != 0;
  }

  // read low 64 bits
  std::uint64_t to_u64(const BigInt& x){
    return x.limbs.empty() ? 0 : x.limbs[0];
  }
}
