#include "bigint/bigint.h"

namespace bi
{
  BigInt mul_schoolbook(const BigInt &a, const BigInt &b)
  {
    BigInt r;
    r.limbs.assign(a.limbs.size() + b.limbs.size(), 0);
    for (std::size_t i = 0; i < a.limbs.size(); ++i)
    {
      dlimb_t carry = 0;
      for (std::size_t j = 0; j < b.limbs.size(); ++j)
      {
        dlimb_t cur = (dlimb_t)a.limbs[i] * b.limbs[j] + r.limbs[i + j] + carry;
        r.limbs[i + j] = (limb_t)cur;
        carry = cur >> 64;
      }
      r.limbs[i + b.limbs.size()] += (limb_t)carry;
    }
    r.normalize();
    return r;
  }
}
