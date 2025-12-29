#pragma once
#include "bigint/bigint.h"

namespace bi {

// Exact division: returns q such that x = q*d (assumes remainder is zero).
// This is a simple long-division style implementation adequate for Paillier L(u) = (u-1)/n.
BigInt div_exact(const BigInt& x, const BigInt& d);

}
