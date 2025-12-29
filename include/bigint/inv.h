#pragma once
#include "bigint/bigint.h"

namespace bi {

// Modular inverse of a modulo mod (mod must be odd, gcd(a,mod)=1).
// Returns x such that (a*x) % mod == 1.
BigInt inv_mod_odd(const BigInt& a, const BigInt& mod);

} // namespace bi
