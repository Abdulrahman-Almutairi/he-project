#pragma once 
#include "bigint/bigint.h"

namespace bi {
    // Fast Modular multiplication (a * b) % mod
    std::uint64_t mul_mod_u64(std::uint64_t a, std::uint64_t b, std::uint64_t mod);

    // Fast Modular exponent (base ^ exp) % mod
    std::uint64_t pow_mod_u64(std::uint64_t base, std::uint64_t exp, std::uint64_t mod);
    
    // gcd
    std::uint64_t gcd_u64(std::uint64_t a, std::uint64_t b);

    // lcm
    std::uint64_t lcm_u64(std::uint64_t a, std::uint64_t b);
    
    // mod inverse (return false if no inverse)
    bool modinv_u64(std::uint64_t a, std::uint64_t mod, std::uint64_t& inv_out);
    
    // Modular multiplication (a * b) % mod -- naive
    BigInt mul_mod(const BigInt& a, const BigInt& b, const BigInt& mod);
    // Modular exponent (base ^ exp) % mod  -- naive
    BigInt pow_mod(BigInt base, std::uint64_t exp, const BigInt& mod);

    // big exponent
    bi::BigInt pow_mod_bigexp(BigInt base, const BigInt& exp, const BigInt& mod);

}