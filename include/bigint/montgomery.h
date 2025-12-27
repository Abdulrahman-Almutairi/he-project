#pragma once
#include "bigint/bigint.h"
#include <cstddef>
#include <cstdint>

namespace bi
{

    struct MontCtx
    {
        BigInt n;                // modulus
        std::size_t k = 0;       // number of limbs
        std::uint64_t n0inv = 0; // -n^{-1} mod 2^64
        BigInt r2;               // R^2 mod n, where R = 2^(64k)
    };

    MontCtx mont_ctx(const BigInt &mod);              // build context
    BigInt mont_reduce(const MontCtx &ctx, BigInt t); // optional (not required if using CIOS)
    BigInt mont_mul(const MontCtx &ctx, const BigInt &a, const BigInt &b);
    BigInt mont_pow_u64(const MontCtx &ctx, BigInt base, std::uint64_t exp);
    BigInt to_mont(const MontCtx &ctx, const BigInt &x);   // x*R mod n
    BigInt from_mont(const MontCtx &ctx, const BigInt &x); // x mod n (x in mont form)

} // namespace bi
