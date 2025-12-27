#include "bigint/montgomery.h"
#include "bigint/mod_arith.h"
#include <algorithm>

namespace bi
{

    // Compute modular inverse of n0 mod 2^64, then negate: n0inv = -n0^{-1} mod 2^64
    // Using Newton iteration for inverse modulo 2^64.
    static std::uint64_t inv_mod_2_64(std::uint64_t a)
    {
        // a must be odd
        std::uint64_t x = 1;
        // Each iteration doubles correct bits
        for (int i = 0; i < 6; ++i)
        {
            x *= 2 - a * x;
        }
        return x;
    }

    static void ensure_k(BigInt &x, std::size_t k)
    {
        if (x.limbs.size() < k)
            x.limbs.resize(k, 0);
    }

    static BigInt mod_simple(BigInt x, const BigInt &mod)
    {
        // Temporary fallback: repeated subtraction only for small overflows.
        // With Montgomery, we should rarely use this.
        while (cmp(x, mod) >= 0)
            sub_inplace(x, mod);
        return x;
    }

    // Compute R^2 mod n by repeated doubling mod n (simple but fine at init-time).
    static BigInt compute_r2(const BigInt &n, std::size_t k)
    {
        // R = 2^(64k). We want R^2 mod n = 2^(128k) mod n.
        BigInt x = BigInt::from_u64(1);
        // Do 128*k limb shifts worth of doubling
        for (std::size_t i = 0; i < 128 * k; ++i)
        {
            // x = (x * 2) mod n
            x = mod_simple(add(x, x), n); 
        }
        return x;
    }

    // If you don't already have add_mod, implement a tiny one locally:
    static BigInt add_mod_local(BigInt a, const BigInt &b, const BigInt &mod)
    {
        // a += b
        std::size_t nlimbs = std::max(a.limbs.size(), b.limbs.size());
        a.limbs.resize(nlimbs, 0);
        std::uint64_t carry = 0;
        for (std::size_t i = 0; i < nlimbs; ++i)
        {
            unsigned __int128 sum = (unsigned __int128)a.limbs[i] + (unsigned __int128)(i < b.limbs.size() ? b.limbs[i] : 0) + carry;
            a.limbs[i] = (std::uint64_t)sum;
            carry = (std::uint64_t)(sum >> 64);
        }
        if (carry)
            a.limbs.push_back(carry);
        a.normalize();
        if (cmp(a, mod) >= 0)
            sub_inplace(a, mod);
        return a;
    }

    MontCtx mont_ctx(const BigInt &mod)
    {
        MontCtx ctx;
        ctx.n = mod;
        ctx.n.normalize();
        ctx.k = ctx.n.limbs.size();
        // modulus must be odd
        std::uint64_t n0 = ctx.n.limbs.empty() ? 0 : ctx.n.limbs[0];
        std::uint64_t inv = inv_mod_2_64(n0);
        ctx.n0inv = (std::uint64_t)(0 - inv); // -n0^{-1} mod 2^64

        // Compute r2 = R^2 mod n (init cost)
        // NOTE: compute_r2 uses add_mod; we use local helper:
        BigInt x = BigInt::from_u64(1);
        for (std::size_t i = 0; i < 128 * ctx.k; ++i)
        {
            x = add_mod_local(x, x, ctx.n);
        }
        ctx.r2 = x;
        return ctx;
    }

    BigInt mont_mul(const MontCtx &ctx, const BigInt &a_in, const BigInt &b_in)
    {
        // CIOS Montgomery multiplication.
        // Input a,b must be < n and typically in Montgomery form.
        const std::size_t k = ctx.k;

        BigInt a = a_in;
        a.normalize();
        BigInt b = b_in;
        b.normalize();

        ensure_k(a, k);
        ensure_k(b, k);

        // t has k+1 limbs (accumulator)
        std::vector<std::uint64_t> t(k + 1, 0);

        for (std::size_t i = 0; i < k; ++i)
        {
            // t += a_i * b
            unsigned __int128 carry = 0;
            for (std::size_t j = 0; j < k; ++j)
            {
                unsigned __int128 prod = (unsigned __int128)a.limbs[i] * (unsigned __int128)b.limbs[j];
                unsigned __int128 sum = (unsigned __int128)t[j] + prod + carry;
                t[j] = (std::uint64_t)sum;
                carry = (sum >> 64);
            }
            // add carry to t[k]
            unsigned __int128 sumk = (unsigned __int128)t[k] + carry;
            t[k] = (std::uint64_t)sumk;

            // m = t0 * n0inv mod 2^64
            std::uint64_t m = (std::uint64_t)((unsigned __int128)t[0] * (unsigned __int128)ctx.n0inv);

            // t += m * n
            carry = 0;
            for (std::size_t j = 0; j < k; ++j)
            {
                unsigned __int128 prod = (unsigned __int128)m * (unsigned __int128)ctx.n.limbs[j];
                unsigned __int128 sum = (unsigned __int128)t[j] + prod + carry;
                t[j] = (std::uint64_t)sum;
                carry = (sum >> 64);
            }
            sumk = (unsigned __int128)t[k] + carry;
            t[k] = (std::uint64_t)sumk;

            // shift t right by one limb: drop t[0]
            for (std::size_t j = 0; j < k; ++j)
                t[j] = t[j + 1];
            t[k] = 0;
        }

        BigInt res;
        res.limbs.assign(t.begin(), t.begin() + k);
        res.normalize();

        // final conditional subtraction
        if (cmp(res, ctx.n) >= 0)
            sub_inplace(res, ctx.n);
        return res;
    }

    BigInt to_mont(const MontCtx &ctx, const BigInt &x)
    {
        // x * R mod n = mont_mul(x, R^2) because mont_mul does *R^{-1}
        return mont_mul(ctx, x, ctx.r2);
    }

    BigInt from_mont(const MontCtx &ctx, const BigInt &x)
    {
        // x * 1 mod n = mont_mul(x, 1)
        return mont_mul(ctx, x, BigInt::from_u64(1));
    }

    BigInt mont_pow_u64(const MontCtx &ctx, BigInt base, std::uint64_t exp)
    {
        BigInt a = to_mont(ctx, base);
        BigInt one = to_mont(ctx, BigInt::from_u64(1));
        BigInt acc = one;

        while (exp > 0)
        {
            if (exp & 1)
                acc = mont_mul(ctx, acc, a);
            a = mont_mul(ctx, a, a);
            exp >>= 1;
        }
        return from_mont(ctx, acc);
    }

} // namespace bi
