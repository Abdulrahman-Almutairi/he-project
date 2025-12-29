#include "bigint/mod_arith.h"
#include "bigint/montgomery.h"
#include <stdexcept>

namespace bi
{
    // static int cmp(const BigInt &a, const BigInt &b)
    // {
    //     if (a.limbs.size() != b.limbs.size())
    //         return (a.limbs.size() < b.limbs.size()) ? -1 : 1;
    //     for (std::size_t i = a.limbs.size(); i-- > 0;)
    //     {
    //         if (a.limbs[i] != b.limbs[i])
    //             return (a.limbs[i] < b.limbs[i]) ? -1 : 1;
    //     }
    //     return 0;
    // }

    static const MontCtx &mont_ctx_cached(const BigInt &mod)
    {
        static thread_local bool inited = false;
        static thread_local MontCtx ctx;
        static thread_local BigInt last_mod;

        BigInt m = mod;
        m.normalize();

        if (!inited || cmp(last_mod, m) != 0)
        {
            ctx = mont_ctx(m);
            last_mod = m;
            inited = true;
        }
        return ctx;
    }

    // extern bool sub_inplace(BigInt &a, const BigInt &b);

    // naive modular reduction by subtraction
    static void mod_reduce(BigInt &x, const BigInt &mod)
    {
        // keep subtracting mod until x < mod
        // FIXME: replace with Montgomery
        throw std::runtime_error("mod_reduce called: Montgomery path not used");
    }

    // 64-bit fast modular multiplication using 128-bit intermediate
    std::uint64_t mul_mod_u64(std::uint64_t a, std::uint64_t b, std::uint64_t mod)
    {
        // (a * b) % mod computed via 128-bit product then remainder
        unsigned __int128 prod = (unsigned __int128)a * (unsigned __int128)b;
        return (std::uint64_t)(prod % mod);
    }

    std::uint64_t pow_mod_u64(std::uint64_t base, std::uint64_t exp, std::uint64_t mod)
    {
        std::uint64_t res = 1 % mod;
        std::uint64_t x = base % mod;

        while (exp > 0)
        {
            if (exp & 1)
                res = mul_mod_u64(res, x, mod);
            x = mul_mod_u64(x, x, mod);
            exp >>= 1;
        }
        return res;
    }

    BigInt mul_mod(const BigInt &a, const BigInt &b, const BigInt &mod)
    {
        // fast path if everything fits in one limb
        if (a.limbs.size() <= 1 && b.limbs.size() <= 1 && mod.limbs.size() == 1)
        {
            std::uint64_t av = a.limbs.empty() ? 0 : a.limbs[0];
            std::uint64_t bv = b.limbs.empty() ? 0 : b.limbs[0];
            std::uint64_t mv = mod.limbs[0];
            BigInt r = BigInt::from_u64(mul_mod_u64(av, bv, mv));
            return r;
        }

        // Montgomery path if modulus is odd
        if (!mod.limbs.empty() && (mod.limbs[0] & 1ULL))
        {
            const MontCtx &ctx = mont_ctx_cached(mod);

            auto A = to_mont(ctx, a);
            auto B = to_mont(ctx, b);
            auto C = mont_mul(ctx, A, B);
            return from_mont(ctx, C);
        }

        // generic (slow) fallback
        BigInt prod = mul_schoolbook(a, b);
        // static void mod_reduce(BigInt&, const BigInt&); // forward if needed
        mod_reduce(prod, mod);
        return prod;
    }

    BigInt pow_mod(BigInt base, std::uint64_t exp, const BigInt &mod)
    {
        // fast path if base and mod fit in one limb
        if (base.limbs.size() <= 1 && mod.limbs.size() == 1)
        {
            std::uint64_t bv = base.limbs.empty() ? 0 : base.limbs[0];
            std::uint64_t mv = mod.limbs[0];
            BigInt r = BigInt::from_u64(pow_mod_u64(bv, exp, mv));
            return r;
        }

        if (!mod.limbs.empty() && (mod.limbs[0] & 1ULL))
        {
            const MontCtx &ctx = mont_ctx_cached(mod);
            return mont_pow_u64(ctx, base, exp);
        }

        // generic (slow) fallback
        BigInt result = BigInt::from_u64(1);

        while (exp > 0)
        {
            if (exp & 1)
            {
                result = mul_mod(result, base, mod);
            }
            base = mul_mod(base, base, mod);
            exp >>= 1;
        }

        return result;
    }

    std::uint64_t gcd_u64(std::uint64_t a, std::uint64_t b)
    {
        while (b)
        {
            auto t = a % b;
            a = b;
            b = t;
        }
        return a;
    }

    std::uint64_t lcm_u64(std::uint64_t a, std::uint64_t b)
    {
        if (a == 0 || b == 0)
        {
            return 0;
        }
        return (a / gcd_u64(a, b)) * b;
    }

    static std::int64_t egcd_i64(std::int64_t a, std::int64_t b, std::int64_t &x, std::int64_t &y)
    {
        if (b == 0)
        {
            x = 1;
            y = 0;
            return a;
        }
        std::int64_t x1, y1;
        auto g = egcd_i64(b, a % b, x1, y1);
        x = y1;
        y = x1 - (a / b) * y1;
        return g;
    }

    bool modinv_u64(std::uint64_t a, std::uint64_t mod, std::uint64_t &inv_out)
    {
        if (mod == 0)
            return false;
        std::int64_t x, y;
        auto g = egcd_i64((std::int64_t)(a % mod), (std::int64_t)mod, x, y);
        if (g != 1 && g != -1)
            return false;
        std::int64_t inv = x % (std::int64_t)mod;
        inv_out = (std::uint64_t)inv;
        return true;
    }

    static bool get_bit(const bi::BigInt &e, std::size_t bit)
    {
        std::size_t limb = bit / 64;
        std::size_t off = bit % 64;
        if (limb >= e.limbs.size())
            return false;
        return (e.limbs[limb] >> off) & 1ULL;
    }

    static std::size_t bitlen(const bi::BigInt &e)
    {
        if (e.limbs.empty())
            return 0;
        std::size_t i = e.limbs.size() - 1;
        std::uint64_t top = e.limbs[i];
        std::size_t bits = i * 64;
        while (top)
        {
            top >>= 1;
            bits++;
        }
        return bits;
    }

    bi::BigInt pow_mod_bigexp(bi::BigInt base, const bi::BigInt &exp, const bi::BigInt &mod)
    {
        bi::BigInt result = bi::BigInt::from_u64(1);

        std::size_t bl = bitlen(exp);
        for (std::size_t i = 0; i < bl; ++i)
        {
            if (get_bit(exp, i))
                result = bi::mul_mod(result, base, mod);
            base = bi::mul_mod(base, base, mod);
        }
        return result;
    }
}