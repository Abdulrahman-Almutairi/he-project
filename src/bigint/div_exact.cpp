#include "bigint/div_exact.h"
#include <stdexcept>
#include <vector>
#include <cstdint>
#include <algorithm>

namespace bi {

static inline int clz64(std::uint64_t x) {
#if defined(__GNUG__)
    return x ? __builtin_clzll(x) : 64;
#else
    // fallback
    int n = 0;
    while (x && (x >> 63) == 0) { x <<= 1; ++n; }
    return x ? n : 64;
#endif
}

// shift left by s bits (0<=s<64)
static BigInt shl_bits(const BigInt& a, unsigned s) {
    if (s == 0) return a;
    BigInt out;
    out.limbs.resize(a.limbs.size());
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < a.limbs.size(); ++i) {
        std::uint64_t x = a.limbs[i];
        out.limbs[i] = (x << s) | carry;
        carry = (s == 0) ? 0 : (x >> (64 - s));
    }
    if (carry) out.limbs.push_back(carry);
    out.normalize();
    return out;
}

// shift right by s bits (0<=s<64)
static BigInt shr_bits(const BigInt& a, unsigned s) {
    if (s == 0) return a;
    BigInt out = a;
    std::uint64_t carry = 0;
    for (std::size_t i = out.limbs.size(); i-- > 0;) {
        std::uint64_t x = out.limbs[i];
        out.limbs[i] = (x >> s) | (carry << (64 - s));
        carry = x & ((1ULL << s) - 1ULL);
    }
    out.normalize();
    return out;
}

// u -= v*qhat at offset j, returns true if underflow happened
static bool submul_at(std::vector<std::uint64_t>& u, const std::vector<std::uint64_t>& v,
                      std::size_t j, std::uint64_t qhat)
{
    unsigned __int128 borrow = 0;
    unsigned __int128 carry = 0;

    const std::size_t n = v.size();
    for (std::size_t i = 0; i < n; ++i) {
        unsigned __int128 prod = (unsigned __int128)qhat * (unsigned __int128)v[i] + carry;
        std::uint64_t p_lo = (std::uint64_t)prod;
        carry = prod >> 64;

        unsigned __int128 sub = (unsigned __int128)u[j + i] - (unsigned __int128)p_lo - borrow;
        u[j + i] = (std::uint64_t)sub;
        borrow = (sub >> 127); // 1 if underflow (wrap)
    }

    unsigned __int128 sub = (unsigned __int128)u[j + n] - (unsigned __int128)carry - borrow;
    u[j + n] = (std::uint64_t)sub;
    return (sub >> 127) != 0; // underflow?
}

static void add_at(std::vector<std::uint64_t>& u, const std::vector<std::uint64_t>& v, std::size_t j)
{
    unsigned __int128 carry = 0;
    const std::size_t n = v.size();
    for (std::size_t i = 0; i < n; ++i) {
        unsigned __int128 sum = (unsigned __int128)u[j + i] + (unsigned __int128)v[i] + carry;
        u[j + i] = (std::uint64_t)sum;
        carry = sum >> 64;
    }
    unsigned __int128 sum = (unsigned __int128)u[j + n] + carry;
    u[j + n] = (std::uint64_t)sum;
}

// Returns q = x / d, assuming exact division
BigInt div_exact(const BigInt& x_in, const BigInt& d_in)
{
    BigInt x = x_in; x.normalize();
    BigInt d = d_in; d.normalize();

    if (d.limbs.empty()) throw std::runtime_error("div_exact: divide by zero");
    if (x.limbs.empty()) return BigInt::from_u64(0);
    if (cmp(x, d) < 0) throw std::runtime_error("div_exact: x < d but expected exact division");

    // Single-limb divisor (fast)
    if (d.limbs.size() == 1) {
        std::uint64_t dv = d.limbs[0];
        if (dv == 0) throw std::runtime_error("div_exact: divide by zero");
        BigInt q;
        q.limbs.resize(x.limbs.size());
        unsigned __int128 rem = 0;
        for (std::size_t i = x.limbs.size(); i-- > 0;) {
            unsigned __int128 cur = (rem << 64) | x.limbs[i];
            std::uint64_t qi = (std::uint64_t)(cur / dv);
            rem = cur % dv;
            q.limbs[i] = qi;
        }
        q.normalize();
        if (rem != 0) throw std::runtime_error("div_exact: non-zero remainder");
        return q;
    }

    // Knuth division D, base B=2^64
    std::vector<std::uint64_t> v = d.limbs;
    std::vector<std::uint64_t> u = x.limbs;

    const std::size_t n = v.size();
    std::size_t m = (u.size() >= n) ? (u.size() - n) : 0;

    // Normalize so top limb of v has highest bit set
    unsigned s = (unsigned)clz64(v.back());
    BigInt vn; vn.limbs = v; vn.normalize();
    BigInt un; un.limbs = u; un.normalize();
    vn = shl_bits(vn, s);
    un = shl_bits(un, s);
    v = vn.limbs;
    u = un.limbs;

    // Ensure u has one extra limb
    if (u.size() == x_in.limbs.size()) u.push_back(0);
    if (u.size() < n + m + 1) u.resize(n + m + 1, 0);

    std::vector<std::uint64_t> q(m + 1, 0);

    for (std::size_t jj = m + 1; jj-- > 0;) {
        std::size_t j = jj;

        // Estimate qhat using top two limbs
        unsigned __int128 numerator = ((unsigned __int128)u[j + n] << 64) | u[j + n - 1];
        std::uint64_t v1 = v[n - 1];
        std::uint64_t v2 = v[n - 2];

        std::uint64_t qhat = (std::uint64_t)(numerator / v1);
        std::uint64_t rhat = (std::uint64_t)(numerator % v1);

        // Fix qhat if too big
        while (qhat == ~0ULL ||
               (unsigned __int128)qhat * (unsigned __int128)v2 >
                   (((unsigned __int128)rhat << 64) | u[j + n - 2])) {
            qhat--;
            unsigned __int128 rr = (unsigned __int128)rhat + v1;
            rhat = (std::uint64_t)rr;
            if (rr >> 64) break;
        }

        // Multiply-subtract
        bool under = submul_at(u, v, j, qhat);
        if (under) {
            // qhat was 1 too large
            qhat--;
            add_at(u, v, j);
        }

        q[j] = qhat;
    }

    // remainder = (u[0..n-1] >> s)
    BigInt rem; rem.limbs.assign(u.begin(), u.begin() + n);
    rem.normalize();
    rem = shr_bits(rem, s);
    rem.normalize();
    if (!rem.limbs.empty()) throw std::runtime_error("div_exact: non-zero remainder");

    BigInt out; out.limbs = q;
    out.normalize();
    return out;
}

} // namespace bi
