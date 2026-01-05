#include "bfv/bfv.h"
#include "bfv/ntt.h"
#include <stdexcept>

namespace bfv
{
    Ciphertext3 mul_raw(const Params &p, const Ciphertext &x, const Ciphertext &y)
    {
        Ciphertext3 r(p);

        // r0 = x0*y0
        r.c0 = mul_negacyclic_ntt(x.c0, y.c0);

        // r1 = x0*y1 + x1*y0
        Poly x0y1 = mul_negacyclic_ntt(x.c0, y.c1);
        Poly x1y0 = mul_negacyclic_ntt(x.c1, y.c0);
        r.c1 = add(x0y1, x1y0);

        // r2 = x1*y1
        r.c2 = mul_negacyclic_ntt(x.c1, y.c1);

        return r;
    }

    static inline std::uint64_t div_round_u64(std::uint64_t x, std::uint64_t d)
    {
        // nearest integer rounding: floor((x + d/2)/d)
        return (x + (d / 2)) / d;
    }

    static inline std::int64_t center_u64(std::uint64_t x, std::uint64_t q)
    {
        // Map x in [0,q) into the centered representative in (-q/2, q/2]
        // For odd q, q/2 is floor(q/2).
        if (x > q / 2)
            return (std::int64_t)x - (std::int64_t)q;
        return (std::int64_t)x;
    }

    static inline std::uint64_t mod_q_from_i64(std::int64_t x, std::uint64_t q)
    {
        std::int64_t r = x % (std::int64_t)q;
        if (r < 0)
            r += (std::int64_t)q;
        return (std::uint64_t)r;
    }

    static inline std::int64_t div_round_i64(std::int64_t x, std::int64_t d)
    {
        // Nearest integer rounding for signed values
        return (x >= 0) ? (x + d / 2) / d
                        : (x - d / 2) / d;
    }

    Ciphertext3 rescale_delta(const Params &p, const Ciphertext3 &ct3)
    {
        if (p.t == 0)
            throw std::runtime_error("t must be nonzero");

        const std::uint64_t Delta = p.q / p.t;
        if (Delta == 0)
            throw std::runtime_error("Delta became 0 (q < t)");

        const std::int64_t Delta_i = (std::int64_t)Delta;

        Ciphertext3 r(p);
        for (std::size_t i = 0; i < p.N; ++i)
        {
            // Center coefficients first, then divide with rounding, then map back mod q
            std::int64_t c0c = center_u64(ct3.c0[i], p.q);
            std::int64_t c1c = center_u64(ct3.c1[i], p.q);
            std::int64_t c2c = center_u64(ct3.c2[i], p.q);

            std::int64_t d0 = div_round_i64(c0c, Delta_i);
            std::int64_t d1 = div_round_i64(c1c, Delta_i);
            std::int64_t d2 = div_round_i64(c2c, Delta_i);

            r.c0[i] = mod_q_from_i64(d0, p.q);
            r.c1[i] = mod_q_from_i64(d1, p.q);
            r.c2[i] = mod_q_from_i64(d2, p.q);
        }

        return r;
    }

    Ciphertext3 mul(const Params &p, const Ciphertext &x, const Ciphertext &y)
    {
        auto raw = mul_raw(p, x, y);
        return rescale_delta(p, raw);
    }
}
