#include "bfv/noise.h"
#include "bfv/ntt.h"
#include <cstdlib> // llabs
#include <stdexcept>

namespace bfv
{
    static inline std::int64_t center_u64(std::uint64_t x, std::uint64_t q)
    {
        if (x > q / 2)
            return (std::int64_t)x - (std::int64_t)q;
        return (std::int64_t)x;
    }

    static inline std::uint64_t add_mod_u64(std::uint64_t a, std::uint64_t b, std::uint64_t q)
    {
        std::uint64_t s = a + b;
        if (s >= q)
            s -= q;
        return s;
    }

    static inline std::uint64_t sub_mod_u64(std::uint64_t a, std::uint64_t b, std::uint64_t q)
    {
        return (a >= b) ? (a - b) : (a + q - b);
    }

    static NoiseStats stats_from_error_poly(const Params &p, const Poly &e_centered)
    {
        NoiseStats st;
        for (std::size_t i = 0; i < p.N; ++i)
        {
            std::int64_t ei = e_centered.a[i];
            std::int64_t ae = std::llabs(ei);
            if (ae > st.inf)
                st.inf = ae;
            st.l2_sq += (std::uint64_t)(ei * ei);
        }
        return st;
    }

    NoiseStats noise_const(const Params &p, const SecretKey &sk, const Ciphertext &ct, std::uint64_t m)
    {
        // v = c0 + c1*s mod q
        Poly c1s = mul_negacyclic_ntt(p, ct.c1, sk.s);
        Poly v(p);
        for (std::size_t i = 0; i < p.N; ++i)
            v[i] = add_mod(ct.c0[i], c1s[i], p.q);

        // expected = Δ*m at coeff 0, 0 elsewhere
        const std::uint64_t Delta = p.q / p.t;
        std::uint64_t exp0 = (Delta * (m % p.t)) % p.q;

        // e = center(v - expected)
        Poly e(p);
        for (std::size_t i = 0; i < p.N; ++i)
        {
            std::uint64_t expi = (i == 0) ? exp0 : 0;
            std::uint64_t diff = sub_mod_u64(v[i], expi, p.q);
            e.a[i] = center_u64(diff, p.q);
        }

        return stats_from_error_poly(p, e);
    }

    NoiseStats noise3_const(const Params &p, const SecretKey &sk, const Ciphertext3 &ct3, std::uint64_t m)
    {
        // v = c0 + c1*s + c2*s^2 mod q
        Poly c1s = mul_negacyclic_ntt(p, ct3.c1, sk.s);
        Poly s2 = mul_negacyclic_ntt(p, sk.s, sk.s);
        Poly c2s2 = mul_negacyclic_ntt(p, ct3.c2, s2);

        Poly v(p);
        for (std::size_t i = 0; i < p.N; ++i)
        {
            auto t0 = add_mod(ct3.c0[i], c1s[i], p.q);
            v[i] = add_mod(t0, c2s2[i], p.q);
        }

        // expected = Δ^2*m at coeff 0 for raw multiply (scale grows)
        const std::uint64_t Delta = p.q / p.t;
        std::uint64_t Delta2 = (Delta * Delta) % p.q;
        std::uint64_t exp0 = (Delta2 * (m % p.t)) % p.q;

        Poly e(p);
        for (std::size_t i = 0; i < p.N; ++i)
        {
            std::uint64_t expi = (i == 0) ? exp0 : 0;
            std::uint64_t diff = sub_mod_u64(v[i], expi, p.q);
            e.a[i] = center_u64(diff, p.q);
        }

        return stats_from_error_poly(p, e);
    }
}
