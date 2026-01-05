#include "bfv/bfv.h"
#include "bfv/ntt.h"
#include <stdexcept>

namespace bfv
{
    static inline std::uint64_t mod_q(std::int64_t x, std::uint64_t q)
    {
        std::int64_t r = x % (std::int64_t)q;
        if (r < 0)
            r += (std::int64_t)q;
        return (std::uint64_t)r;
    }

    Poly encode_const(const Params &p, std::uint64_t m)
    {
        Poly r(p);
        r[0] = m % p.t;
        return r;
    }

    std::uint64_t decode_const(const Params &p, const Poly &m)
    {
        return m.a.empty() ? 0 : (m[0] % p.t);
    }

    Ciphertext encrypt_const(const Params &p, const PublicKey &pk, std::uint64_t m, std::mt19937_64 &rng)
    {
        if (p.N != 8 || p.q != 97)
            throw std::runtime_error("toy encrypt expects Params{N=8,q=97}");

        // message poly in R_t (constant)
        Poly mp = encode_const(p, m);

        // scale into R_q: Δ * m  (Δ = floor(q/t))
        const std::uint64_t Delta = p.q / p.t;

        Poly scaled(p);
        for (std::size_t i = 0; i < p.N; ++i)
            scaled[i] = 0;
        scaled[0] = (mp[0] * Delta) % p.q;

        // sample u, e1, e2 small
        Poly u = sample_ternary(p, rng);
        Poly e1 = sample_small_error(p, rng);
        Poly e2 = sample_small_error(p, rng);

        // c0 = b*u + e1 + scaled
        Poly bu = mul_negacyclic_ntt(pk.b, u);
        Poly c0(p);
        for (std::size_t i = 0; i < p.N; ++i)
        {
            std::uint64_t t0 = add_mod(bu[i], e1[i], p.q);
            c0[i] = add_mod(t0, scaled[i], p.q);
        }

        // c1 = a*u + e2
        Poly au = mul_negacyclic_ntt(pk.a, u);
        Poly c1(p);
        for (std::size_t i = 0; i < p.N; ++i)
            c1[i] = add_mod(au[i], e2[i], p.q);

        Ciphertext ct(p);
        ct.c0 = c0;
        ct.c1 = c1;
        return ct;
    }

    Ciphertext encrypt_const_noiseless(const Params &p, const PublicKey &pk, std::uint64_t m, std::mt19937_64 &rng)
    {
        if (p.N != 8 || p.q != 97)
            throw std::runtime_error("toy encrypt expects Params{N=8,q=97}");

        // message poly in R_t (constant)
        Poly mp = encode_const(p, m);

        const std::uint64_t Delta = p.q / p.t;

        Poly scaled(p);
        for (std::size_t i = 0; i < p.N; ++i)
            scaled[i] = 0;
        scaled[0] = (mp[0] * Delta) % p.q;

        // sample u (keep randomness), but set e1=e2=0
        Poly u = sample_ternary(p, rng);

        Poly zero(p); // all zeros

        // c0 = b*u + scaled
        Poly bu = mul_negacyclic_ntt(pk.b, u);
        Poly c0(p);
        for (std::size_t i = 0; i < p.N; ++i)
            c0[i] = add_mod(bu[i], scaled[i], p.q);

        // c1 = a*u
        Poly au = mul_negacyclic_ntt(pk.a, u);
        Poly c1(p);
        for (std::size_t i = 0; i < p.N; ++i)
            c1[i] = au[i];

        Ciphertext ct(p);
        ct.c0 = c0;
        ct.c1 = c1;
        return ct;
    }

    std::uint64_t decrypt_const(const Params &p, const SecretKey &sk, const Ciphertext &ct)
    {
        if (p.N != 8 || p.q != 97)
            throw std::runtime_error("toy decrypt expects Params{N=8,q=97}");

        // v = c0 + c1*s mod q
        Poly c1s = mul_negacyclic_ntt(ct.c1, sk.s);
        Poly v(p);
        for (std::size_t i = 0; i < p.N; ++i)
            v[i] = add_mod(ct.c0[i], c1s[i], p.q);

        // Recover m0 by rounding v0/Δ to nearest integer mod t
        const std::uint64_t Delta = p.q / p.t;
        auto center = [&](std::uint64_t x) -> std::int64_t
        {
            // map [0,q) to (-q/2, q/2]
            if (x > p.q / 2)
                return (std::int64_t)x - (std::int64_t)p.q;
            return (std::int64_t)x;
        };

        std::int64_t v0c = center(v[0]);

        // nearest rounding of v0c / Delta
        std::int64_t m_hat_i = (v0c >= 0)
                                   ? (v0c + (std::int64_t)Delta / 2) / (std::int64_t)Delta
                                   : (v0c - (std::int64_t)Delta / 2) / (std::int64_t)Delta;

        std::uint64_t m_hat = (std::uint64_t)((m_hat_i % (std::int64_t)p.t + (std::int64_t)p.t) % (std::int64_t)p.t);
        return m_hat;
    }

}
