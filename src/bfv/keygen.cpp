#include "bfv/bfv.h"
#include "bfv/ntt.h"
#include <stdexcept>

namespace bfv
{
    Poly sample_uniform_q(const Params &p, std::mt19937_64 &rng)
    {
        Poly r(p);
        std::uniform_int_distribution<std::uint64_t> dist(0, p.q - 1);
        for (std::size_t i = 0; i < p.N; ++i)
            r[i] = dist(rng);
        return r;
    }

    Poly sample_ternary(const Params &p, std::mt19937_64 &rng)
    {
        Poly r(p);
        // {-1,0,1} mapped into mod q as {q-1, 0, 1}
        std::uniform_int_distribution<int> dist(0, 2);
        for (std::size_t i = 0; i < p.N; ++i)
        {
            int v = dist(rng) - 1; // -1,0,1
            if (v == -1)
                r[i] = p.q - 1;
            else
                r[i] = (std::uint64_t)v;
        }
        return r;
    }

    Poly sample_small_error(const Params &p, std::mt19937_64 &rng)
    {
        // For now, reuse ternary as noise (classic starter choice).
        return sample_ternary(p, rng);
    }

    KeyPair keygen(const Params &p, std::mt19937_64 &rng)
    {
        // For now, we only support the toy NTT params, to keep correctness tight.
        if (p.N != 8 || p.q != 97)
            throw std::runtime_error("toy keygen expects Params{N=8,q=97,...} (generalization comes next)");

        KeyPair kp;
        kp.sk.s = sample_ternary(p, rng);

        Poly a = sample_uniform_q(p, rng);
        Poly e = sample_small_error(p, rng);

        // Compute a*s using fast negacyclic multiplication (NTT)
        Poly as = mul_negacyclic_ntt(a, kp.sk.s);

        // b = -(a*s + e) mod q
        Poly b(p);
        for (std::size_t i = 0; i < p.N; ++i)
        {
            std::uint64_t t = add_mod(as[i], e[i], p.q);
            b[i] = (t == 0) ? 0 : (p.q - t);
        }

        kp.pk.a = a;
        kp.pk.b = b;
        return kp;
    }

    KeyPair keygen_noiseless(const Params &p, std::mt19937_64 &rng)
    {
        if (p.N != 8 || p.q != 97)
            throw std::runtime_error("toy keygen_noiseless expects Params{N=8,q=97}");

        KeyPair kp;
        kp.sk.s = sample_ternary(p, rng);

        Poly a = sample_uniform_q(p, rng);

        // e = 0
        Poly e(p);
        for (std::size_t i = 0; i < p.N; ++i)
            e[i] = 0;

        Poly as = mul_negacyclic_ntt(a, kp.sk.s);

        // b = -(a*s + e) = -(a*s)
        Poly b(p);
        for (std::size_t i = 0; i < p.N; ++i)
        {
            std::uint64_t t0 = as[i]; // e=0
            b[i] = (t0 == 0) ? 0 : (p.q - t0);
        }

        kp.pk.a = a;
        kp.pk.b = b;
        return kp;
    }
}
