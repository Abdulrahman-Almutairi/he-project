#pragma once
#include "bfv/params.h"
#include <cstdint>
#include <vector>

namespace bfv
{
    // Polynomial in R_q = Z_q[x] / (x^N + 1)
    // Stored as N coefficients modulo q (uint64_t), little_endian by degree.
    struct Poly
    {
        const Params *p = nullptr;
        std::vector<std::uint64_t> a; // size N

        Poly() = default;
        explicit Poly(const Params &params) : p(&params), a(params.N, 0) {}

        std::size_t N() const { return p ? p->N : 0; }
        std::uint64_t q() const { return p ? p->q : 0; }
        std::uint64_t t() const { return p ? p->t : 0; }

        std::uint64_t &operator[](std::size_t i) { return a[i]; }
        const std::uint64_t &operator[](std::size_t i) const { return a[i]; }

        void set_params(const Params &params)
        {
            p = &params;
            a.assign(params.N, 0);
        }
    };

    // --- coefficient helpers mod q ---

    inline std::uint64_t add_mod(std::uint64_t x, std::uint64_t y, std::uint64_t q)
    {
        unsigned __int128 s = (unsigned __int128)x + (unsigned __int128)y;
        std::uint64_t r = (std::uint64_t)s;
        if (r >= q || (s >> 64))
            r %= q;
        return r;
    }

    inline std::uint64_t sub_mod(std::uint64_t x, std::uint64_t y, std::uint64_t q)
    {
        return (x >= y) ? (x - y) : (q - (y - x));
    }

    inline std::uint64_t mul_mod(std::uint64_t x, std::uint64_t y, std::uint64_t q)
    {
        unsigned __int128 p = (unsigned __int128)x * (unsigned __int128)y;
        return (std::uint64_t)(p % q);
    }

    // -------- polynomial ops in R_q (coefficient-wise) --------

    // c = a + b (mod q)
    Poly add(const Poly &a, const Poly &b);

    // c = a - b (mod q)
    Poly sub(const Poly &a, const Poly &b);

    // Multiply in R_q = Z_q[x]/(x^N+1).
    // For now, we implement a slow schoolbook negacyclic convolution.
    // Later we will replace with NTT-based multiplication.
    Poly mul_negacyclic_schoolbook(const Poly &a, const Poly &b);
}
