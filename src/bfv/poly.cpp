#include "bfv/poly.h"
#include <stdexcept>

namespace bfv
{
    static void ensure_same_params(const Poly &x, const Poly &y)
    {
        if (!x.p || !y.p)
            throw std::runtime_error("Poly has null params");
        if (x.p != y.p)
            throw std::runtime_error("Poly params mismatch (must reference same Params object)");
    }

    Poly add(const Poly &x, const Poly &y)
    {
        ensure_same_params(x, y);
        Poly r(*x.p);
        const std::uint64_t q = x.q();
        for (std::size_t i = 0; i < x.N(); ++i)
            r[i] = add_mod(x[i], y[i], q);
        return r;
    }

    Poly sub(const Poly &x, const Poly &y)
    {
        ensure_same_params(x, y);
        Poly r(*x.p);
        const std::uint64_t q = x.q();
        for (std::size_t i = 0; i < x.N(); ++i)
            r[i] = sub_mod(x[i], y[i], q);
        return r;
    }

    // Negacyclic schoolbook: (a*b) mod (x^N + 1)
    // (x^N == -1) => terms that wrap around get a minus sign.
    Poly mul_negacyclic_schoolbook(const Poly &x, const Poly &y)
    {
        ensure_same_params(x, y);
        Poly r(*x.p);
        const std::uint64_t q = x.q();
        const std::size_t N = x.N();

        for (std::size_t i = 0; i < N; ++i)
        {
            for (std::size_t j = 0; j < N; ++j)
            {
                std::size_t k = i + j;
                std::uint64_t prod = mul_mod(x[i], y[j], q);

                if (k < N)
                {
                    r[k] = add_mod(r[k], prod, q);
                }
                else
                {
                    // wrap with sign flip: x^(k) = x^(k-N) * x^N = -x^(k-N)
                    k -= N;
                    r[k] = sub_mod(r[k], prod, q);
                }
            }
        }
        return r;
    }

}
