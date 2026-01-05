#include "bfv/bfv.h"
#include "bfv/ntt.h"
#include <cassert>
#include <iostream>

int main()
{
    // Keep toy params where enc/dec passed
    bfv::Params P{8, 97, 2};
    std::mt19937_64 rng(123);

    auto kp = bfv::keygen_noiseless(P, rng);

    for (std::uint64_t a = 0; a < P.t; ++a)
    {
        for (std::uint64_t b = 0; b < P.t; ++b)
        {
            auto ca = bfv::encrypt_const_noiseless(P, kp.pk, a, rng);
            auto cb = bfv::encrypt_const_noiseless(P, kp.pk, b, rng);

            auto prod3 = bfv::mul_raw(P, ca, cb);

            // NOTE: we cannot decrypt Ciphertext3 yet (needs relinearization OR extended decrypt),
            // so we do a temporary check by manually evaluating:
            //
            // v = c0 + c1*s + c2*s^2  (mod q), then decode.
            //
            // This is the "true" decrypt for 3-part ciphertexts.
            auto c1s = bfv::mul_negacyclic_ntt(prod3.c1, kp.sk.s);
            auto s2 = bfv::mul_negacyclic_ntt(kp.sk.s, kp.sk.s);
            auto c2s2 = bfv::mul_negacyclic_ntt(prod3.c2, s2);

            bfv::Poly v(P);
            for (std::size_t i = 0; i < P.N; ++i)
            {
                auto t0 = bfv::add_mod(prod3.c0[i], c1s[i], P.q);
                v[i] = bfv::add_mod(t0, c2s2[i], P.q);
            }

            // Decode constant term with same rounding as decrypt_const
            const std::uint64_t Delta = P.q / P.t;
            auto center = [&](std::uint64_t x) -> std::int64_t
            {
                if (x > P.q / 2)
                    return (std::int64_t)x - (std::int64_t)P.q;
                return (std::int64_t)x;
            };

            std::int64_t v0c = center(v[0]);
            std::int64_t Delta_i = (std::int64_t)(P.q / P.t);
            std::int64_t Delta2_i = Delta_i * Delta_i;

            std::int64_t m_hat_i = (v0c >= 0)
                                       ? (v0c + Delta2_i / 2) / Delta2_i
                                       : (v0c - Delta2_i / 2) / Delta2_i;

            std::uint64_t m_hat = (std::uint64_t)((m_hat_i % (std::int64_t)P.t + (std::int64_t)P.t) % (std::int64_t)P.t);
            m_hat %= P.t;

            assert(m_hat == ((a * b) % P.t));
        }
    }

    std::cout << "BFV ciphertext-ciphertext multiply (toy, 3-part decrypt) passed.\n";
    return 0;
}
