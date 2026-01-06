#include "bfv/bfv.h"
#include "bfv/relin.h"
#include "bfv/noise.h"
#include <fstream>
#include <iostream>

int main()
{
    bfv::Params P{1024, 12289, 2};
    std::mt19937_64 rng(123);

    auto kp = bfv::keygen(P, rng);
    auto rlk = bfv::relin_keygen_noiseless(P, kp.sk, rng, 4);

    std::ofstream out("bfv_noise.csv");
    out << "op,a,b,noise_inf,noise_l2_sq\n";

    for (std::uint64_t a = 0; a < P.t; ++a)
    {
        for (std::uint64_t b = 0; b < P.t; ++b)
        {
            auto ca = bfv::encrypt_const(P, kp.pk, a, rng);
            auto cb = bfv::encrypt_const(P, kp.pk, b, rng);

            auto na = bfv::noise_const(P, kp.sk, ca, a);
            out << "enc," << a << "," << b << "," << na.inf << "," << na.l2_sq << "\n";

            auto cc_add = bfv::add(P, ca, cb);
            auto nadd = bfv::noise_const(P, kp.sk, cc_add, (a + b) % P.t);
            out << "add," << a << "," << b << "," << nadd.inf << "," << nadd.l2_sq << "\n";

            auto ct3 = bfv::mul_raw(P, ca, cb);
            auto nraw = bfv::noise3_const(P, kp.sk, ct3, (a * b) % P.t);
            out << "mul_raw," << a << "," << b << "," << nraw.inf << "," << nraw.l2_sq << "\n";

            auto cc = bfv::relinearize(P, ct3, rlk);
            auto nrelin = bfv::noise_const(P, kp.sk, cc, (a * b) % P.t);
            out << "relin," << a << "," << b << "," << nrelin.inf << "," << nrelin.l2_sq << "\n";
        }
    }

    std::cout << "Wrote bfv_noise.csv (noise stats).\n";
    return 0;
}
