#pragma once
#include <cstddef>
#include <cstdint>

namespace bfv
{
    struct Params
    {
        std::size_t N;     // polynomial degree (power of two)
        std::uint64_t q;   // ciphertext modulus
        std::uint64_t t;   // plaintext modulus

        // convenience
        std::uint64_t delta() const
        {
            return q / t;
        }
    };
}