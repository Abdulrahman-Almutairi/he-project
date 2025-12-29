#pragma once
#include "bigint/bigint.h"
#include <cstddef>
#include <cstdint>
#include <random>

namespace bi
{

    // Miller–Rabin probable prime test.
    // returns true => probably prime, false => definitely composite.
    bool is_probable_prime(const BigInt &n, std::mt19937_64 &rng, int rounds = 20);

    // Random odd BigInt with exactly `bits` bits (top bit set).
    BigInt random_odd_bits(std::size_t bits, std::mt19937_64 &rng);

    // Generate a probable prime with `bits` bits using Miller–Rabin.
    BigInt generate_prime(std::size_t bits, std::mt19937_64 &rng, int rounds = 20);

} 