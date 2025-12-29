#include "bigint/prime.h"
#include "bigint/mod_arith.h" // mul_mod, pow_mod, pow_mod_bigexp (we use pow_mod_bigexp)
#include <algorithm>
#include <vector>

namespace bi
{

    // ---------- Small helpers on BigInt ----------

    static bool is_zero(const BigInt &x)
    {
        return x.limbs.empty();
    }

    static bool is_one(const BigInt &x)
    {
        return x.limbs.size() == 1 && x.limbs[0] == 1;
    }

    static bool is_even(const BigInt &x)
    {
        if (x.limbs.empty())
            return true;
        return (x.limbs[0] & 1ULL) == 0;
    }

    static BigInt add_u64(BigInt x, std::uint64_t v)
    {
        std::uint64_t carry = v;
        std::size_t i = 0;
        while (carry)
        {
            if (i >= x.limbs.size())
                x.limbs.push_back(0);
            unsigned __int128 sum = (unsigned __int128)x.limbs[i] + carry;
            x.limbs[i] = (std::uint64_t)sum;
            carry = (std::uint64_t)(sum >> 64);
            ++i;
        }
        x.normalize();
        return x;
    }

    static BigInt sub_u64(BigInt x, std::uint64_t v)
    {
        // assume x >= v
        std::uint64_t borrow = v;
        std::size_t i = 0;
        while (borrow)
        {
            std::uint64_t xi = (i < x.limbs.size()) ? x.limbs[i] : 0;
            std::uint64_t newv = xi - borrow;
            // borrow occurs if xi < borrow
            borrow = (xi < borrow) ? 1 : 0;
            x.limbs[i] = newv;
            ++i;
        }
        x.normalize();
        return x;
    }

    static void shr1_inplace(BigInt &x)
    {
        // x >>= 1
        std::uint64_t carry = 0;
        for (std::size_t i = x.limbs.size(); i-- > 0;)
        {
            std::uint64_t new_carry = x.limbs[i] & 1ULL;
            x.limbs[i] = (x.limbs[i] >> 1) | (carry << 63);
            carry = new_carry;
        }
        x.normalize();
    }

    static std::size_t bitlen(const BigInt &x)
    {
        if (x.limbs.empty())
            return 0;
        std::size_t i = x.limbs.size() - 1;
        std::uint64_t top = x.limbs[i];
        std::size_t bits = i * 64;
        while (top)
        {
            top >>= 1;
            ++bits;
        }
        return bits;
    }

    // Rejection-sample uniform random in [0, limit] where limit >= 0
    static BigInt random_below_or_equal(const BigInt &limit, std::mt19937_64 &rng)
    {
        BigInt L = limit;
        L.normalize();
        if (L.limbs.empty())
            return BigInt::from_u64(0);

        const std::size_t k = L.limbs.size();
        const std::size_t top_bits = (bitlen(L) - 1) % 64 + 1;
        const std::uint64_t top_mask = (top_bits == 64) ? ~0ULL : ((1ULL << top_bits) - 1);

        while (true)
        {
            BigInt x;
            x.limbs.resize(k);
            for (std::size_t i = 0; i < k; ++i)
                x.limbs[i] = rng();
            x.limbs[k - 1] &= top_mask;
            x.normalize();
            if (cmp(x, L) <= 0)
                return x;
        }
    }

    // Returns random a in [low, high], assuming low <= high
    static BigInt random_range(const BigInt &low, const BigInt &high, std::mt19937_64 &rng)
    {
        // sample t in [0, high-low], return low + t
        BigInt diff = high;
        sub_inplace(diff, low);
        BigInt t = random_below_or_equal(diff, rng);
        // return low + t
        BigInt res = low;
        // res += t (simple add)
        const std::size_t n = std::max(res.limbs.size(), t.limbs.size());
        res.limbs.resize(n, 0);
        std::uint64_t carry = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
            unsigned __int128 sum = (unsigned __int128)res.limbs[i] + (unsigned __int128)(i < t.limbs.size() ? t.limbs[i] : 0) + carry;
            res.limbs[i] = (std::uint64_t)sum;
            carry = (std::uint64_t)(sum >> 64);
        }
        if (carry)
            res.limbs.push_back(carry);
        res.normalize();
        return res;
    }

    // ---------- Miller–Rabin ----------

    static bool trial_div_small_primes(const BigInt &n)
    {
        // Quick filter: if n mod p == 0 for small p, composite (unless n==p).
        static const std::uint32_t primes[] = {
            3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

        // if n fits in u64, do direct
        if (n.limbs.size() <= 1)
        {
            std::uint64_t v = n.limbs.empty() ? 0 : n.limbs[0];
            for (auto p : primes)
            {
                if (v == p)
                    return false; // not composite
                if (v % p == 0)
                    return true; // composite
            }
            return false;
        }

        // for big n, compute n mod p using limb method
        for (auto p : primes)
        {
            std::uint64_t r = 0;
            for (std::size_t i = n.limbs.size(); i-- > 0;)
            {
                // r = (r * 2^64 + limb) mod p
                // since p fits in 32 bits, we can do this safely with __int128
                unsigned __int128 acc = (unsigned __int128)r << 64;
                acc += n.limbs[i];
                r = (std::uint64_t)(acc % p);
            }
            if (r == 0)
                return true; // divisible by p => composite (n != p since n is huge)
        }
        return false;
    }

    bool is_probable_prime(const BigInt &n_in, int rounds)
    {
        BigInt n = n_in;
        n.normalize();

        // Handle small n
        if (is_zero(n))
            return false;
        if (n.limbs.size() == 1)
        {
            std::uint64_t v = n.limbs[0];
            if (v < 2)
                return false;
            if (v == 2 || v == 3)
                return true;
            if ((v & 1ULL) == 0)
                return false;
        }

        // even => composite
        if (is_even(n))
            return false;

        // small prime trial division (fast composite filter)
        if (trial_div_small_primes(n))
            return false;

        // Write n-1 = 2^s * d with d odd
        BigInt n_minus_1 = n;
        sub_inplace(n_minus_1, BigInt::from_u64(1));

        BigInt d = n_minus_1;
        std::size_t s = 0;
        while (!is_zero(d) && is_even(d))
        {
            shr1_inplace(d);
            ++s;
        }

        // Setup RNG (deterministic seed for reproducibility in tests/bench)
        std::mt19937_64 rng(123456789ULL);

        // Bounds for a: [2, n-2]
        BigInt two = BigInt::from_u64(2);
        BigInt n_minus_2 = n;
        sub_inplace(n_minus_2, two);

        for (int r = 0; r < rounds; ++r)
        {
            // pick random a in [2, n-2]
            BigInt a = random_range(two, n_minus_2, rng);

            // x = a^d mod n
            BigInt x = pow_mod_bigexp(a, d, n);

            // If x == 1 or x == n-1 => pass round
            if (is_one(x) || cmp(x, n_minus_1) == 0)
                continue;

            bool witness = true; // assume composite unless we hit n-1
            for (std::size_t i = 1; i < s; ++i)
            {
                x = mul_mod(x, x, n);
                if (cmp(x, n_minus_1) == 0)
                {
                    witness = false;
                    break;
                }
                if (is_one(x))
                    return false; // non-trivial sqrt of 1 => composite
            }
            if (witness)
                return false; // composite found
        }

        return true; // probably prime
    }

    bool is_probable_prime(const BigInt &n_in, std::mt19937_64 &rng, int rounds)
    {
        BigInt n = n_in;
        n.normalize();

        // Handle small n
        if (is_zero(n))
            return false;
        if (n.limbs.size() == 1)
        {
            std::uint64_t v = n.limbs[0];
            if (v < 2)
                return false;
            if (v == 2 || v == 3)
                return true;
            if ((v & 1ULL) == 0)
                return false;
        }

        // even => composite
        if (is_even(n))
            return false;

        // small prime trial division (fast composite filter)
        if (trial_div_small_primes(n))
            return false;

        // Write n-1 = 2^s * d with d odd
        BigInt n_minus_1 = n;
        sub_inplace(n_minus_1, BigInt::from_u64(1));

        BigInt d = n_minus_1;
        std::size_t s = 0;
        while (!is_zero(d) && is_even(d))
        {
            shr1_inplace(d);
            ++s;
        }

        // Setup RNG (deterministic seed for reproducibility in tests/bench)
        // std::mt19937_64 rng(123456789ULL);

        // Bounds for a: [2, n-2]
        BigInt two = BigInt::from_u64(2);
        BigInt n_minus_2 = n;
        sub_inplace(n_minus_2, two);

        for (int r = 0; r < rounds; ++r)
        {
            // pick random a in [2, n-2]
            BigInt a = random_range(two, n_minus_2, rng);

            // x = a^d mod n
            BigInt x = pow_mod_bigexp(a, d, n);

            // If x == 1 or x == n-1 => pass round
            if (is_one(x) || cmp(x, n_minus_1) == 0)
                continue;

            bool witness = true; // assume composite unless we hit n-1
            for (std::size_t i = 1; i < s; ++i)
            {
                x = mul_mod(x, x, n);
                if (cmp(x, n_minus_1) == 0)
                {
                    witness = false;
                    break;
                }
                if (is_one(x))
                    return false; // non-trivial sqrt of 1 => composite
            }
            if (witness)
                return false; // composite found
        }

        return true; // probably prime
    }

    // ---------- Prime generation ----------

    BigInt random_odd_bits(std::size_t bits, std::mt19937_64 &rng)
    {
        if (bits < 2)
            bits = 2;

        const std::size_t limbs = (bits + 63) / 64;
        const std::size_t top_bits = (bits - 1) % 64 + 1;

        BigInt x;
        x.limbs.resize(limbs);
        for (std::size_t i = 0; i < limbs; ++i)
            x.limbs[i] = rng();

        // Mask top limb to desired bit length
        std::uint64_t mask = (top_bits == 64) ? ~0ULL : ((1ULL << top_bits) - 1);
        x.limbs[limbs - 1] &= mask;

        // Set top bit to ensure exact bit length
        x.limbs[limbs - 1] |= (1ULL << (top_bits - 1));

        // Make odd
        x.limbs[0] |= 1ULL;

        x.normalize();
        return x;
    }

    BigInt generate_prime(std::size_t bits, std::mt19937_64& rng, int rounds)
    {
        while (true)
        {
            BigInt cand = random_odd_bits(bits, rng);
            if (is_probable_prime(cand, rng, rounds))
                return cand;
        }
    }

} // namespace bi
