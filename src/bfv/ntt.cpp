#include "bfv/ntt.h"
#include <stdexcept>
#include <vector>
#include <cstdint>
#include <algorithm>

namespace bfv
{
  static std::uint64_t modexp(std::uint64_t a, std::uint64_t e, std::uint64_t mod)
  {
    std::uint64_t r = 1 % mod;
    while (e)
    {
      if (e & 1)
        r = (std::uint64_t)((__uint128_t)r * a % mod);
      a = (std::uint64_t)((__uint128_t)a * a % mod);
      e >>= 1;
    }
    return r;
  }

  static std::uint64_t modinv_prime(std::uint64_t a, std::uint64_t mod)
  {
    // mod must be prime
    return modexp(a, mod - 2, mod);
  }

  static void bit_reverse(std::vector<std::uint64_t> &a)
  {
    std::size_t n = a.size();
    for (std::size_t i = 1, j = 0; i < n; ++i)
    {
      std::size_t bit = n >> 1;
      for (; j & bit; bit >>= 1)
        j ^= bit;
      j |= bit;
      if (i < j)
        std::swap(a[i], a[j]);
    }
  }

  // Find a primitive root modulo prime q (generator of multiplicative group)
  static std::uint64_t primitive_root(std::uint64_t q)
  {
    // factor q-1 (trial division is fine for q up to ~64-bit for our demo primes)
    std::uint64_t phi = q - 1;
    std::vector<std::uint64_t> factors;

    std::uint64_t x = phi;
    for (std::uint64_t p = 2; p * p <= x; ++p)
    {
      if (x % p == 0)
      {
        factors.push_back(p);
        while (x % p == 0)
          x /= p;
      }
    }
    if (x > 1)
      factors.push_back(x);

    for (std::uint64_t g = 2; g < q; ++g)
    {
      bool ok = true;
      for (auto f : factors)
      {
        if (modexp(g, phi / f, q) == 1)
        {
          ok = false;
          break;
        }
      }
      if (ok)
        return g;
    }
    throw std::runtime_error("primitive_root: failed");
  }

  // Compute psi: primitive 2N-th root where psi^N = -1 mod q
  static std::uint64_t find_psi(std::size_t N, std::uint64_t q)
  {
    std::size_t twoN = 2 * N;
    if ((q - 1) % twoN != 0)
      throw std::runtime_error("NTT requires q ≡ 1 (mod 2N)");

    std::uint64_t g = primitive_root(q);
    std::uint64_t psi = modexp(g, (q - 1) / (std::uint64_t)twoN, q);

    // verify psi^(2N)=1 and psi^N = q-1 (-1)
    if (modexp(psi, (std::uint64_t)twoN, q) != 1)
      throw std::runtime_error("psi not 2N root");
    if (modexp(psi, (std::uint64_t)N, q) != (q - 1))
      throw std::runtime_error("psi^N != -1");

    return psi;
  }

  // Cyclic NTT size N with primitive N-th root "root"
  static void ntt_cyclic(std::vector<std::uint64_t> &a, std::size_t N, std::uint64_t q, std::uint64_t root)
  {
    bit_reverse(a);

    for (std::size_t len = 2; len <= N; len <<= 1)
    {
      std::size_t half = len >> 1;
      std::uint64_t wlen = modexp(root, (std::uint64_t)(N / len), q);

      for (std::size_t i = 0; i < N; i += len)
      {
        std::uint64_t w = 1;
        for (std::size_t j = 0; j < half; ++j)
        {
          std::uint64_t u = a[i + j];
          std::uint64_t v = (std::uint64_t)((__uint128_t)a[i + j + half] * w % q);

          std::uint64_t t0 = u + v;
          if (t0 >= q)
            t0 -= q;
          a[i + j] = t0;

          std::uint64_t t1 = (u >= v) ? (u - v) : (u + q - v);
          a[i + j + half] = t1;

          w = (std::uint64_t)((__uint128_t)w * wlen % q);
        }
      }
    }
  }

  static void intt_cyclic(std::vector<std::uint64_t> &a, std::size_t N, std::uint64_t q, std::uint64_t root)
  {
    std::uint64_t inv_root = modinv_prime(root, q);

    bit_reverse(a);

    for (std::size_t len = 2; len <= N; len <<= 1)
    {
      std::size_t half = len >> 1;
      std::uint64_t wlen = modexp(inv_root, (std::uint64_t)(N / len), q);

      for (std::size_t i = 0; i < N; i += len)
      {
        std::uint64_t w = 1;
        for (std::size_t j = 0; j < half; ++j)
        {
          std::uint64_t u = a[i + j];
          std::uint64_t v = (std::uint64_t)((__uint128_t)a[i + j + half] * w % q);

          std::uint64_t t0 = u + v;
          if (t0 >= q)
            t0 -= q;
          a[i + j] = t0;

          std::uint64_t t1 = (u >= v) ? (u - v) : (u + q - v);
          a[i + j + half] = t1;

          w = (std::uint64_t)((__uint128_t)w * wlen % q);
        }
      }
    }

    std::uint64_t invN = modinv_prime((std::uint64_t)N, q);
    for (auto &x : a)
      x = (std::uint64_t)((__uint128_t)x * invN % q);
  }

  void ntt(const Params &p, Poly &a)
  {
    if (p.N == 0 || (p.N & (p.N - 1)) != 0)
      throw std::runtime_error("NTT requires N power of two");
    if (a.N() != p.N || a.q() != p.q)
      throw std::runtime_error("NTT params mismatch");

    std::uint64_t psi = find_psi(p.N, p.q);
    std::uint64_t omega = (std::uint64_t)((__uint128_t)psi * psi % p.q); // psi^2 has order N

    // pre-twist: a[i] *= psi^i
    for (std::size_t i = 0; i < p.N; ++i)
    {
      std::uint64_t t = modexp(psi, (std::uint64_t)i, p.q);
      a[i] = (std::uint64_t)((__uint128_t)a[i] * t % p.q);
    }

    ntt_cyclic(a.a, p.N, p.q, omega);
  }

  void intt(const Params &p, Poly &a)
  {
    if (p.N == 0 || (p.N & (p.N - 1)) != 0)
      throw std::runtime_error("NTT requires N power of two");
    if (a.N() != p.N || a.q() != p.q)
      throw std::runtime_error("INTT params mismatch");

    std::uint64_t psi = find_psi(p.N, p.q);
    std::uint64_t omega = (std::uint64_t)((__uint128_t)psi * psi % p.q);

    intt_cyclic(a.a, p.N, p.q, omega);

    // post-twist by psi^{-i}
    std::uint64_t inv_psi = modinv_prime(psi, p.q);
    for (std::size_t i = 0; i < p.N; ++i)
    {
      std::uint64_t t = modexp(inv_psi, (std::uint64_t)i, p.q);
      a[i] = (std::uint64_t)((__uint128_t)a[i] * t % p.q);
    }
  }

  Poly mul_negacyclic_ntt(const Params &p, const Poly &a, const Poly &b)
  {
    Poly A = a;
    Poly B = b;

    ntt(p, A);
    ntt(p, B);

    for (std::size_t i = 0; i < p.N; ++i)
      A[i] = (std::uint64_t)((__uint128_t)A[i] * B[i] % p.q);

    intt(p, A);
    return A;
  }
}
