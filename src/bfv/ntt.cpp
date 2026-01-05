#include "bfv/ntt.h"
#include <vector>
#include <stdexcept>

namespace bfv
{
  // Toy fixed parameters for first correctness test
  static constexpr std::uint64_t q = 97;
  static constexpr std::size_t N = 8;

  // psi: primitive 2N-th root of unity (order 16) mod 97
  // For q=97, N=8, psi=8 works (order 16).
  static constexpr std::uint64_t psi = 8;

  // omega = psi^2 is primitive N-th root (order 8)
  static constexpr std::uint64_t omega = (psi * psi) % q;

  static std::uint64_t modexp(std::uint64_t a, std::uint64_t e)
  {
    std::uint64_t r = 1;
    while (e)
    {
      if (e & 1) r = (r * a) % q;
      a = (a * a) % q;
      e >>= 1;
    }
    return r;
  }

  static std::uint64_t modinv(std::uint64_t a)
  {
    // q is prime in this toy setup
    return modexp(a, q - 2);
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

  // Cyclic NTT of size N with root = omega (order N)
  static void ntt_cyclic(std::vector<std::uint64_t> &a, std::uint64_t root)
  {
    bit_reverse(a);

    for (std::size_t len = 2; len <= N; len <<= 1)
    {
      std::size_t half = len >> 1;
      std::uint64_t wlen = modexp(root, N / len);

      for (std::size_t i = 0; i < N; i += len)
      {
        std::uint64_t w = 1;
        for (std::size_t j = 0; j < half; ++j)
        {
          std::uint64_t u = a[i + j];
          std::uint64_t v = (a[i + j + half] * w) % q;

          a[i + j] = (u + v) % q;
          a[i + j + half] = (u + q - v) % q;

          w = (w * wlen) % q;
        }
      }
    }
  }

  // Inverse cyclic NTT
  static void intt_cyclic(std::vector<std::uint64_t> &a, std::uint64_t root)
  {
    // inverse root
    std::uint64_t inv_root = modinv(root);

    bit_reverse(a);

    for (std::size_t len = 2; len <= N; len <<= 1)
    {
      std::size_t half = len >> 1;
      std::uint64_t wlen = modexp(inv_root, N / len);

      for (std::size_t i = 0; i < N; i += len)
      {
        std::uint64_t w = 1;
        for (std::size_t j = 0; j < half; ++j)
        {
          std::uint64_t u = a[i + j];
          std::uint64_t v = (a[i + j + half] * w) % q;

          a[i + j] = (u + v) % q;
          a[i + j + half] = (u + q - v) % q;

          w = (w * wlen) % q;
        }
      }
    }

    // multiply by invN
    std::uint64_t invN = modinv(N);
    for (auto &x : a)
      x = (x * invN) % q;
  }

  // Negacyclic NTT:
  // pre-twist by psi^i, cyclic NTT with omega=psi^2
  void ntt(Poly &a)
  {
    if (a.N() != N || a.q() != q)
      throw std::runtime_error("NTT params mismatch (toy expects N=8,q=97)");

    // pre-twist: a[i] *= psi^i
    for (std::size_t i = 0; i < N; ++i)
    {
      std::uint64_t t = modexp(psi, i);
      a[i] = (a[i] * t) % q;
    }

    ntt_cyclic(a.a, omega);
  }

  // inverse negacyclic NTT:
  // inverse cyclic NTT, post-twist by psi^{-i}
  void intt(Poly &a)
  {
    if (a.N() != N || a.q() != q)
      throw std::runtime_error("NTT params mismatch (toy expects N=8,q=97)");

    intt_cyclic(a.a, omega);

    // post-twist: a[i] *= psi^{-i}
    std::uint64_t inv_psi = modinv(psi);
    for (std::size_t i = 0; i < N; ++i)
    {
      std::uint64_t t = modexp(inv_psi, i);
      a[i] = (a[i] * t) % q;
    }
  }

  Poly mul_negacyclic_ntt(const Poly &a, const Poly &b)
  {
    Poly A = a;
    Poly B = b;

    ntt(A);
    ntt(B);

    for (std::size_t i = 0; i < N; ++i)
      A[i] = (A[i] * B[i]) % q;

    intt(A);
    return A;
  }
}
