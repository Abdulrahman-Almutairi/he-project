#include "bfv/relin.h"
#include "bfv/ntt.h"
#include <stdexcept>

namespace bfv
{
  static std::uint64_t mod_q_u64(std::int64_t x, std::uint64_t q)
  {
    std::int64_t r = x % (std::int64_t)q;
    if (r < 0) r += (std::int64_t)q;
    return (std::uint64_t)r;
  }

  static Poly poly_scalar_mul(const Params &p, const Poly &a, std::uint64_t k)
  {
    Poly r(p);
    k %= p.q;
    for (std::size_t i = 0; i < p.N; ++i)
      r[i] = (a[i] * k) % p.q;
    return r;
  }

  // Noiseless RLWE encryption of a *known* plaintext poly m:
  // returns (b,a) such that b + a*s = m (mod q).
  static Ciphertext encrypt_poly_noiseless_under_s(const Params &p, const SecretKey &sk,
                                                   const Poly &m, std::mt19937_64 &rng)
  {
    Poly a = sample_uniform_q(p, rng);
    Poly as = mul_negacyclic_ntt(p, a, sk.s);

    Poly b(p);
    for (std::size_t i = 0; i < p.N; ++i)
    {
      // b = m - a*s (mod q)
      b[i] = mod_q_u64((std::int64_t)m[i] - (std::int64_t)as[i], p.q);
    }

    Ciphertext ct(p);
    ct.c0 = b;
    ct.c1 = a;  // by convention: ciphertext is (c0,c1) where c0 + c1*s decrypts
    return ct;
  }

  static std::size_t digits_needed(std::uint64_t q, std::uint64_t base)
  {
    // minimal L s.t. base^L >= q
    std::size_t L = 0;
    std::uint64_t cur = 1;
    while (cur < q)
    {
      // careful: q is tiny for toy; for real we’ll do BigInt
      cur *= base;
      ++L;
    }
    return L;
  }

  RelinKey relin_keygen_noiseless(const Params &p, const SecretKey &sk, std::mt19937_64 &rng,
                                  std::uint64_t base)
  {
    if (base < 2) throw std::runtime_error("relin base must be >= 2");

    RelinKey rk;
    rk.base = base;
    rk.L = digits_needed(p.q, base);
    rk.rk.resize(rk.L, Ciphertext(p));

    // s^2 in the ring
    Poly s2 = mul_negacyclic_ntt(p, sk.s, sk.s);

    // rk[i] encrypts (base^i * s^2)
    std::uint64_t powB = 1 % p.q;
    for (std::size_t i = 0; i < rk.L; ++i)
    {
      Poly m = poly_scalar_mul(p, s2, powB); // m = (B^i)*s^2
      rk.rk[i] = encrypt_poly_noiseless_under_s(p, sk, m, rng);

      powB = (powB * (base % p.q)) % p.q;
    }

    return rk;
  }

  std::vector<Poly> decompose_base_B(const Params &p, const Poly &c2, std::uint64_t base, std::size_t L)
  {
    std::vector<Poly> d;
    d.reserve(L);
    for (std::size_t i = 0; i < L; ++i)
      d.emplace_back(p);

    for (std::size_t j = 0; j < p.N; ++j)
    {
      std::uint64_t x = c2[j] % p.q;
      for (std::size_t i = 0; i < L; ++i)
      {
        d[i][j] = x % base;   // digit in [0, base)
        x /= base;
      }
    }

    return d;
  }

  Ciphertext relinearize(const Params &p, const Ciphertext3 &ct3, const RelinKey &rk)
  {
    if (rk.rk.empty())
      throw std::runtime_error("relinearize: relin key is empty");
    if (rk.L != rk.rk.size())
      throw std::runtime_error("relinearize: bad relin key sizes");

    // Decompose c2
    auto digits = decompose_base_B(p, ct3.c2, rk.base, rk.L);

    // Start with (c0,c1)
    Ciphertext out(p);
    out.c0 = ct3.c0;
    out.c1 = ct3.c1;

    // out += sum_i digits[i] * rk[i]
    // Here digits[i] is a plaintext poly; multiply ciphertext polys by plaintext poly in ring.
    for (std::size_t i = 0; i < rk.L; ++i)
    {
      const Ciphertext &ki = rk.rk[i];
      const Poly &di = digits[i];

      Poly t0 = mul_negacyclic_ntt(p, ki.c0, di);
      Poly t1 = mul_negacyclic_ntt(p, ki.c1, di);

      for (std::size_t j = 0; j < p.N; ++j)
      {
        out.c0[j] = add_mod(out.c0[j], t0[j], p.q);
        out.c1[j] = add_mod(out.c1[j], t1[j], p.q);
      }
    }

    return out;
  }
}
