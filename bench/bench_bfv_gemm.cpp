#include "bfv/bfv.h"
#include "bfv/eval.h"
#include "bfv/relin.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

static double ms_since(const std::chrono::steady_clock::time_point &t0,
                       const std::chrono::steady_clock::time_point &t1)
{
  return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

int main()
{
  // Baseline params (works + NTT friendly)
  bfv::Params P{1024, 12289, 2};

  // Matrix sizes (start small, then increase)
  const std::size_t A_rows = 32;
  const std::size_t A_cols = 32;
  const std::size_t B_cols = 32;

  std::mt19937_64 rng(123);

  // ---------- Keygen ----------
  auto t0 = std::chrono::steady_clock::now();
  auto kp = bfv::keygen(P, rng);
  auto rlk = bfv::relin_keygen_noiseless(P, kp.sk, rng, 4);
  auto t1 = std::chrono::steady_clock::now();
  double ms_keygen = ms_since(t0, t1);

  // ---------- Make plaintext matrices ----------
  std::vector<std::uint64_t> A_plain(A_rows * A_cols);
  std::vector<std::uint64_t> B_plain(A_cols * B_cols);

  for (auto &x : A_plain) x = rng() % P.t;
  for (auto &x : B_plain) x = rng() % P.t;

  // ---------- Encrypt ----------
  t0 = std::chrono::steady_clock::now();
  std::vector<bfv::Ciphertext> enc_A;
  std::vector<bfv::Ciphertext> enc_B;
  enc_A.reserve(A_plain.size());
  enc_B.reserve(B_plain.size());

  for (auto x : A_plain) enc_A.push_back(bfv::encrypt_const(P, kp.pk, x, rng));
  for (auto x : B_plain) enc_B.push_back(bfv::encrypt_const(P, kp.pk, x, rng));
  t1 = std::chrono::steady_clock::now();
  double ms_encrypt = ms_since(t0, t1);

  // ---------- GEMM ----------
  t0 = std::chrono::steady_clock::now();
  auto enc_C = bfv::enc_gemm_cc(P, enc_A, A_rows, A_cols, enc_B, B_cols, rlk);
  t1 = std::chrono::steady_clock::now();
  double ms_gemm = ms_since(t0, t1);

  // ---------- Decrypt ----------
  t0 = std::chrono::steady_clock::now();
  std::vector<std::uint64_t> C_dec(enc_C.size());
  for (std::size_t i = 0; i < enc_C.size(); ++i)
    C_dec[i] = bfv::decrypt_const(P, kp.sk, enc_C[i]);
  t1 = std::chrono::steady_clock::now();
  double ms_decrypt = ms_since(t0, t1);

  // ---------- Output CSV ----------
  std::ofstream out("bfv_gemm_bench.csv", std::ios::app);
  out << "N,q,t,A_rows,A_cols,B_cols,ms_keygen,ms_encrypt,ms_gemm,ms_decrypt\n";
  out << P.N << "," << P.q << "," << P.t << ","
      << A_rows << "," << A_cols << "," << B_cols << ","
      << ms_keygen << "," << ms_encrypt << "," << ms_gemm << "," << ms_decrypt << "\n";

  std::cout << "Wrote bfv_gemm_bench.csv\n";
  std::cout << "keygen(ms): " << ms_keygen << "\n";
  std::cout << "encrypt(ms): " << ms_encrypt << "\n";
  std::cout << "gemm(ms): " << ms_gemm << "\n";
  std::cout << "decrypt(ms): " << ms_decrypt << "\n";
  return 0;
}
