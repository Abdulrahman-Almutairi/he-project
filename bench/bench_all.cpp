#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "paillier/paillier.h"
#include "bigint/bigint.h"

using Clock = std::chrono::steady_clock;

struct Stats {
  double mean_ms = 0.0;
  double min_ms  = 0.0;
  double max_ms  = 0.0;
};

template <class Fn>
Stats time_trials(Fn&& fn, std::size_t trials) {
  double sum = 0.0;
  double mn  = std::numeric_limits<double>::infinity();
  double mx  = 0.0;

  for (std::size_t i = 0; i < trials; ++i) {
    auto t0 = Clock::now();
    fn();
    auto t1 = Clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    sum += ms;
    if (ms < mn) mn = ms;
    if (ms > mx) mx = ms;
  }

  Stats s;
  s.mean_ms = sum / (double)trials;
  s.min_ms  = mn;
  s.max_ms  = mx;
  return s;
}

static std::uint64_t parse_u64(const std::string& s, std::uint64_t def) {
  try {
    return (std::uint64_t)std::stoull(s);
  } catch (...) {
    return def;
  }
}

static std::size_t parse_sz(const std::string& s, std::size_t def) {
  try {
    return (std::size_t)std::stoull(s);
  } catch (...) {
    return def;
  }
}

static std::string get_arg(int argc, char** argv, const std::string& key, const std::string& def = "") {
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == key && i + 1 < argc) return argv[i + 1];
  }
  return def;
}

static bool has_flag(int argc, char** argv, const std::string& flag) {
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == flag) return true;
  }
  return false;
}

static void emit_csv(const std::string& timestamp,
                     const std::string& commit,
                     const std::string& build,
                     std::uint64_t key_bits,
                     const std::string& op,
                     std::size_t rows,
                     std::size_t cols,
                     std::size_t b_cols,
                     std::size_t trials,
                     const Stats& s) {
  std::cout
    << timestamp << ","
    << commit << ","
    << build << ","
    << key_bits << ","
    << op << ","
    << rows << ","
    << cols << ","
    << b_cols << ","
    << trials << ","
    << s.mean_ms << ","
    << s.min_ms << ","
    << s.max_ms
    << "\n";
}

int main(int argc, char** argv) {
  // CLI args
  const bool csv = has_flag(argc, argv, "--csv");
  const std::string commit = get_arg(argc, argv, "--commit", "no-git");
  const std::string build  = get_arg(argc, argv, "--build",  "unknown");
  const std::string ts     = get_arg(argc, argv, "--timestamp", "unknown");

  std::uint64_t key_bits = parse_u64(get_arg(argc, argv, "--key_bits", "2048"), 2048);
  std::size_t trials     = parse_sz(get_arg(argc, argv, "--trials", "20"), 20);

  std::size_t dot_n   = parse_sz(get_arg(argc, argv, "--dot_n", "256"), 256);
  std::size_t rows    = parse_sz(get_arg(argc, argv, "--rows",  "64"), 64);
  std::size_t cols    = parse_sz(get_arg(argc, argv, "--cols",  "64"), 64);
  std::size_t b_cols  = parse_sz(get_arg(argc, argv, "--b_cols","64"), 64);

  // Keygen (toy version is fast and uses small primes)
  auto kp = he::keygen((unsigned)key_bits);
  std::uint64_t n = bi::to_u64(kp.pk.n);

  // RNG
  std::mt19937_64 gen(12345);
  std::uniform_int_distribution<std::uint64_t> dist_m(0, n ? (n - 1) : 0);
  std::uniform_int_distribution<std::uint64_t> dist_small(0, 100);

  // Prepare some sample plaintexts/ciphertexts
  std::uint64_t m = dist_m(gen);
  auto c = he::encrypt(kp.pk, m);

  auto c1 = he::encrypt(kp.pk, dist_m(gen));
  auto c2 = he::encrypt(kp.pk, dist_m(gen));
  std::uint64_t k = dist_small(gen);

  // dot vectors
  std::vector<bi::BigInt> enc_a(dot_n);
  std::vector<std::uint64_t> b(dot_n);
  for (std::size_t i = 0; i < dot_n; ++i) {
    auto ai = dist_m(gen);
    enc_a[i] = he::encrypt(kp.pk, ai);
    b[i] = dist_small(gen);
  }

  // GEMV: A rows x cols encrypted, x cols plaintext
  std::vector<bi::BigInt> enc_A(rows * cols);
  std::vector<std::uint64_t> x(cols);
  for (std::size_t i = 0; i < rows * cols; ++i) enc_A[i] = he::encrypt(kp.pk, dist_m(gen));
  for (std::size_t j = 0; j < cols; ++j) x[j] = dist_small(gen);

  // GEMM: B is plaintext A_cols x B_cols
  std::vector<std::uint64_t> B(cols * b_cols);
  for (std::size_t i = 0; i < cols * b_cols; ++i) B[i] = dist_small(gen);

  // Benchmarks
  auto st_encrypt = time_trials([&](){
    volatile auto tmp = he::encrypt(kp.pk, dist_m(gen));
    (void)tmp;
  }, trials);

  auto st_decrypt = time_trials([&](){
    volatile auto tmp = he::decrypt_u64(kp.pk, kp.sk, c);
    (void)tmp;
  }, trials);

  auto st_add = time_trials([&](){
    volatile auto tmp = he::add(kp.pk, c1, c2);
    (void)tmp;
  }, trials);

  auto st_scale = time_trials([&](){
    volatile auto tmp = he::scale(kp.pk, c1, k);
    (void)tmp;
  }, trials);

  auto st_dot = time_trials([&](){
    volatile auto tmp = he::enc_dot_cp(kp.pk, enc_a, b);
    (void)tmp;
  }, trials);

  auto st_gemv = time_trials([&](){
    volatile auto tmp = he::enc_gemv_cp(kp.pk, enc_A, rows, cols, x);
    (void)tmp;
  }, trials);

  auto st_gemm = time_trials([&](){
    volatile auto tmp = he::enc_gemm_cp(kp.pk, enc_A, rows, cols, B, b_cols);
    (void)tmp;
  }, trials);

  // Output
  if (!csv) {
    std::cout << "Run with --csv to print CSV lines.\n";
    std::cout << "n=" << n << " trials=" << trials << " rows=" << rows << " cols=" << cols << " b_cols=" << b_cols << " dot_n=" << dot_n << "\n";
    auto pr = [&](const char* name, const Stats& s){
      std::cout << name << ": mean=" << s.mean_ms << "ms min=" << s.min_ms << "ms max=" << s.max_ms << "ms\n";
    };
    pr("encrypt", st_encrypt);
    pr("decrypt", st_decrypt);
    pr("add",     st_add);
    pr("scale",   st_scale);
    pr("dot",     st_dot);
    pr("gemv",    st_gemv);
    pr("gemm",    st_gemm);
    return 0;
  }

  emit_csv(ts, commit, build, key_bits, "encrypt", 0,0,0, trials, st_encrypt);
  emit_csv(ts, commit, build, key_bits, "decrypt", 0,0,0, trials, st_decrypt);
  emit_csv(ts, commit, build, key_bits, "add",     0,0,0, trials, st_add);
  emit_csv(ts, commit, build, key_bits, "scale",   0,0,0, trials, st_scale);
  emit_csv(ts, commit, build, key_bits, "dot_cp",  dot_n,0,0, trials, st_dot);
  emit_csv(ts, commit, build, key_bits, "gemv_cp", rows,cols,0, trials, st_gemv);
  emit_csv(ts, commit, build, key_bits, "gemm_cp", rows,cols,b_cols, trials, st_gemm);

  return 0;
}
