#!/usr/bin/env bash
set -euo pipefail

BUILD=${1:-Release}

TS=$(date -Iseconds)
COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "no-git")

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=$BUILD
cmake --build build -j

./build/bench/bench_all \
  --csv \
  --commit "$COMMIT" \
  --build "$BUILD" \
  --timestamp "$TS" \
  --trials 20 \
  --dot_n 256 \
  --rows 64 --cols 64 --b_cols 64 \
  >> bench_results/results.csv

echo "Saved results to bench_results/results.csv"
