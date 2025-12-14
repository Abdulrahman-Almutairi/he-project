how to benchmark:
1- git commit -m "YOUR_CHANGES"
2- ./tools/run_bench.sh Release
3- comment insights here.



03949bc:
- first implementation of GEMV, GEMM operations.
- we can consider this a baseline.
- Very Basic.


3dbf6f2:
- removing row allocations mattered in GEMV. faster.
- GEMV is memory/allocation sensative.
- Did not affect GEMM however.

0a23706:
- adding blocking to GEMM implementation (Tiling).
- minor speedup.
- best config so far (16, 8, 16) -> (rows tile, cols tile, inner (k-dimension) tile)

da9e0fc:
- added openMp to GEMM, huge speed up.
- OMP num of threads for benchmarking: 8
- Pragma omp parallel loop collapse(2)
