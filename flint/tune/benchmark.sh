#! /bin/bash
# Run benchmarks to determine the optimal parameters for tuning fixed-based exponentiation.

set -e

cd "$(dirname "$0")/.."

# Build the project
cmake --build build

./build/benchmark --benchmark_filter='^fmpz_bench::modexp::' --benchmark_repetitions=5 --benchmark_out=tune/bench_out.json

echo Generated benchmark data for tuning fixed-based exponentiation.
