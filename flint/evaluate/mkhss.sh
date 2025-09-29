#! /bin/bash

set -euo pipefail

cd "$(dirname "$0")/.."

# Build the project
cmake --build build

# Benchmark on mkhss::b::* or baseline_mkhss::b::*
./build/benchmark --benchmark_filter='.*mkhss::b::.*' --benchmark_out=evaluate/result.json

echo Successfully evaluated MKHSS
