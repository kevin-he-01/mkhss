#! /bin/bash

cd "$(dirname "$0")/.."
set -euxo pipefail

COUNTERS='task-clock,context-switches,cpu-migrations,page-faults,cycles,stalled-cycles-frontend,instructions,branches,branch-misses,cache-references,cache-misses'

sudo perf stat -e "$COUNTERS" -o ./perf_exp/flint_stat.txt -- ./build/fuzzy_pake 8 9 5 2 2 --alice-only
sudo perf stat -e "$COUNTERS" -o ./perf_exp/rns_stat.txt -- ./build/rns_fuzzy_pake 8 9 5 2 2 --alice-only
sudo perf record -o ./perf_exp/perf.data.flint -- ./build/fuzzy_pake 8 9 5 2 2 --alice-only
sudo perf record -o ./perf_exp/perf.data.rns -- ./build/rns_fuzzy_pake 8 9 5 2 2 --alice-only

# Signal done
echo $'\a'
sleep 1
echo $'\a'
sleep 1
echo $'\a'
sleep 1
echo $'\a'
sleep 1
echo $'\a'
