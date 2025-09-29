#! /bin/bash

set -euo pipefail

set -o noclobber

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 L W b Q T"
    exit 1
fi
L="$1"
W="$2"
b="$3"
Q="$4"
T="$5"

# if [[ -f /sys/devices/system/cpu/cpufreq/policy0/scaling_governor ]]; then
#     grep -q "performance" /sys/devices/system/cpu/cpufreq/policy0/scaling_governor || \
#     (
#         echo "WARNING: CPU scaling governor is not set to 'performance'. This may affect benchmark results." >&2
#         sleep 1
#     )
# else
#     # Non Linux systems or systems without cpufreq
#     echo "WARNING: Unable to check CPU scaling governor. Make sure frequency scaling is off." >&2
#     sleep 1
# fi

cd "$(dirname "$0")/.."

# Build the project
# cmake --build build

# Benchmark our ANIKE scheme

mkdir -p evaluate_rns/anike

./build/rns_fuzzy_pake "$L" "$W" "$b" "$T" "$Q" 2>&1 >evaluate_rns/anike/tmp.txt | tee evaluate_rns/anike/rns_fuzzy_pake_${L}_${W}_${b}_${T}_${Q}_stderr.txt

mv evaluate_rns/anike/tmp.txt evaluate_rns/anike/rns_fuzzy_pake_${L}_${W}_${b}_${T}_${Q}_stdout.txt
