#! /bin/bash

set -euo pipefail
if [[ -f /sys/devices/system/cpu/cpufreq/policy0/scaling_governor ]]; then
    grep -q "performance" /sys/devices/system/cpu/cpufreq/policy0/scaling_governor || \
    (
        echo "WARNING: CPU scaling governor is not set to 'performance'. This may affect benchmark results." >&2
        sleep 1
    )
else
    # Non Linux systems or systems without cpufreq
    echo "WARNING: Unable to check CPU scaling governor. Make sure frequency scaling is off." >&2
    sleep 1
fi
cd "$(dirname "$0")/.."
# Build the project
cmake --build build

set +e

./evaluate/fuzzy_pake_ours.sh 8 9 5 2 2
./evaluate/fuzzy_pake_baseline.sh 8 9 5 2 2
./evaluate/fuzzy_pake_ours.sh 10 8 16 1 1
./evaluate/fuzzy_pake_baseline.sh 10 8 16 1 1
./evaluate/fuzzy_pake_ours.sh 12 10 8 3 3
./evaluate/fuzzy_pake_baseline.sh 12 10 8 3 3
