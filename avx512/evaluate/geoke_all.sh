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

./evaluate/location_ours.sh 2 32
./evaluate/location_baseline.sh 2 32
./evaluate/location_ours.sh 3 48
./evaluate/location_baseline.sh 3 48
./evaluate/location_ours.sh 4 64
./evaluate/location_baseline.sh 4 64
