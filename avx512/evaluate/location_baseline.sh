#! /bin/bash

set -euo pipefail

# set -o noclobber

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 D L"
    exit 1
fi
D="$1"
L="$2"

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

# Benchmark our ANIKE scheme

mkdir -p evaluate/anike

./build/baseline_location_anike "$L" "$D" 2>&1 >evaluate/anike/tmp.txt | tee evaluate/anike/baseline_location_${L}_${D}_stderr.txt

mv evaluate/anike/tmp.txt evaluate/anike/baseline_location_${L}_${D}_stdout.txt
