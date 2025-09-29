#! /bin/bash

cd "$(dirname -- "${BASH_SOURCE[0]}")" || exit 1

set -euo pipefail

# Build the project
cd ..
cmake --build build
cd "$(dirname -- "${BASH_SOURCE[0]}")" || exit 1

./fuzzy_pake_flint.sh 8 9 5 2 2
./fuzzy_pake_rns.sh 8 9 5 2 2
./fuzzy_pake_flint.sh 10 8 16 1 1
./fuzzy_pake_rns.sh 10 8 16 1 1
./fuzzy_pake_flint.sh 12 10 8 3 3
./fuzzy_pake_rns.sh 12 10 8 3 3

./geoke_flint.sh 2 32
./geoke_rns.sh 2 32
./geoke_flint.sh 3 48
./geoke_rns.sh 3 48
./geoke_flint.sh 4 64
./geoke_rns.sh 4 64

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
