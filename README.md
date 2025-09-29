# Concretely-Efficient Multi-Key Homomorphic Secret Sharing and Applications

This is the code associated with the paper "Concretely-Efficient Multi-Key Homomorphic Secret Sharing and Applications", to appear at IEEE S&P 2026.

Required:

- Libraries: FLINT, GMP (with C++ support), Google Benchmark Tool, OpenSSL
- Build tool: CMake
- Hardware: For `avx512/`, processor must be x86_64 and support AVX512-IFMA.

Organization:

- `flint/`: A version that uses FLINT for modular exponentiations.
  - Tested on an Ubuntu 22.04 machine on a AMD Ryzen 7 7735U processor equipped with 32 GB of memory (called _laptop_ in the paper).
- `avx512/`: A version that uses Langowski and Devadas for modular exponentiations, which _requires_ a CPU with support of the AVX512-IFMA instruction set.
  - Tested on an AWS bare metal instance of type `c7i.metal-24xl`, an x86_64 machine running Amazon Linux 2023 with support for the AVX512-IFMA instruction set.

## Build instructions

We only tested on the Clang compiler, so if you have multiple compilers installed (e.g., GCC), first do:

```bash
export CC=clang
export CXX=clang++
```

Then, to build:

```bash
# To build FLINT
cd flint/
mkdir build
cmake -S . -B build
cmake --build build
# To build Langowski and Devadas (if your machine supports AVX512-IFMA)
cd ../avx512/
mkdir build
cmake -S . -B build
cmake --build build
```

0. (optional) The choice of the fixed-base exponentiation algorithm (comb vs windowing) depends on the relative speeds of multiplication and squarings, so you may want to tune these hardware-specific values. See section **Tuning instructions** below.
1. To run the benchmarks we used to produce the comparison with the Couteau et al. baseline, run:

    ```bash
    ./flint/evaluate/mkhss.sh # For MKHSS
    ./flint/evaluate/geoke_all.sh # For Geolocation-based key exchange
    ./flint/evaluate/pake_all.sh # For fuzzy PAKE
    ```

    The commands will generate benchmark data in `flint/evaluate/result.json` and `flint/evaluate/anike/`.

2. To run the benchmarks comparing FLINT versus Langowski and Devadas, run:

    ```bash
    ./avx512/evaluate_rns/all.sh
    ```

    The commands will generate benchmark data in `avx512/evaluate_rns/anike`.

3. To view the results in table format (LaTeX), run:

    ```bash
    # Our work vs. Couteau et al.
    python3 flint/evaluate/mkhss_table_runtime.py
    python3 flint/evaluate/fpake_table_runtime.py
    python3 flint/evaluate/geoke_table_runtime.py
    # FLINT vs. AVX512
    python3 avx512/evaluate/rns_fpake_table_runtime.py
    python3 avx512/evaluate/rns_geoke_table_runtime.py
    ```

## Tuning instructions

The choice of the fixed-base exponentiation algorithm (comb vs windowing) depends on the relative speeds of multiplication and squarings, so you may want to tune these hardware-specific values. To do so, run

```bash
# For FLINT
./flint/tune/benchmark.sh
# For Langowski and Devadas
./avx512/tune/benchmark.sh
./avx512/tune/rns_benchmark.sh
```

Then, run these commands and manually paste the data into the initializer list for the variable `precomputed_tune_table` in various files:

```bash
python3 flint/tune/tune.py # Paste output to flint/include/fixed_base_exp.hh
python3 avx512/tune/tune.py flint_bench_out.json # Paste output to avx512/include/fixed_base_exp.hh
python3 avx512/tune/tune.py rns_bench_out.json # Paste output to avx512/include/rns_fbe.hh
```
