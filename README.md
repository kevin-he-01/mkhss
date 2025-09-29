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
