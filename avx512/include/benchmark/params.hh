#ifndef BENCHMARK_PARAMS_HH
#define BENCHMARK_PARAMS_HH

#include <benchmark/benchmark.h>
#include <gmpxx.h>
#include <iostream>
#include <unordered_map>

#include "util.hh"
#include "mkhss.hh"
#include "flint/fmpz.h"
#include "flint/fmpz_mod.h"

// Parameters
constexpr bool SHORT_EXPONENT_ASSUMPTION = true;
constexpr int SEC_PARAM = 128; // Lambda. Computational security parameter (against offline attacks)
constexpr int STAT_SEC_PARAM = TAU; // Tau. Statistical security parameter (for correctness)

constexpr int LOG2_B_EVALUATE_MKHSS = 1;
constexpr int STRIDE = 200;

constexpr int N_WIDTH = 3072; // Needed for 128 bits of security against factoring

extern mpz_class test_n;
extern mpz_class test_n_squared;

// Pregenerate moduli for optimal MKHSS performance
extern std::unordered_map<int, mpz_class> moduli_table;

void print_benchmark_parameters();
void initialize_moduli_table();

inline void get_modulus(int n_width, mpz_class &modulus) {
    auto it = moduli_table.find(n_width);
    if (it != moduli_table.end()) {
        modulus = it->second;
    } else {
        std::cerr << "Error: No modulus found for n_width = " << n_width << std::endl;
        throw std::runtime_error("Modulus not found");
    }
}

#endif // BENCHMARK_PARAMS_HH
