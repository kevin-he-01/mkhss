#include <benchmark/benchmark.h>
#include <gmpxx.h>
#include <iostream>

#include "util.hh"

namespace primegen {
    inline void PrimeGen(benchmark::State& state) {
        mpz_class p, q;
        for (auto _ : state) {
            // Generate two safe primes of the specified bit length
            // The bit length is given by state.range(0)
            // This will be used to benchmark the prime generation time
            generate_safe_prime(p.get_mpz_t(), state.range(0) / 2);
            generate_safe_prime(q.get_mpz_t(), state.range(0) / 2);
            benchmark::DoNotOptimize(p);
            benchmark::DoNotOptimize(q);
        }
    }
}

BENCHMARK(primegen::PrimeGen)
    ->Arg(3072)->Arg(3520)->Arg(3968); // Example argument for 3072-bit primes

