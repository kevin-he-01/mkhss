#include <benchmark/benchmark.h>
#include "nike.hh"

namespace nike_bench {
    using namespace nike;

    static void KeyGen(benchmark::State& state) {
        crs crs;
        setup(crs);
        private_key sk;
        public_key pk;
        for (auto _ : state) {
            keygen(crs, sk, pk);
            benchmark::DoNotOptimize(sk);
            benchmark::DoNotOptimize(pk);
        }
    }

    static void KeyDer(benchmark::State& state) {
        crs crs;
        setup(crs);
        private_key sk, sk_B;
        public_key pk, pk_B;
        uint8_t shared_key[SHARED_KEY_LEN];
        keygen(crs, sk, pk);
        keygen(crs, sk_B, pk_B);
        for (auto _ : state) {
            keyder(crs, sk, pk_B, shared_key);
            benchmark::DoNotOptimize(shared_key);
        }
    }
}

BENCHMARK(nike_bench::KeyGen);
BENCHMARK(nike_bench::KeyDer);
