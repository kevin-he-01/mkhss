//////// Expected usage ////////

// To define the benchmark for the baseline (See benchmark/baseline_mkhss.cc):

// #define BASELINE
// #include "benchmark/mkhss_template.cc"

// To define the benchmark for the optimized MKHSS (See benchmark/mkhss.cc):

// #include "benchmark/mkhss_template.cc"

//////// This file is to reduce code duplication between the two benchmarks ////////

#include <benchmark/benchmark.h>
#include <gmpxx.h>
#include <iostream>
#include <optional>
#include <unordered_map>

#include "util.hh"
#include "flint/fmpz.h"
#include "flint/fmpz_mod.h"
#include "benchmark/params.hh"

#ifdef BASELINE

#include "baseline_mkhss.hh"
#define mkhss baseline_mkhss

#else

#include "mkhss.hh"

#endif

namespace mkhss::b {
    using namespace mkhss;

    struct config {
        int w;
        int n_width;
    };

    std::unordered_map<int, config> moduli_table_ours = {
        { 1024, {1, 3968} }, // For log_b = 1024, w = 1 and n_width = 3968
        { 2048, {2, 3520} },  // For log_b = 2048, w = 2 and n_width = 3520
        { 4096, {4, 3296} }  // For log_b = 4096, w = 4 and n_width = 3296
    };
    
    std::unordered_map<int, config> moduli_table_baseline = {
        { 4096, {4, 4224} } // For log_b = 4096, w = 4 and n_width = 4224
    };

    // inline int get_w(int log_b) {
    //     if constexpr (mkhss::IS_BASELINE) {
    //         return 3;
    //     } else {
    //         if (log_b == 1024) {
    //             return 2; // For log_b = 1024, use w = 2
    //         } else if (log_b == 2048) {
    //             return 3; // For log_b = 2048, use w = 3
    //         } else {
    //             return 1;
    //         }
    //     }
    // }

    inline void benchmark_create_parameters(parameters &params,
                                            int log2_b) {
        std::unordered_map<int, config> moduli_table = mkhss::IS_BASELINE ? moduli_table_baseline : moduli_table_ours;
        config cfg;
        if (moduli_table.find(log2_b) != moduli_table.end()) {
            cfg = moduli_table[log2_b];
        } else {
            assert(log2_b <= (mkhss::IS_BASELINE ? 2048 : 512));
            // Default config
            cfg.w = mkhss::IS_BASELINE ? 3 : 1; // Default w for baseline is 3, for optimized MKHSS is 1
            cfg.n_width = 3072; // For 128 bit of security against factoring
        }
        // int w = get_w(log2_b);
        create_parameters(params, 
                        SHORT_EXPONENT_ASSUMPTION,
                        STAT_SEC_PARAM,
                        SEC_PARAM,
                        cfg.n_width,
                        log2_b,
                        cfg.w);
    }

    inline void benchmark_setup_crs(crs &crs, const parameters &params) {
        // Use precomputed moduli for the given n_width
        if (moduli_table.find(params.n_width) == moduli_table.end()) {
            std::cerr << "Error: No modulus found for n_width = " << params.n_width << std::endl;
            throw std::runtime_error("Modulus not found");
        }
        setup(crs, params, moduli_table[params.n_width].get_mpz_t());
    }

    static void Setup(benchmark::State& state) {
        initialize_moduli_table();

        parameters params;
        crs crs;

        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        for (auto _ : state) {
            benchmark_setup_crs(crs, params);
            benchmark::DoNotOptimize(crs);
        }
    }

    static void Keygen(benchmark::State& state) {
        initialize_moduli_table();

        parameters params;
        private_key sk;
        public_key pk;
        crs crs;

        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        benchmark_setup_crs(crs, params);
        for (auto _ : state) {
            keygen(crs, pk, sk);
            benchmark::DoNotOptimize(pk);
            benchmark::DoNotOptimize(sk);
        }
    }

    constexpr int x_int = 1; // The integer to be shared

    static void Share(benchmark::State& state) {
        initialize_moduli_table();

        parameters params;
        private_key sk;
        public_key pk;
        private_share x_self;
        public_share x_other;
        FMPZ x;
        crs crs;

        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        benchmark_setup_crs(crs, params);
        keygen(crs, pk, sk);
        gen_random_fmpz(x, log_b);
        for (auto _ : state) {
            share(crs, pk, x, x_self, x_other);
            benchmark::DoNotOptimize(x_self);
            benchmark::DoNotOptimize(x_other);
        }
    }

    // static void ShareAmortized(benchmark::State& state) {
    //     parameters params;
    //     private_key sk;
    //     public_key pk;
    //     private_share x_self;
    //     public_share x_other;
    //     FMPZ x;
    //     crs crs;

    //     int log_b = state.range(0);
    //     int num_shares = state.range(1); // Number of shares to precompute
    //     benchmark_create_parameters(params, log_b);
    //     benchmark_setup_crs(crs, params);
    //     keygen(crs, pk, sk);
    //     gen_random_fmpz(x, log_b);
    //     for (auto _ : state) {
    //         precompute_share(crs, pk, num_shares); // Precompute for a specific number of shares
    //         for (int i = 0; i < num_shares; ++i) {
    //             // Generate a new share for each iteration
    //             share(crs, pk, x, x_self, x_other);
    //         }
    //         benchmark::DoNotOptimize(x_self);
    //         benchmark::DoNotOptimize(x_other);
    //     }
    // }

    static void SyncShareSelf(benchmark::State& state) {
        initialize_moduli_table();

        parameters params;
        private_key sk_alice;
        public_key pk_alice;
        private_key sk_bob;
        public_key pk_bob;
        private_share x_self;
        public_share x_other;
        input_share x_synced;
        FMPZ x(x_int);
        crs crs;
        rms_context ctx;

        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        benchmark_setup_crs(crs, params);

        keygen(crs, pk_alice, sk_alice);
        keygen(crs, pk_bob, sk_bob);
        // gen_random_fmpz(x, test_params.log2_b);
        share(crs, pk_alice, x, x_self, x_other);
        rms_bootstrap(crs, party_id::ALICE, sk_alice, pk_bob, ctx);
        for (auto _ : state) {
            sync_share_self(crs, ctx, sk_alice, pk_bob, x_self, x_synced);
            benchmark::DoNotOptimize(x_synced);
        }
    }

    static void SyncShareOther(benchmark::State& state) {
        initialize_moduli_table();

        parameters params;
        private_key sk_alice;
        public_key pk_alice;
        private_key sk_bob;
        public_key pk_bob;
        private_share x_self;
        public_share x_other;
        input_share x_synced;
        FMPZ x(x_int);
        crs crs;

        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        benchmark_setup_crs(crs, params);

        keygen(crs, pk_alice, sk_alice);
        keygen(crs, pk_bob, sk_bob);
        // gen_random_fmpz(x, test_params.log2_b);
        share(crs, pk_alice, x, x_self, x_other);
        for (auto _ : state) {
            sync_share_other(crs, sk_bob, pk_alice, x_other, x_synced);
            benchmark::DoNotOptimize(x_synced);
        }
    }

    static void RMSBootstrapAlice(benchmark::State& state) {
        initialize_moduli_table();

        parameters params;
        private_key sk_0;
        public_key pk_0;
        private_key sk_1;
        public_key pk_1;
        private_share x_self;
        public_share x_other;
        input_share x_synced;
        FMPZ x(x_int);
        crs crs;
        rms_context ctx;

        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        benchmark_setup_crs(crs, params);

        keygen(crs, pk_0, sk_0);
        keygen(crs, pk_1, sk_1);
        // gen_random_fmpz(x, test_params.log2_b);
        share(crs, pk_0, x, x_self, x_other);
        sync_share_other(crs, sk_1, pk_0, x_other, x_synced);
        for (auto _ : state) {
            rms_bootstrap(crs, party_id::ALICE, sk_1, pk_0, ctx);
            benchmark::DoNotOptimize(ctx);
        }
    }

    static void RMSBootstrapBob(benchmark::State& state) {
        initialize_moduli_table();

        parameters params;
        private_key sk_0;
        public_key pk_0;
        private_key sk_1;
        public_key pk_1;
        private_share x_self;
        public_share x_other;
        input_share x_synced;
        FMPZ x(x_int);
        crs crs;
        rms_context ctx;

        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        benchmark_setup_crs(crs, params);

        keygen(crs, pk_0, sk_0);
        keygen(crs, pk_1, sk_1);
        // gen_random_fmpz(x, test_params.log2_b);
        share(crs, pk_0, x, x_self, x_other);
        sync_share_other(crs, sk_1, pk_0, x_other, x_synced);
        for (auto _ : state) {
            rms_bootstrap(crs, party_id::BOB, sk_1, pk_0, ctx);
            benchmark::DoNotOptimize(ctx);
        }
    }

    static void RMSISub(benchmark::State& state) {
        initialize_moduli_table();

        parameters params;
        private_key sk_alice;
        public_key pk_alice;
        private_key sk_bob;
        public_key pk_bob;

        private_share x_self;
        public_share x_other;
        input_share x_synced;

        private_share y_self;
        public_share y_other;
        input_share y_synced;

        input_share x_minus_y;
        // input_share x_plus_y;

        FMPZ x(x_int), y(x_int);
        crs crs;
        // memory_share one, x_mem, x2_mem;
    
        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        benchmark_setup_crs(crs, params);
    
        keygen(crs, pk_alice, sk_alice);
        keygen(crs, pk_bob, sk_bob);
        // gen_random_fmpz(x, test_params.log2_b);
        share(crs, pk_alice, x, x_self, x_other);
        sync_share_other(crs, sk_bob, pk_alice, x_other, x_synced);

        share(crs, pk_bob, y, y_self, y_other);
        sync_share_other(crs, sk_alice, pk_bob, y_other, y_synced);
        // rms_bootstrap(crs, sk_bob, pk_alice, one);
    
        // RMSContext ctx(crs);
        // ctx.mult(0xFEEDFACEULL, x_synced, one, x_mem);
    
        for (auto _ : state) {
            rms_isub(crs, x_synced, y_synced, x_minus_y);
            // ctx.mult(0xFEEDF00DULL, x_synced, x_mem, x2_mem);
            // benchmark::DoNotOptimize(x2_mem);
        }
    }

    static void RMSIAdd(benchmark::State& state) {
        initialize_moduli_table();

        parameters params;
        private_key sk_alice;
        public_key pk_alice;
        private_key sk_bob;
        public_key pk_bob;

        private_share x_self;
        public_share x_other;
        input_share x_synced;

        private_share y_self;
        public_share y_other;
        input_share y_synced;

        // input_share x_minus_y;
        input_share x_plus_y;

        FMPZ x(x_int), y(x_int);
        crs crs;
        // memory_share one, x_mem, x2_mem;
    
        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        benchmark_setup_crs(crs, params);
    
        keygen(crs, pk_alice, sk_alice);
        keygen(crs, pk_bob, sk_bob);
        // gen_random_fmpz(x, test_params.log2_b);
        share(crs, pk_alice, x, x_self, x_other);
        sync_share_other(crs, sk_bob, pk_alice, x_other, x_synced);

        share(crs, pk_bob, y, y_self, y_other);
        sync_share_other(crs, sk_alice, pk_bob, y_other, y_synced);
        // rms_bootstrap(crs, sk_bob, pk_alice, one);
    
        // RMSContext ctx(crs);
        // ctx.mult(0xFEEDFACEULL, x_synced, one, x_mem);
    
        for (auto _ : state) {
            rms_iadd(crs, x_synced, y_synced, x_plus_y);
            benchmark::DoNotOptimize(x_plus_y);
        }
    }

    static void RMSMAdd(benchmark::State& state) {
        initialize_moduli_table();

        parameters params;
        private_key sk_alice;
        public_key pk_alice;
        private_key sk_bob;
        public_key pk_bob;

        private_share x_self;
        public_share x_other;
        input_share x_synced;

        private_share y_self;
        public_share y_other;
        input_share y_synced;

        FMPZ x(x_int), y(x_int);
        crs crs;
        rms_context ctx;
        memory_share x_mem, y_mem, x_plus_y_mem;
    
        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        benchmark_setup_crs(crs, params);
    
        keygen(crs, pk_alice, sk_alice);
        keygen(crs, pk_bob, sk_bob);
        // gen_random_fmpz(x, test_params.log2_b);
        share(crs, pk_alice, x, x_self, x_other);
        sync_share_other(crs, sk_bob, pk_alice, x_other, x_synced);

        share(crs, pk_bob, y, y_self, y_other);
        sync_share_other(crs, sk_alice, pk_bob, y_other, y_synced);
        rms_bootstrap(crs, party_id::BOB, sk_bob, pk_alice, ctx);
    
        ctx.convert(crs, 0xFEEDFACEULL, x_synced, x_mem);
        ctx.convert(crs, 0xFEEDF00DULL, y_synced, y_mem);
    
        for (auto _ : state) {
            rms_madd(x_mem, y_mem, x_plus_y_mem);
            benchmark::DoNotOptimize(x_plus_y_mem);
        }
    }

    static void RMSMult(benchmark::State& state) {
        initialize_moduli_table();

        parameters params;
        private_key sk_alice;
        public_key pk_alice;
        private_key sk_bob;
        public_key pk_bob;
        private_share x_self;
        public_share x_other;
        input_share x_synced;
        FMPZ x(x_int);
        crs crs;
        rms_context ctx;
        memory_share x_mem, x2_mem;
    
        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        benchmark_setup_crs(crs, params);
    
        keygen(crs, pk_alice, sk_alice);
        keygen(crs, pk_bob, sk_bob);
        // gen_random_fmpz(x, test_params.log2_b);
        share(crs, pk_alice, x, x_self, x_other);
        sync_share_other(crs, sk_bob, pk_alice, x_other, x_synced);
        rms_bootstrap(crs, party_id::BOB, sk_bob, pk_alice, ctx);
    
        ctx.convert(crs, 0xFEEDFACEULL, x_synced, x_mem);
    
        for (auto _ : state) {
            ctx.mult(crs, 0xFEEDF00DULL, x_synced, x_mem, x2_mem);
            benchmark::DoNotOptimize(x2_mem);
        }
    }

    static void RMSConvert(benchmark::State& state) {
        initialize_moduli_table();
        
        // Convert an input share to a memory share
        // This is usually slightly cheaper than a multiplication, because the initial subtractive share of 1 has only one bit:
        // <1>_A = 0 and <1>_B = 1
        // and exponentiation by 0 or 1 is as fast as multiplication
        parameters params;
        private_key sk_alice;
        public_key pk_alice;
        private_key sk_bob;
        public_key pk_bob;
        private_share x_self;
        public_share x_other;
        input_share x_synced;
        FMPZ x(x_int);
        crs crs;
        rms_context ctx;
        memory_share one, x_mem;
    
        int log_b = state.range(0);
        benchmark_create_parameters(params, log_b);
        benchmark_setup_crs(crs, params);
    
        keygen(crs, pk_alice, sk_alice);
        keygen(crs, pk_bob, sk_bob);
        // gen_random_fmpz(x, test_params.log2_b);
        share(crs, pk_alice, x, x_self, x_other);
        sync_share_other(crs, sk_bob, pk_alice, x_other, x_synced);
        rms_bootstrap(crs, party_id::BOB, sk_bob, pk_alice, ctx);
    
        for (auto _ : state) {
            ctx.convert(crs, 0xFEEDF00DULL, x_synced, x_mem);
            benchmark::DoNotOptimize(x_mem);
        }
    }
}

#ifdef BASELINE
#define LOG_B_RANGE {1, 4096}
#else
#define LOG_B_RANGE {1, 64, 128, 256, 512, 1024, 2048, 4096}
#endif

#define BENCHMARK_MKHSS(FuncName) BENCHMARK(FuncName)->ArgsProduct({LOG_B_RANGE})
#define BENCHMARK_MKHSS_RAW(FuncName) BENCHMARK(FuncName)

BENCHMARK_MKHSS(mkhss::b::Setup);
BENCHMARK_MKHSS(mkhss::b::Keygen);
BENCHMARK_MKHSS(mkhss::b::Share);
// BENCHMARK_MKHSS_RAW(mkhss::b::ShareAmortized)->ArgsProduct({LOG_B_RANGE, benchmark::CreateRange(2, 512, 2)}); // Amortized sharing with different number of shares
BENCHMARK_MKHSS(mkhss::b::RMSBootstrapAlice);
BENCHMARK_MKHSS(mkhss::b::RMSBootstrapBob);
BENCHMARK_MKHSS(mkhss::b::SyncShareSelf);
BENCHMARK_MKHSS(mkhss::b::SyncShareOther);
BENCHMARK_MKHSS(mkhss::b::RMSIAdd);
BENCHMARK_MKHSS(mkhss::b::RMSISub);
BENCHMARK_MKHSS(mkhss::b::RMSMAdd);
BENCHMARK_MKHSS(mkhss::b::RMSMult);
BENCHMARK_MKHSS(mkhss::b::RMSConvert);

#undef LOG_B_RANGE
