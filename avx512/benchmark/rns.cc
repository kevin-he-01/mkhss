#include <benchmark/benchmark.h>
#include <optional>

#include "util.hh"
#include "rns_helper.hh"
#include "params.hh"
#include "rns_modexp.hh"
#include "rns_fbe.hh"

namespace rns {
    typedef MontInt<N2_LIMBS> MontIntN2;
    typedef RNSCtx<N2_LIMBS> RNSCtxN2;
    typedef std::optional<RNSCtxN2> RNSObj;

    mpz_class test_n(TEST_N_HEX);

    static void get_rns_obj(RNSObj &ctx, FMPZ &modulus) {
        modulus.set_mpz(test_n.get_mpz_t());
        modulus.mul(modulus, modulus);
        PrecomputeIO io(N2_PRECOMPUTED);
        ctx.emplace(io, modulus);
    }

    static void ToMont(benchmark::State &state) {
        RNSObj ctx;
        FMPZ n2;
        get_rns_obj(ctx, n2);
        MontIntN2 a_mont;
        FMPZ a;
        gen_random_under_fmpz(a, n2);
        for (auto _ : state) {
            benchmark::DoNotOptimize(a);
            ctx->to_mont(a_mont, a);
            benchmark::DoNotOptimize(a_mont);
        }
    }

    static void ToMontAVX(benchmark::State &state) {
        RNSObj ctx;
        FMPZ n2;
        get_rns_obj(ctx, n2);
        MontIntN2 a_mont;
        FMPZ a;
        gen_random_under_fmpz(a, n2);
        for (auto _ : state) {
            benchmark::DoNotOptimize(a);
            ctx->to_mont_avx(a_mont, a);
            benchmark::DoNotOptimize(a_mont);
        }
    }

    static void MulMont(benchmark::State &state) {
        RNSObj ctx;
        FMPZ n2;
        get_rns_obj(ctx, n2);
        MontIntN2 a_mont, b_mont, c_mont;
        FMPZ a, b;
        gen_random_under_fmpz(a, n2);
        gen_random_under_fmpz(b, n2);
        ctx->to_mont(a_mont, a);
        ctx->to_mont(b_mont, b);
        for (auto _ : state) {
            benchmark::DoNotOptimize(a_mont);
            benchmark::DoNotOptimize(b_mont);
            ctx->mul_mont(c_mont, a_mont, b_mont);
            benchmark::DoNotOptimize(c_mont);
        }
    }

    static void MulMontSeq(benchmark::State &state) {
        RNSObj ctx;
        FMPZ n2;
        get_rns_obj(ctx, n2);
        MontIntN2 a_mont, b_mont;
        FMPZ a, b;
        gen_random_under_fmpz(a, n2);
        gen_random_under_fmpz(b, n2);
        ctx->to_mont(a_mont, a);
        ctx->to_mont(b_mont, b);
        benchmark::DoNotOptimize(a_mont);
        for (auto _ : state) {
            benchmark::DoNotOptimize(b_mont);
            ctx->mul_mont(a_mont, a_mont, b_mont);
            benchmark::DoNotOptimize(a_mont);
        }
    }

    static void FromMont(benchmark::State &state) {
        RNSObj ctx;
        FMPZ n2;
        get_rns_obj(ctx, n2);
        MontIntN2 a_mont;
        FMPZ a;
        gen_random_under_fmpz(a, n2);
        ctx->to_mont(a_mont, a);
        FMPZ result;
        for (auto _ : state) {
            benchmark::DoNotOptimize(a_mont);
            ctx->from_mont(result, a_mont);
            benchmark::DoNotOptimize(result);
        }
    }

    static void FromMontAVX(benchmark::State &state) {
        RNSObj ctx;
        FMPZ n2;
        get_rns_obj(ctx, n2);
        MontIntN2 a_mont;
        FMPZ a;
        gen_random_under_fmpz(a, n2);
        ctx->to_mont(a_mont, a);
        FMPZ result;
        for (auto _ : state) {
            benchmark::DoNotOptimize(a_mont);
            ctx->from_mont_avx(result, a_mont);
            benchmark::DoNotOptimize(result);
        }
    }
}

namespace rns::modexp {
    static void ModExpFmpz(benchmark::State& state) {
        // FMPZ base;
        // FMPZ exp;
        // FMPZ result;
        // FMPZ modulus = test_n;
        // modulus.pow(modulus, state.range(1) + 1); // Set modulus to n^{w+1}

        // // Generate random base and exponent
        // gen_random_under_fmpz(base, modulus);
        // gen_random_fmpz(exp, state.range(0)); // Random exponent of given length
        // for (auto _ : state) {
        //     benchmark::DoNotOptimize(base);
        //     benchmark::DoNotOptimize(exp);
        //     // Perform modular exponentiation
        //     result.powm(base, exp, modulus);
        //     benchmark::DoNotOptimize(result);
        // }
        RNSObj ctx;
        FMPZ base;
        FMPZ exp;
        MontIntN2 result;
        FMPZ modulus;
        get_rns_obj(ctx, modulus);
        
        gen_random_under_fmpz(base, modulus);
        gen_random_fmpz(exp, state.range(0)); // Random exponent of given length
        for (auto _ : state) {
            benchmark::DoNotOptimize(base);
            benchmark::DoNotOptimize(exp);
            rns_modexp::sliding_window(result, base, exp, modulus, ctx);
            benchmark::DoNotOptimize(result);
        }
    }

    static void ModExpFixedBasePrecompute(benchmark::State& state) {
        RNSObj ctx;
        FMPZ base;
        FMPZ exp;
        FMPZ result;
        FMPZ modulus;
        get_rns_obj(ctx, modulus);

        // Generate random base and exponent
        gen_random_under_fmpz(base, modulus);
        gen_random_fmpz(exp, state.range(0)); // Random exponent of given length
        for (auto _ : state) {
            // Initialize fixed base exponentiation context 
            rns_fbe::FixedBaseExp fixed_base_exp;
            fixed_base_exp.precompute(base, modulus, state.range(0), state.range(1), ctx); // Precompute powers of base
            benchmark::DoNotOptimize(fixed_base_exp);
        }
    }

    static void ModExpFixedBase(benchmark::State& state) {
        RNSObj ctx;
        FMPZ base;
        FMPZ exp;
        FMPZ result;
        FMPZ modulus;
        get_rns_obj(ctx, modulus);

        // Generate random base and exponent
        gen_random_under_fmpz(base, modulus);
        gen_random_fmpz(exp, state.range(0)); // Random exponent of given length
        // Initialize fixed base exponentiation context 
        rns_fbe::FixedBaseExp fixed_base_exp;
        fixed_base_exp.precompute(base, modulus, state.range(0), state.range(1), ctx); // Precompute powers of base
        for (auto _ : state) {
            benchmark::DoNotOptimize(base);
            benchmark::DoNotOptimize(exp);
            // Perform modular exponentiation
            // result.powm(base, exp, modulus);
            fixed_base_exp.compute(exp, result); // Use precomputed powers for exponentiation
            benchmark::DoNotOptimize(result);
        }
    }
}

static void RNSFixedBaseArgs_i(benchmark::internal::Benchmark* b, int64_t i) {
  b->Args({129, i})->Args({256, i})->Args({899, i});
}

static void RNSFixedBaseArgs(benchmark::internal::Benchmark* b) {
  RNSFixedBaseArgs_i(b, -1);
  for (int64_t i = 1; i <= 9; ++i) {
    RNSFixedBaseArgs_i(b, i);
  }
}

BENCHMARK(rns::ToMont);
BENCHMARK(rns::ToMontAVX);
BENCHMARK(rns::MulMont);
BENCHMARK(rns::MulMontSeq);
BENCHMARK(rns::FromMont);
BENCHMARK(rns::FromMontAVX);

BENCHMARK(rns::modexp::ModExpFmpz)->Arg(129)->Arg(256)->Arg(899);
BENCHMARK(rns::modexp::ModExpFixedBasePrecompute)->Apply(RNSFixedBaseArgs);
BENCHMARK(rns::modexp::ModExpFixedBase)->Apply(RNSFixedBaseArgs);
