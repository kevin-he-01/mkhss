#include <benchmark/benchmark.h>
#include <gmpxx.h>
#include <iostream>

#include "util.hh"
#include "flint/fmpz.h"
#include "flint/fmpz_mod.h"
#include "benchmark/params.hh"
#include "fixed_base_exp.hh"

namespace fmpz_bench {
    // static void ModMulMPZ(benchmark::State& state) {
    //     mpz_class a, b, c;
    //     mpz_class n_w_1 = test_n;
    //     mpz_pow_ui(n_w_1.get_mpz_t(), n_w_1.get_mpz_t(), state.range(0) + 1); // n^{w+1}
    //     gen_random_under(Z(a), Z(n_w_1));
    //     gen_random_under(Z(b), Z(n_w_1));

    //     for (auto _ : state) {
    //         benchmark::DoNotOptimize(a);
    //         benchmark::DoNotOptimize(b);
    //         mpz_mul(c.get_mpz_t(), Z(a), Z(b));
    //         mpz_mod(c.get_mpz_t(), c.get_mpz_t(), Z(n_w_1));
    //         benchmark::DoNotOptimize(c);
    //     }
    // }
    
    // static void MulMPZ(benchmark::State& state) {
    //     mpz_class a, b, c;
    //     mpz_class n_w_1 = test_n;
    //     mpz_pow_ui(n_w_1.get_mpz_t(), n_w_1.get_mpz_t(), state.range(0) + 1); // n^{w+1}
    //     gen_random_under(Z(a), Z(n_w_1));
    //     gen_random_under(Z(b), Z(n_w_1));

    //     for (auto _ : state) {
    //         benchmark::DoNotOptimize(a);
    //         benchmark::DoNotOptimize(b);
    //         mpz_mul(c.get_mpz_t(), Z(a), Z(b));
    //         benchmark::DoNotOptimize(c);
    //     }
    // }
    
    // static void SqrMPZ(benchmark::State& state) {
    //     mpz_class a, b, c;
    //     mpz_class n_w_1 = test_n;
    //     mpz_pow_ui(n_w_1.get_mpz_t(), n_w_1.get_mpz_t(), state.range(0) + 1); // n^{w+1}
    //     gen_random_under(Z(a), Z(n_w_1));
    //     gen_random_under(Z(b), Z(n_w_1));

    //     for (auto _ : state) {
    //         benchmark::DoNotOptimize(a);
    //         mpz_mul(c.get_mpz_t(), Z(a), Z(a));
    //         benchmark::DoNotOptimize(c);
    //     }
    // }
    
    // static void ModMPZ(benchmark::State& state) {
    //     // mpz_class c, n2 = test_n_squared * test_n_squared;
    //     // gen_random_under(Z(c), Z(n2));
    //     mpz_class a, b, c, result;
    //     mpz_class n_w_1 = test_n;
    //     mpz_pow_ui(n_w_1.get_mpz_t(), n_w_1.get_mpz_t(), state.range(0) + 1); // n^{w+1}
    //     gen_random_under(Z(a), Z(n_w_1));
    //     gen_random_under(Z(b), Z(n_w_1));

    //     mpz_mul(c.get_mpz_t(), Z(a), Z(b));
    //     for (auto _ : state) {
    //         benchmark::DoNotOptimize(c);
    //         mpz_mod(result.get_mpz_t(), c.get_mpz_t(), Z(test_n_squared));
    //         benchmark::DoNotOptimize(result);
    //     }
    // }

    static void MulFmpzFlexSize(benchmark::State& state) {
        FMPZ a, b, c;
        // FMPZ modulus;
        gen_random_fmpz(a, state.range(0));
        gen_random_fmpz(b, state.range(1));
    
        for (auto _ : state) {
            benchmark::DoNotOptimize(a);
            benchmark::DoNotOptimize(b);
            c.mul(a, b);
            benchmark::DoNotOptimize(c);
        }
    }
    
    static void MulFmpz(benchmark::State& state) {
        FMPZ a, b, c;
        FMPZ modulus = test_n;
        modulus.pow(modulus, state.range(0) + 1); // Set modulus to n^{w+1}
        gen_random_under_fmpz(a, modulus);
        gen_random_under_fmpz(b, modulus);
    
        for (auto _ : state) {
            benchmark::DoNotOptimize(a);
            benchmark::DoNotOptimize(b);
            c.mul(a, b);
            benchmark::DoNotOptimize(c);
        }
    }

    static void MulSecondFmpz(benchmark::State& state) {
        FMPZ a, b, c;
        FMPZ modulus = test_n;
        modulus.pow(modulus, state.range(0) + 1); // Set modulus to n^{w+1}
        gen_random_under_fmpz(a, modulus);
        gen_random_under_fmpz(b, modulus);
    
        a.mul(a, b);
        for (auto _ : state) {
            benchmark::DoNotOptimize(a);
            benchmark::DoNotOptimize(b);
            c.mul(a, b);
            benchmark::DoNotOptimize(c);
        }
    }
    
    static void SqrFmpz(benchmark::State& state) {
        FMPZ a, c;
        FMPZ modulus = test_n;
        modulus.pow(modulus, state.range(0) + 1); // Set modulus to n^{w+1}
        gen_random_under_fmpz(a, modulus);
    
        for (auto _ : state) {
            benchmark::DoNotOptimize(a);
            c.mul(a, a);
            benchmark::DoNotOptimize(c);
        }
    }

    static void SqrSecondFmpz(benchmark::State& state) {
        FMPZ a, b, c;
        FMPZ modulus = test_n;
        modulus.pow(modulus, state.range(0) + 1); // Set modulus to n^{w+1}
        gen_random_under_fmpz(a, modulus);
    
        b.mul(a, a);
        for (auto _ : state) {
            benchmark::DoNotOptimize(a);
            benchmark::DoNotOptimize(b);
            c.mul(b, a);
            benchmark::DoNotOptimize(c);
        }
    }

    static void PrecomputedModFmpz(benchmark::State& state) {
        FMPZModCtx ctx;
        FMPZ modulus = test_n;
        modulus.pow(modulus, state.range(0) + 1); // Set modulus to n^{w+1}
        FMPZ a, b, c, result;
        gen_random_under_fmpz(a, modulus);
        gen_random_under_fmpz(b, modulus);
        ctx.set_modulus(modulus);
        c.mul(a, b);
        for (auto _ : state) {
            benchmark::DoNotOptimize(c);
            ctx.mod(result, c);
            benchmark::DoNotOptimize(result);
        }
    }

    static void PrecomputedModFmpzFlexMod(benchmark::State& state) {
        FMPZModCtx ctx;
        FMPZ modulus;
        gen_random_fmpz(modulus, state.range(0)); // Generate a random modulus of given length
        FMPZ a, b, c, result;
        gen_random_under_fmpz(a, modulus);
        gen_random_under_fmpz(b, modulus);
        ctx.set_modulus(modulus);
        c.mul(a, b);
        for (auto _ : state) {
            benchmark::DoNotOptimize(c);
            ctx.mod(result, c);
            benchmark::DoNotOptimize(result);
        }
    }

    static void PrecomputedModSecondFmpz(benchmark::State& state) {
        FMPZModCtx ctx;
        FMPZ modulus = test_n;
        modulus.pow(modulus, state.range(0) + 1); // Set modulus to n^{w+1}
        FMPZ a, b, c, result;
        gen_random_under_fmpz(a, modulus);
        gen_random_under_fmpz(b, modulus);
        ctx.set_modulus(modulus);
        c.mul(a, b);
        c.mul(c, b);
        for (auto _ : state) {
            benchmark::DoNotOptimize(c);
            ctx.mod(result, c);
            benchmark::DoNotOptimize(result);
        }
    }
    
    // static void PrecomputedModFmpz(benchmark::State& state) {
    //     fmpz_mod_ctx_t ctx;
    //     fmpz_t n;
    //     fmpz_init(n);
    //     fmpz_set_mpz(n, Z(test_n_squared));
    //     fmpz_mod_ctx_init(ctx, n);
    
    //     mpz_class a_mpz, b_mpz;
    //     fmpz_t a, b, c;
    //     gen_random_under(Z(a_mpz), Z(test_n_squared));
    //     gen_random_under(Z(b_mpz), Z(test_n_squared));
    //     fmpz_init(a);
    //     fmpz_init(b);
    //     fmpz_init(c);
    //     fmpz_set_mpz(a, Z(a_mpz));
    //     fmpz_set_mpz(b, Z(b_mpz));
    //     fmpz_mul(a, a, b);
    
    //     for (auto _ : state) {
    //         benchmark::DoNotOptimize(a);
    //         // benchmark::DoNotOptimize(b);
    //         fmpz_mod_set_fmpz(c, a, ctx);
    //         benchmark::DoNotOptimize(c);
    //     }
    
    //     fmpz_mod_ctx_clear(ctx);
    //     fmpz_clear(a);
    //     fmpz_clear(b);
    //     fmpz_clear(c);
    //     fmpz_clear(n);
    // }
    
    // static void ModFmpz(benchmark::State& state) {
    //     mpz_class a_mpz, b_mpz, c_mpz;
    //     fmpz_t a, b, c, n;
    //     fmpz_init(a);
    //     fmpz_init(b);
    //     fmpz_init(c);
    //     fmpz_init(n);
    
    //     fmpz_set_mpz(n, Z(test_n_squared));
    //     gen_random_under(Z(a_mpz), Z(test_n_squared));
    //     gen_random_under(Z(b_mpz), Z(test_n_squared));
    //     fmpz_set_mpz(a, Z(a_mpz));
    //     fmpz_set_mpz(b, Z(b_mpz));
    //     fmpz_mul(a, a, b);
    
    //     for (auto _ : state) {
    //         benchmark::DoNotOptimize(a);
    //         fmpz_mod(c, a, n);
    //         benchmark::DoNotOptimize(c);
    //     }
    
    //     fmpz_clear(a);
    //     fmpz_clear(b);
    //     fmpz_clear(c);
    //     fmpz_clear(n);
    // }
    
    // static void ModMulFmpzMod(benchmark::State& state) {
    //     fmpz_mod_ctx_t ctx;
    //     fmpz_t n;
    //     fmpz_init(n);
    //     fmpz_set_mpz(n, Z(test_n_squared));
    //     fmpz_mod_ctx_init(ctx, n);
    
    //     mpz_class a_mpz, b_mpz;
    //     fmpz_t a, b, c;
    //     gen_random_under(Z(a_mpz), Z(test_n_squared));
    //     gen_random_under(Z(b_mpz), Z(test_n_squared));
    //     fmpz_init(a);
    //     fmpz_init(b);
    //     fmpz_init(c);
    //     fmpz_set_mpz(a, Z(a_mpz));
    //     fmpz_set_mpz(b, Z(b_mpz));
    
    //     for (auto _ : state) {
    //         benchmark::DoNotOptimize(a);
    //         benchmark::DoNotOptimize(b);
    //         fmpz_mod_mul(c, a, b, ctx);
    //         benchmark::DoNotOptimize(c);
    //     }
    
    //     fmpz_mod_ctx_clear(ctx);
    //     fmpz_clear(a);
    //     fmpz_clear(b);
    //     fmpz_clear(c);
    //     fmpz_clear(n);
    // }
}

namespace fmpz_bench::modexp {
    static void ModExpFmpz(benchmark::State& state) {
        FMPZ base;
        FMPZ exp;
        FMPZ result;
        FMPZ modulus = test_n;
        modulus.pow(modulus, state.range(1) + 1); // Set modulus to n^{w+1}

        // Generate random base and exponent
        gen_random_under_fmpz(base, modulus);
        gen_random_fmpz(exp, state.range(0)); // Random exponent of given length
        for (auto _ : state) {
            benchmark::DoNotOptimize(base);
            benchmark::DoNotOptimize(exp);
            // Perform modular exponentiation
            result.powm(base, exp, modulus);
            benchmark::DoNotOptimize(result);
        }
    }

    static void ModExpFixedBasePrecompute(benchmark::State& state) {
        FMPZ base;
        FMPZ exp;
        FMPZ result;
        FMPZ modulus = test_n;
        modulus.pow(modulus, state.range(2) + 1); // Set modulus to n^{w+1}

        // Generate random base and exponent
        gen_random_under_fmpz(base, modulus);
        gen_random_fmpz(exp, state.range(0)); // Random exponent of given length
        for (auto _ : state) {
            // Initialize fixed base exponentiation context 
            FixedBaseExp fixed_base_exp;
            fixed_base_exp.precompute(base, modulus, state.range(0), state.range(1)); // Precompute powers of base
            benchmark::DoNotOptimize(fixed_base_exp);
        }
    }

    static void ModExpFixedBase(benchmark::State& state) {
        FMPZ base;
        FMPZ exp;
        FMPZ result;
        FMPZ modulus = test_n;
        modulus.pow(modulus, state.range(2) + 1); // Set modulus to n^{w+1}

        // Generate random base and exponent
        gen_random_under_fmpz(base, modulus);
        gen_random_fmpz(exp, state.range(0)); // Random exponent of given length
        // Initialize fixed base exponentiation context 
        FixedBaseExp fixed_base_exp;
        fixed_base_exp.precompute(base, modulus, state.range(0), state.range(1)); // Precompute powers of base
        for (auto _ : state) {
            benchmark::DoNotOptimize(base);
            benchmark::DoNotOptimize(exp);
            // Perform modular exponentiation
            // result.powm(base, exp, modulus);
            fixed_base_exp.compute(exp, result); // Use precomputed powers for exponentiation
            benchmark::DoNotOptimize(result);
        }
    }

    // static void ModExpFmpzPwrTwoExponent(benchmark::State& state) {
    //     // Exponents of the form 10000000..., will requires a lot of squaring but no multiplication
    //     FMPZ base;
    //     FMPZ exp;
    //     FMPZ result;
    //     FMPZ modulus = test_n;
    //     modulus.pow(modulus, state.range(1) + 1); // Set modulus to n^{w+1}

    //     // Generate random base and 2^k exponent
    //     gen_random_under_fmpz(base, modulus);
    //     fmpz_one_2exp(exp.get_fmpz(), state.range(0)); // Exponent is 2^k
    //     for (auto _ : state) {
    //         benchmark::DoNotOptimize(base);
    //         benchmark::DoNotOptimize(exp);
    //         // Perform modular exponentiation
    //         result.powm(base, exp, modulus);
    //         benchmark::DoNotOptimize(result);
    //     }
    // }
}

// BENCHMARK(fmpz_bench::ModMPZ)->Arg(1)->Arg(3);
// BENCHMARK(fmpz_bench::MulMPZ)->Arg(1)->Arg(3);
// BENCHMARK(fmpz_bench::SqrMPZ)->Arg(1)->Arg(3);
// BENCHMARK(fmpz_bench::ModMulMPZ)->Arg(1)->Arg(3);

BENCHMARK(fmpz_bench::MulFmpz)->Arg(1)->Arg(3);
// BENCHMARK(fmpz_bench::ModFmpz);
BENCHMARK(fmpz_bench::MulSecondFmpz)->Arg(1)->Arg(3);
BENCHMARK(fmpz_bench::PrecomputedModFmpz)->Arg(1)->Arg(3);
BENCHMARK(fmpz_bench::PrecomputedModSecondFmpz)->Arg(1)->Arg(3);
BENCHMARK(fmpz_bench::SqrFmpz)->Arg(1)->Arg(3);
BENCHMARK(fmpz_bench::SqrSecondFmpz)->Arg(1)->Arg(3);
// BENCHMARK(fmpz_bench::ModMulFmpzMod);

// BENCHMARK(fmpz_bench::modexp::ModExpFmpz)->Ranges({{128, 1024}, {1, 3}})->Args({3072,3});
// BENCHMARK(fmpz_bench::modexp::ModExpFixedBasePrecompute)
//     ->ArgsProduct({
//         benchmark::CreateRange(128, 1024, 4),
//         benchmark::CreateDenseRange(1, 4, 1),
//         {1, 3}
//     });
// BENCHMARK(fmpz_bench::modexp::ModExpFixedBase)
//     ->ArgsProduct({
//         benchmark::CreateRange(128, 1024, 4),
//         benchmark::CreateDenseRange(1, 4, 1),
//         {1, 3}
//     });
// BENCHMARK(fmpz_bench::modexp::ModExpFmpzPwrTwoExponent)->Ranges({{128, 1024}, {1, 3}})->Args({3072,3});

static void FixedBaseArgs_i(benchmark::internal::Benchmark* b, int64_t i) {
  b->Args({129, i, 1})->Args({256, i, 1})->Args({899, i, 1})
    ->Args({256, i, 3})->Args({3072, i, 3})->Args({3072 * 3, i, 3});
}

static void FixedBaseArgs(benchmark::internal::Benchmark* b) {
  FixedBaseArgs_i(b, -1);
  for (int64_t i = 1; i <= 9; ++i) {
    FixedBaseArgs_i(b, i);
  }
}

BENCHMARK(fmpz_bench::modexp::ModExpFmpz)->Args({129, 1})->Args({256, 1})->Args({899, 1})
    ->Args({256, 3})->Args({3072, 3})->Args({3072 * 3, 3});
BENCHMARK(fmpz_bench::modexp::ModExpFixedBasePrecompute)->Apply(FixedBaseArgs);
BENCHMARK(fmpz_bench::modexp::ModExpFixedBase)->Apply(FixedBaseArgs);

BENCHMARK(fmpz_bench::MulFmpzFlexSize)
          ->Args({512, 512})
          ->Args({1024, 1024})
          ->Args({2048, 2048})
          ->Args({4096, 4096})
          ->Args({8192, 8192});
// BENCHMARK(fmpz_bench::MulFmpzFlexSize)
//           ->Args({512, 2 * 512})
//           ->Args({1024, 2 * 1024})
//           ->Args({2048, 2 * 2048})
//           ->Args({4096, 2 * 4096})
//           ->Args({8192, 2 * 8192});
// BENCHMARK(fmpz_bench::MulFmpzFlexSize)
//           ->Args({64, 8192})
//           ->Args({128, 8192})
//           ->Args({256, 8192})
//           ->Args({512, 8192})
//           ->Args({1024, 8192})
//           ->Args({2048, 8192})
//           ->Args({4096, 8192})
//           ->Args({8192, 8192});
BENCHMARK(fmpz_bench::PrecomputedModFmpzFlexMod)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)->Arg(8192);
