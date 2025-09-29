#ifndef RNS_MODEXP_HH
#define RNS_MODEXP_HH

// Modular exponentiation algorithms
// - Brauer
// - Straus

#include <iostream>
#include <vector>
#include <optional>

#include "fmpz_class.hh"
#include "fmpz_mod.hh"
#include "rns_helper.hh"
#include "params.hh"

namespace rns_modexp {
    // TODO: port over Straus algorithm to do 2 exponentiations at nearly the cost of one
    typedef MontInt<N2_LIMBS> MontIntMod;
    inline std::optional<RNSCtx<N2_LIMBS>> hack_nullopt = std::nullopt;

    // Taken from GMP
    static inline int
    win_size (mp_bitcnt_t eb)
    {
        int k;
        static mp_bitcnt_t x[] = {7,25,81,241,673,1793,4609,11521,28161,~(mp_bitcnt_t)0};
        for (k = 0; eb > x[k++]; )
            ;
        return k;
    }

    inline void sliding_window_nonnegative(MontIntMod &result, const MontIntMod &base, const FMPZ &exponent, std::optional<RNSCtx<N2_LIMBS>> &ctx = hack_nullopt) {
        assert(exponent.cmp(0) >= 0 && "Exponent must be non-negative");
        assert(ctx.has_value() && "RNSCtx must be provided for sliding window exponentiation");

        if (exponent.cmp(0) == 0) { // fmpz_sizeinbase returns 1 if exponent is 0, leading to uninitialized values
            // result.set(1);
            result.dummy_one = true; // Use the dummy representation of 1
            return;
        }

        // Sliding-window exponentiation algorithm
        MontIntMod current_mont = base;
        std::vector<MontIntMod> precomputed_cache_mont; // Cache for powers of base
        MontIntMod base_sq_mont;
        ctx->mul_mont(base_sq_mont, current_mont, current_mont); // Square the base in Montgomery
        int exponent_num_bits = fmpz_sizeinbase(exponent.get_fmpz(), 2);
        int k = win_size(exponent_num_bits);
        precomputed_cache_mont.reserve((1 << k) / 2); // Reserve space for precomputed powers
        precomputed_cache_mont.push_back(current_mont); // Store the 0-th power (base^1)
        for (int i = 3; i < (1 << k); i += 2) {
            // current.mul(current, base_sq);
            // mod_ctx.mod(current, current);
            // precomputed_cache.push_back(current);
            ctx->mul_mont(current_mont, current_mont, base_sq_mont); // Square the current power in Montgomery
            precomputed_cache_mont.push_back(current_mont);
        }

        MontIntMod a_mont(true); // Start with the dummy representation of 1
        int i = ((int) exponent_num_bits) - 1;
        while (i >= 0) {
            if (!fmpz_tstbit(exponent.get_fmpz(), i)) {
                // a.mul(a, a);
                // mod_ctx.mod(a, a);
                ctx->mul_mont(a_mont, a_mont, a_mont); // Square the accumulated result in Montgomery
                i--;
            } else {
                int l = std::max(0, i - k + 1);
                while (!fmpz_tstbit(exponent.get_fmpz(), l)) {
                    l++;
                }
                for (int j = l; j <= i; j++) {
                    // a.mul(a, a);
                    // mod_ctx.mod(a, a);
                    ctx->mul_mont(a_mont, a_mont, a_mont); // Square the accumulated result in Montgomery
                }
                size_t window_value = 0;
                for (int j = l; j <= i; j++) {
                    window_value |= (fmpz_tstbit(exponent.get_fmpz(), j) << (j - l));
                }
                assert(window_value % 2 == 1 && "Window value must be odd");
                ctx->mul_mont(a_mont, a_mont, precomputed_cache_mont[window_value / 2]); // Multiply by the precomputed power in Montgomery
                i = l - 1; // Move to the next segment
            }
        }

        result = a_mont;
    }

    inline void sliding_window(MontIntMod &result, const FMPZ &base, const FMPZ &exponent, const FMPZ &modulus, std::optional<RNSCtx<N2_LIMBS>> &ctx = hack_nullopt) {
        MontIntMod base_mont;
        if (exponent.cmp(0) < 0) {
            FMPZ negated_exponent;
            fmpz_neg(negated_exponent.get_fmpz(), exponent.get_fmpz());
            FMPZ inverted_base;
            inverted_base.invmod(base, modulus);
            ctx->to_mont_avx(base_mont, inverted_base);
            sliding_window_nonnegative(result, base_mont, negated_exponent, ctx);
        } else {
            ctx->to_mont_avx(base_mont, base);
            sliding_window_nonnegative(result, base_mont, exponent, ctx);
        }
    }

    inline void straus_nonnegative(MontIntMod &result, const std::vector<MontIntMod> &bases, const std::vector<FMPZ> &exponents, std::optional<RNSCtx<N2_LIMBS>> &ctx = hack_nullopt) {
        assert(bases.size() == exponents.size() && "Bases and exponents must have the same size");
        size_t n = bases.size();

        for (const auto &exponent : exponents) {
            assert(exponent.cmp(0) >= 0 && "Exponent must be non-negative");
        }

        // Sliding-window exponentiation algorithm
        std::vector<std::vector<MontIntMod>> precomputed_cache; // Caches for each exponent
        std::vector<int> ks;
        precomputed_cache.reserve(n);
        ks.reserve(n);
        int max_exponent_num_bits = 0;

        for (size_t i = 0; i < n; i++) {
            const MontIntMod &base = bases[i];
            const FMPZ &exponent = exponents[i];
            MontIntMod current(base);
            precomputed_cache.emplace_back();
            std::vector<MontIntMod> &cache = precomputed_cache.back(); // Cache for powers of base
            int k = 0;

            if (exponent.cmp(0) != 0) {
                int exponent_num_bits = fmpz_sizeinbase(exponent.get_fmpz(), 2);

                k = win_size(exponent_num_bits);
                cache.reserve((1 << k) - 1); // Reserve space for precomputed powers
                cache.push_back(current); // Store the 0-th power (base^1)
                for (int j = 2; j < (1 << k); j++) {
                    ctx->mul_mont(current, current, base); // Multiply the current power in Montgomery
                    cache.push_back(current);
                }
                assert(cache.size() == (1 << k) - 1);

                if (exponent_num_bits > max_exponent_num_bits) {
                    max_exponent_num_bits = exponent_num_bits;
                }
            }

            ks.push_back(k);
        }

        MontIntMod a(1);
        for (int window = max_exponent_num_bits - 1; window >= 0; window--) {
            for (size_t i = 0; i < n; i++) {
                const FMPZ &exponent = exponents[i];
                int k = ks[i];
                if (k == 0) continue; // Skip if exponent is 0
                if (window % k == 0) {
                    size_t window_value = 0;
                    for (int j = 0; j < k; j++) {
                        window_value |= (fmpz_tstbit(exponent.get_fmpz(), window + j) << j);
                    }
                    if (window_value != 0) {
                        ctx->mul_mont(a, a, precomputed_cache[i][window_value - 1]); // Multiply by the precomputed power in Montgomery
                    }
                }
            }
            if (window != 0) {
                ctx->mul_mont(a, a, a); // Square the accumulated result in Montgomery
            }
        }

        result = a;
    }

    inline void straus_nonnegative_fmpz(FMPZ &result, const std::vector<FMPZ> &bases, const std::vector<FMPZ> &exponents, std::optional<RNSCtx<N2_LIMBS>> &ctx = hack_nullopt) {
        std::vector<MontIntMod> bases_mont;
        bases_mont.reserve(bases.size());
        for (const auto &base : bases) {
            MontIntMod base_mont;
            ctx->to_mont_avx(base_mont, base);
            bases_mont.push_back(base_mont);
        }
        MontIntMod result_mont;
        straus_nonnegative(result_mont, bases_mont, exponents, ctx);
        ctx->from_mont_avx(result, result_mont);
    }

    inline void straus(FMPZ &result, const std::vector<FMPZ> &bases, const std::vector<FMPZ> &exponents, const FMPZ &modulus, std::optional<RNSCtx<N2_LIMBS>> &ctx = hack_nullopt) {
        assert(bases.size() == exponents.size() && "Bases and exponents must have the same size");
        size_t n = bases.size();

        bool any_negative = false;
        for (const auto &exponent : exponents) {
            if (exponent.cmp(0) < 0) {
                any_negative = true;
                break;
            }
        }

        if (any_negative) {
            std::vector<FMPZ> negated_exponents;
            std::vector<FMPZ> inverted_bases;
            negated_exponents.reserve(exponents.size());
            for (size_t i = 0; i < n; i++) {
                const FMPZ &exponent = exponents[i];
                const FMPZ &base = bases[i];
                if (exponent.cmp(0) < 0) {
                    FMPZ negated_exponent;
                    FMPZ inverted_base;
                    inverted_base.invmod(base, modulus);
                    fmpz_neg(negated_exponent.get_fmpz(), exponent.get_fmpz());
                    inverted_bases.push_back(inverted_base);
                    negated_exponents.push_back(negated_exponent);
                } else {
                    inverted_bases.push_back(base);
                    negated_exponents.push_back(exponent);
                }
            }

            straus_nonnegative_fmpz(result, inverted_bases, negated_exponents, ctx);
        } else {
            straus_nonnegative_fmpz(result, bases, exponents, ctx);
        }
    }
}

#endif
