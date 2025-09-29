#ifndef MODEXP_HH
#define MODEXP_HH

// Modular exponentiation algorithms
// - Brauer
// - Straus

#include <iostream>
#include <vector>
#include "fmpz_class.hh"
#include "fmpz_mod.hh"

namespace modexp {
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

    inline void sliding_window_nonnegative(FMPZ &result, const FMPZ &base, const FMPZ &exponent, const FMPZ &modulus) {
        assert(exponent.cmp(0) >= 0 && "Exponent must be non-negative");

        if (exponent.cmp(0) == 0) { // fmpz_sizeinbase returns 1 if exponent is 0, leading to uninitialized values
            result.set(1);
            return;
        }

        // Sliding-window exponentiation algorithm
        FMPZModCtx mod_ctx(modulus);
        FMPZ current(base);
        std::vector<FMPZ> precomputed_cache; // Cache for powers of base
        FMPZ base_sq;
        base_sq.mul(base, base);
        int exponent_num_bits = fmpz_sizeinbase(exponent.get_fmpz(), 2);
        int k = win_size(exponent_num_bits);
        precomputed_cache.reserve((1 << k) / 2); // Reserve space for precomputed powers
        precomputed_cache.push_back(current); // Store the 0-th power (base^1)
        for (int i = 3; i < (1 << k); i += 2) {
            current.mul(current, base_sq);
            mod_ctx.mod(current, current);
            precomputed_cache.push_back(current);
        }

        FMPZ a(1);
        int i = ((int) exponent_num_bits) - 1;
        while (i >= 0) {
            if (!fmpz_tstbit(exponent.get_fmpz(), i)) {
                a.mul(a, a);
                mod_ctx.mod(a, a);
                i--;
            } else {
                int l = std::max(0, i - k + 1);
                while (!fmpz_tstbit(exponent.get_fmpz(), l)) {
                    l++;
                }
                for (int j = l; j <= i; j++) {
                    a.mul(a, a);
                    mod_ctx.mod(a, a);
                }
                size_t window_value = 0;
                for (int j = l; j <= i; j++) {
                    window_value |= (fmpz_tstbit(exponent.get_fmpz(), j) << (j - l));
                }
                assert(window_value % 2 == 1 && "Window value must be odd");
                a.mul(a, precomputed_cache[window_value / 2]);
                mod_ctx.mod(a, a);
                i = l - 1; // Move to the next segment
            }
        }

        result.set(a);
    }

    inline void naive_brauer_nonnegative(FMPZ &result, const FMPZ &base, const FMPZ &exponent, const FMPZ &modulus) {
        assert(exponent.cmp(0) >= 0 && "Exponent must be non-negative");

        if (exponent.cmp(0) == 0) { // fmpz_sizeinbase returns 1 if exponent is 0, leading to uninitialized values
            result.set(1);
            return;
        }

        // Sliding-window exponentiation algorithm
        FMPZModCtx mod_ctx(modulus);
        FMPZ current(base);
        std::vector<FMPZ> precomputed_cache; // Cache for powers of base
        int exponent_num_bits = fmpz_sizeinbase(exponent.get_fmpz(), 2);
        int k = win_size(exponent_num_bits);
        precomputed_cache.reserve((1 << k) - 1); // Reserve space for precomputed powers
        precomputed_cache.push_back(current); // Store the 0-th power (base^1)
        for (int i = 2; i < (1 << k); i++) {
            current.mul(current, base);
            mod_ctx.mod(current, current);
            precomputed_cache.push_back(current);
        }
        assert(precomputed_cache.size() == (1 << k) - 1);

        FMPZ a(1);
        for (int window = exponent_num_bits - 1; window >= 0; window--) {
            if (window % k == 0) {
                size_t window_value = 0;
                for (int j = 0; j < k; j++) {
                    window_value |= (fmpz_tstbit(exponent.get_fmpz(), window + j) << j);
                }
                if (window_value != 0) {
                    a.mul(a, precomputed_cache[window_value - 1]);
                    mod_ctx.mod(a, a);
                }
            }
            if (window != 0) {
                a.mul(a, a);
                mod_ctx.mod(a, a);
            }
        }

        result.set(a);
    }

    inline void straus_nonnegative(FMPZ &result, const std::vector<FMPZ> &bases, const std::vector<FMPZ> &exponents, const FMPZ &modulus) {
        assert(bases.size() == exponents.size() && "Bases and exponents must have the same size");
        size_t n = bases.size();

        for (const auto &exponent : exponents) {
            assert(exponent.cmp(0) >= 0 && "Exponent must be non-negative");
        }

        // Sliding-window exponentiation algorithm
        FMPZModCtx mod_ctx(modulus);
        std::vector<std::vector<FMPZ>> precomputed_cache; // Caches for each exponent
        std::vector<int> ks;
        precomputed_cache.reserve(n);
        ks.reserve(n);
        int max_exponent_num_bits = 0;

        for (size_t i = 0; i < n; i++) {
            const FMPZ &base = bases[i];
            const FMPZ &exponent = exponents[i];
            FMPZ current(base);
            precomputed_cache.emplace_back();
            std::vector<FMPZ> &cache = precomputed_cache.back(); // Cache for powers of base
            int k = 0;

            if (exponent.cmp(0) != 0) {
                int exponent_num_bits = fmpz_sizeinbase(exponent.get_fmpz(), 2);

                k = win_size(exponent_num_bits);
                cache.reserve((1 << k) - 1); // Reserve space for precomputed powers
                cache.push_back(current); // Store the 0-th power (base^1)
                for (int j = 2; j < (1 << k); j++) {
                    current.mul(current, base);
                    mod_ctx.mod(current, current);
                    cache.push_back(current);
                }
                assert(cache.size() == (1 << k) - 1);

                if (exponent_num_bits > max_exponent_num_bits) {
                    max_exponent_num_bits = exponent_num_bits;
                }
            }

            ks.push_back(k);
        }

        FMPZ a(1);
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
                        a.mul(a, precomputed_cache[i][window_value - 1]);
                        mod_ctx.mod(a, a);
                    }
                }
            }
            if (window != 0) {
                a.mul(a, a);
                mod_ctx.mod(a, a);
            }
        }

        result.set(a);
    }

    inline void sliding_window(FMPZ &result, const FMPZ &base, const FMPZ &exponent, const FMPZ &modulus) {
        if (exponent.cmp(0) < 0) {
            FMPZ negated_exponent;
            fmpz_neg(negated_exponent.get_fmpz(), exponent.get_fmpz());
            FMPZ inverted_result;
            sliding_window_nonnegative(inverted_result, base, negated_exponent, modulus);
            result.invmod(inverted_result, modulus);
        } else {
            sliding_window_nonnegative(result, base, exponent, modulus);
        }
    }

    inline void naive_brauer(FMPZ &result, const FMPZ &base, const FMPZ &exponent, const FMPZ &modulus) {
        if (exponent.cmp(0) < 0) {
            FMPZ negated_exponent;
            fmpz_neg(negated_exponent.get_fmpz(), exponent.get_fmpz());
            FMPZ inverted_result;
            naive_brauer_nonnegative(inverted_result, base, negated_exponent, modulus);
            result.invmod(inverted_result, modulus);
        } else {
            naive_brauer_nonnegative(result, base, exponent, modulus);
        }
    }

    inline void straus(FMPZ &result, const std::vector<FMPZ> &bases, const std::vector<FMPZ> &exponents, const FMPZ &modulus) {
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

            straus_nonnegative(result, inverted_bases, negated_exponents, modulus);
        } else {
            straus_nonnegative(result, bases, exponents, modulus);
        }
    }
}

#endif
