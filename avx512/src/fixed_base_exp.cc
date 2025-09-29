#include "fixed_base_exp.hh"

namespace fixed_base_exp {
    std::vector<int> k_table = {18, 45, 136, 385, 1056, 2863, 7400, 18801, 46540};
    
    size_t fbe_optimal_k(size_t exponent_bits) {
        for (size_t k = 0; k < k_table.size(); ++k) {
            if (exponent_bits < k_table[k]) {
                return k + 1;
            }
        }
        return k_table.size() + 1; // Return the largest window size if num_invocations exceeds all thresholds
    }

    void FixedBaseExp::clear() {
        base.set(0);
        modulus.set(0);
        mod_ctx.clear();
        max_exponent_bits = 0;
        window_size = 0;
        precomputed_powers.clear();
        power_table.clear();

        initialized = false; // Reset the initialized flag
    }

    void FixedBaseExp::downcast(const FixedBaseExp &other, const FMPZ &new_modulus) {
        if (initialized) {
            std::cerr << "Warning: FixedBaseExp already initialized. Clearing previous state." << std::endl;
            clear(); // Clear previous state if already initialized
        }

        // ASSUMES that new_modulus divides other.modulus
        base.set(other.base);
        modulus.set(new_modulus);
        mod_ctx.set_modulus(new_modulus);
        max_exponent_bits = other.max_exponent_bits;
        window_size = other.window_size;
        precomputed_powers.reserve(other.precomputed_powers.size());
        for (const auto &value : other.precomputed_powers) {
            FMPZ &new_value = precomputed_powers.emplace_back();
            mod_ctx.mod(new_value, value); // Reduce modulo the new modulus
        }

        initialized = true; // Mark as initialized after downcasting
    }
    
    void FixedBaseExp::precompute(const FMPZ &base_value, const FMPZ &modulus, size_t exponent_bits, size_t window_length, std::optional<RNSCtx<N2_LIMBS>> &ctx) {
        // std::cerr << "Precomputing with exponent_bits: " << exponent_bits
        //           << ", window_length: " << window_length << std::endl;
        assert(window_length != 0);
        if (initialized) {
            std::cerr << "Warning: FixedBaseExp already initialized. Clearing previous state." << std::endl;
            clear(); // Clear previous state if already initialized
        }
        if (window_length == FIXED_BASE_WINDOW) {
            // Use fixed-base window method
            use_fixed_base_window = true;

            window_length = fbe_optimal_k(exponent_bits);
            // std::cerr << "exponent_bits: " << exponent_bits
            //           << ", window_length: " << window_length << std::endl;
        
            this->modulus.set(modulus);
            mod_ctx.set_modulus(modulus);
            mod_ctx.mod(base, base_value); // Ensure base is reduced modulo the modulus
            max_exponent_bits = exponent_bits;
            window_size = window_length;
        
            // Precompute powers
            size_t precomputed_count = exponent_bits / window_length + (exponent_bits % window_length != 0);
            size_t t = precomputed_count - 1;
            precomputed_powers.reserve(precomputed_count);
            precomputed_powers.emplace_back(base); // Store the 0-th power (base^1)
            FMPZ current = base; // Start with the base
            // Compute g^{2^window_length}, g^{2^{2*window_length}}, ... (g is implicit)
            for (size_t i = 0; i < t; ++i) {
                for (size_t j = 0; j < window_length; ++j) {
                    current.mul(current, current);
                    mod_ctx.mod(current, current); // Reduce modulo the modulus
                }
                precomputed_powers.emplace_back(current); // Store the current power
            }
        
            // Print out precomputed powers for debugging
            // for (size_t i = 0; i < precomputed_powers.size(); ++i) {
            //     std::cout << "Element " << i << ": ";
            //     for (const auto &power : precomputed_powers[i]) {
            //         std::cout << power << " ";
            //     }
            //     std::cout << std::endl;
            // }
        } else {
            // Use fixed-base comb method
            use_fixed_base_window = false;

            this->modulus.set(modulus);
            mod_ctx.set_modulus(modulus);
            mod_ctx.mod(base, base_value); // Ensure base is reduced modulo the modulus
            max_exponent_bits = exponent_bits;
            window_size = window_length;

            // Precompute powers
            size_t num_segments = exponent_bits / window_length + (exponent_bits % window_length != 0);
            power_table.reserve(num_segments);
            FMPZ current_base = base;
            for (size_t i = 0; i < num_segments; ++i) {
                // current_base.set(base);
                // fmpz_one_2exp(tmp.get_fmpz(), i * window_length); // 2^{i * window_length}
                // current_base.powm(current_base, tmp, modulus); // base^(2^{i * window_length})
                auto &segment = power_table.emplace_back();
                segment.reserve((1 << window_length) - 1); // 2^window_length entries

                segment.emplace_back(current_base); // Start with base^(2^{i * window_length})
                for (size_t j = 2; j < (1 << window_length); ++j) {
                    FMPZ &current_power = segment.back();
                    FMPZ &next_power = segment.emplace_back();
                    if (j % 2 == 0) {
                        // Even index, squaring is faster than multiplying by the base
                        FMPZ &j_over_2_th = segment[j / 2 - 1]; // j/2-th power 
                        next_power.mul(j_over_2_th, j_over_2_th);
                    } else {
                        // Odd index, multiply by the base
                        next_power.mul(current_power, current_base);
                    }
                    mod_ctx.mod(next_power, next_power);
                }
                current_base.mul(current_base, segment.back()); // Update current_base to the last computed power
                mod_ctx.mod(current_base, current_base); // Ensure it is reduced modulo the modulus
            }

            // Print out power table for debugging
            // for (size_t i = 0; i < power_table.size(); ++i) {
            //     std::cout << "Segment " << i << ": ";
            //     for (const auto &power : power_table[i]) {
            //         std::cout << power << " ";
            //     }
            //     std::cout << std::endl;
            // }
        }

        initialized = true; // Mark as initialized after precomputation
    }
    
    void FixedBaseExp::compute_positive(const FMPZ &exponent, FMPZ &result) const {
        assert(initialized && "FixedBaseExp not initialized. Call precompute() first.");

        FMPZ tmp;

        if (use_fixed_base_window) {
            // Use fixed-base window method
            // Process the exponent in segments
            assert(exponent.cmp(0) >= 0 && "Exponent must be non-negative");
            size_t exponent_bits = fmpz_sizeinbase(exponent.get_fmpz(), 2);
            if (exponent_bits > max_exponent_bits) {
                std::cerr << "Exponent bits " << exponent_bits << " exceed precomputed range of " << max_exponent_bits << std::endl;
                assert(false && "Exponent exceeds precomputed range");
            }

            std::vector<ulong> segment_values;
            segment_values.reserve(precomputed_powers.size());

            for (size_t i = 0; i < precomputed_powers.size(); ++i) {
                size_t segment_start = i * window_size;
                size_t segment_end = segment_start + window_size;
                // TODO: probably a more efficient way to do this (currently quadratic in width of exponent)
                fmpz_fdiv_q_2exp(tmp.get_fmpz(), exponent.get_fmpz(), segment_start);
                FMPZ segment_value_fmpz;
                fmpz_fdiv_r_2exp(segment_value_fmpz.get_fmpz(), tmp.get_fmpz(), window_size);

                ulong segment_value = fmpz_get_ui(segment_value_fmpz.get_fmpz());
                segment_values.push_back(segment_value);
            }

            FMPZ a(1), b(1);
            for (ulong j = (1 << window_size) - 1; j >= 1; --j) {
                for (size_t i = 0; i < segment_values.size(); ++i) {
                    if (segment_values[i] == j) {
                        b.mul(b, precomputed_powers[i]); // Multiply by the precomputed power
                        mod_ctx.mod(b, b); // Reduce modulo the modulus
                        // TODO: micro optimization: don't multiply if b is already 1
                    }
                }
                a.mul(a, b); // Multiply the accumulated result
                mod_ctx.mod(a, a); // Reduce modulo the modulus
            }

            result.set(a); // Set the final result
        } else {
            // Use fixed-base comb method
            FMPZ a(1);

            // Process the exponent in segments
            assert(exponent.cmp(0) >= 0 && "Exponent must be non-negative");
            size_t exponent_bits = fmpz_sizeinbase(exponent.get_fmpz(), 2);
            if (exponent_bits > max_exponent_bits) {
                std::cerr << "Exponent bits " << exponent_bits << " exceed precomputed range of " << max_exponent_bits << std::endl;
                assert(false && "Exponent exceeds precomputed range");
            }
            for (size_t i = 0; i < power_table.size(); ++i) {
                size_t segment_start = i * window_size;
                size_t segment_end = segment_start + window_size;
                // TODO: probably a more efficient way to do this (currently quadratic in width of exponent)
                fmpz_fdiv_q_2exp(tmp.get_fmpz(), exponent.get_fmpz(), segment_start);
                FMPZ segment_value_fmpz;
                fmpz_fdiv_r_2exp(segment_value_fmpz.get_fmpz(), tmp.get_fmpz(), window_size);

                if (segment_value_fmpz.cmp(0) > 0) {
                    ulong segment_value = fmpz_get_ui(segment_value_fmpz.get_fmpz());
                    // std::cout << "Segment " << i << ": segment_value = " << segment_value << std::endl;
                    // assert(segment_value < (1UL << window_size) && "Segment value exceeds precomputed range");
                    const FMPZ &power = power_table[i][segment_value - 1]; // -1 because we skip the 0-th power
                    // Multiply the result by the precomputed power
                    a.mul(a, power);
                    mod_ctx.mod(a, a); // Modulo reduction
                }
            }

            result.set(a); // Set the final result
        }
    }

    void FixedBaseExp::compute(const FMPZ &exponent, FMPZ &result) const {
        if (exponent.cmp(0) < 0) {
            FMPZ negated_exponent;
            fmpz_neg(negated_exponent.get_fmpz(), exponent.get_fmpz());
            FMPZ inverted_result;
            compute_positive(negated_exponent, inverted_result);
            // inverted_result.invmod(inverted_result, modulus);
            mod_ctx.invmod(inverted_result, inverted_result);
            result.set(inverted_result);
        } else {
            compute_positive(exponent, result);
        }
    }

}
