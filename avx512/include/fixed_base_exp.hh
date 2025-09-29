#ifndef FIXED_BASE_EXP_HH
#define FIXED_BASE_EXP_HH

// Fixed-based optimization for modular exponentiation

#include <vector>
#include <optional>
#include <map>

#include "fmpz_class.hh"
#include "fmpz_mod.hh"
#include "rns_helper.hh"
#include "params.hh"

namespace fixed_base_exp {
    constexpr size_t FIXED_BASE_WINDOW = -1; // Special value to indicate fixed-base comb method
    inline std::map<int, std::vector<std::pair<int, size_t>>> precomputed_tune_table = { // TODO: tune this using FLINT data
        // w = 1
{129, {{956, 8}, {393, 7}, {194, 6}, {49, 5}, {27, 4}, {2, -1}}},
{256, {{690, 8}, {335, 7}, {144, 6}, {61, 5}, {2, -1}}},
{899, {{811, 8}, {356, 7}, {151, 6}, {2, -1}}},
    };

    inline std::optional<RNSCtx<N2_LIMBS>> hack_nullopt = std::nullopt;

    inline size_t get_optimal_window_size(size_t exponent_bits, int num_invocations) {
        // std::cerr << "DEBUG: exponent_bits: " << exponent_bits << ", num_invocations: " << num_invocations << std::endl;
        if (precomputed_tune_table.find(exponent_bits) != precomputed_tune_table.end()) {
            for (const auto &pair : precomputed_tune_table[exponent_bits]) {
                if (num_invocations >= pair.first) {
                    // std::cerr << "Using window size: " << pair.second << " for num_invocations: " << num_invocations << std::endl;
                    return pair.second;
                }
            }
            // If we reach here, it means num_invocations is smaller than all thresholds
            return 0; // No need to precompute
        } else {
            std::cerr << "No precomputed tuning data for exponent_bits: " << exponent_bits << std::endl;
            assert (false && "No precomputed tuning data for this exponent_bits");
        }
        // // Find the optimal window size based on the number of invocations
        // for (size_t i = 0; i < tune_table.size(); ++i) {
        //     if (num_invocations < tune_table[i]) {
        //         // std::cerr << "Using window size: " << i << " for num_invocations: " << num_invocations << std::endl;
        //         return i;
        //     }
        // }
        // #ifndef NO_WARNING
        // if (!warned) {
        //     warned = true;
        //     // Print a warning only once if num_invocations exceeds all thresholds
        //     // This is to avoid spamming the console with warnings
        //     std::cerr << "Warning: num_invocations " << num_invocations << " exceeds all thresholds, using largest window size." << std::endl;
        //     // TODO: extrapolate
        // }
        // #endif
        // return tune_table.size(); // Return the largest window size if num_invocations exceeds all thresholds
    }

    // inline size_t get_optimal_window_size(int num_invocations) {
    //     return num_invocations >= 2 ? 1 : 0;
    // }

    class FixedBaseExp {
    private:
        FMPZ base;
        FMPZ modulus; // Modulus for the exponentiation
        FMPZModCtx mod_ctx; // Context for modular arithmetic
        size_t max_exponent_bits = 0; // Maximum exponent bits for precomputation
        size_t window_size = 0; // Length of each segment in bits
        std::vector<FMPZ> precomputed_powers; // Fixed-base window only: Precomputed powers of the base
        std::vector<std::vector<FMPZ>> power_table; // Fixed-based comb only: A table of precomputed powers
        bool use_fixed_base_window = false; // Whether to use fixed-base window method
        bool initialized = false; // Flag to check if precomputation has been done
    
        public:
    
        void clear();
    
        void downcast(const FixedBaseExp &other, const FMPZ &new_modulus);
    
        // Precompute for all exponents up to a certain bit length
        void precompute(const FMPZ &base_value, const FMPZ &modulus, size_t exponent_bits, size_t window_length = 1, std::optional<RNSCtx<N2_LIMBS>> &ctx = hack_nullopt);
    
        void compute_positive(const FMPZ &exponent, FMPZ &result) const;
    
        void compute(const FMPZ &exponent, FMPZ &result) const;
    
        inline bool is_initialized() const {
            return initialized;
        }
    };
}

#endif
