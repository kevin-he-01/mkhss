#ifndef RNS_FIXED_BASE_EXP_HH
#define RNS_FIXED_BASE_EXP_HH

// Fixed-based optimization for modular exponentiation

#include <vector>
#include <optional>
#include <map>

#include "fmpz_class.hh"
#include "fmpz_mod.hh"
#include "rns_helper.hh"
#include "params.hh"

namespace rns_fbe {
    constexpr size_t FIXED_BASE_WINDOW = -1; // Special value to indicate fixed-base comb method
    inline std::map<int, std::vector<std::pair<int, size_t>>> precomputed_tune_table = { // TODO: tune this using RNS data
        // w = 1
{129, {{915, 8}, {626, 7}, {180, 6}, {76, 5}, {36, 4}, {2, -1}}},
{256, {{850, 8}, {441, 7}, {177, 6}, {87, 5}, {63, 4}, {2, -1}}},
{899, {{1003, 8}, {461, 7}, {177, 6}, {2, -1}}},
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
        RNSCtx<N2_LIMBS> *ctx;
        // FMPZ modulus; // Modulus for the exponentiation
        FMPZModCtx mod_ctx; // Context for modular arithmetic
        size_t max_exponent_bits = 0; // Maximum exponent bits for precomputation
        size_t window_size = 0; // Length of each segment in bits
        std::vector<MontInt<N2_LIMBS>> precomputed_powers_mont; // Fixed-base window only: Precomputed powers of the base
        std::vector<std::vector<MontInt<N2_LIMBS>>> power_table_mont; // Fixed-based comb only: A table of precomputed powers
        bool use_fixed_base_window = false; // Whether to use fixed-base window method
        bool initialized = false; // Flag to check if precomputation has been done
    
        public:
    
        void clear();
    
        // void downcast(const FixedBaseExp &other, const FMPZ &new_modulus);
    
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
