#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/ulong_extras.h>
#include "fmpz_class.hh"

class MultiCRTCtx {
private:
    /* Data needed by multi CRT functions */
    fmpz_comb_t comb;
    fmpz_comb_temp_t comb_temp;

    ulong * primes;

    static_assert(sizeof(ulong) == 8, "This class is implemented assuming 64-bit system");

public:
    inline MultiCRTCtx(uint64_t *primes_in, int num_primes) {
        primes = new ulong[num_primes];

        memcpy(primes, primes_in, sizeof(ulong) * num_primes);

        fmpz_comb_init(comb, primes, num_primes);
        fmpz_comb_temp_init(comb_temp, comb);
    }

    // Delete copy constructor and assignment operator
    MultiCRTCtx(const MultiCRTCtx&) = delete;
    MultiCRTCtx& operator=(const MultiCRTCtx&) = delete;
    // Delete move constructor and assignment operator
    MultiCRTCtx(MultiCRTCtx&&) = delete;
    MultiCRTCtx& operator=(MultiCRTCtx&&) = delete;

    inline ~MultiCRTCtx() {
        fmpz_comb_temp_clear(comb_temp);
        fmpz_comb_clear(comb);
        delete[] primes;
    }

    // Convert integer to RNS form, i.e., reduce modulo each prime
    inline void to_rns(uint64_t *out_rns_residues, const FMPZ &input) {
        /* Reduce modulo all primes */
        fmpz_multi_mod_ui(out_rns_residues, input.get_fmpz(), comb, comb_temp);
    }

    // Convert RNS form back to integer
    inline void from_rns(FMPZ &out_integer, const uint64_t *in_rns_residues) {
         /* Reconstruct */
        fmpz_multi_CRT_ui(out_integer.get_fmpz(), in_rns_residues, comb, comb_temp, 0);
    }
};
