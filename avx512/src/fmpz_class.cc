#include "fmpz_class.hh"

// #ifdef USE_RNS_MULTIPLICATION
#include "rns_helper.hh"
#include "rns_modexp.hh"

std::optional<RNSCtx<N2_LIMBS>> *rns_ctx;

#ifdef USE_RNS_MULTIPLICATION
void FMPZ::powm(const FMPZ &base, const FMPZ &exp, const FMPZ &mod) {
    if (mod.is_modulus) {
        assert(rns_ctx != nullptr && "RNSCtx pointer not set");
        rns_modexp::MontIntMod result;
        rns_modexp::sliding_window(result, base, exp, mod, *rns_ctx);
        (*rns_ctx)->from_mont_avx(*this, result);
    } else {
        std::cerr << "FMPZ::powm: Modulus is not set, using fmpz_powm" << std::endl;
        fmpz_powm(value, base.value, exp.value, mod.value);
    }
}
#else
void FMPZ::powm(const FMPZ &base, const FMPZ &exp, const FMPZ &mod) {
    fmpz_powm(value, base.value, exp.value, mod.value);
}
#endif
