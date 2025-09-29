#ifndef FMPZ_MOD_CTX_HH
#define FMPZ_MOD_CTX_HH

#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include "fmpz_class.hh"

// fmpz_mod_ctx_t allow precomputing inverses for repeated mod reduction w.r.t. a fixed modulus.

// Wrapper around fmpz_mod_ctx_t to simplify memory management
// and provide a more C++-like interface.
class FMPZModCtx {
private:
    fmpz_mod_ctx_t ctx;
    FMPZ modulus;
    bool initialized = false;
    
public:
    inline FMPZModCtx() {}
    
    inline FMPZModCtx(const FMPZ &mod) {
        fmpz_mod_ctx_init(ctx, mod.get_fmpz());
        modulus.set(mod);
        initialized = true;
    }

    inline FMPZModCtx(const FMPZModCtx &other) {
        std::cerr << "Warning: not tested copy constructor for FMPZModCtx" << std::endl;
        if (other.initialized) {
            fmpz_mod_ctx_init(ctx, other.modulus.get_fmpz());
            modulus.set(other.modulus);
            initialized = true;
        } else {
            initialized = false;
        }
    }

    inline ~FMPZModCtx() {
        if (initialized)
            fmpz_mod_ctx_clear(ctx);
    }

    inline void clear() {
        if (initialized) {
            fmpz_mod_ctx_clear(ctx);
            initialized = false;
        }
    }

    inline void set_modulus(const FMPZ &mod) {
        if (!initialized) {
            fmpz_mod_ctx_init(ctx, mod.get_fmpz());
            initialized = true;
        } else {
            fmpz_mod_ctx_set_modulus(ctx, mod.get_fmpz());
        }
    }

    // Delete copy constructor and assignment operator
    // FMPZModCtx(const FMPZModCtx&) = delete;
    FMPZModCtx& operator=(const FMPZModCtx&) = delete;
    // Delete move constructor and assignment operator
    FMPZModCtx(FMPZModCtx&&) = delete;
    FMPZModCtx& operator=(FMPZModCtx&&) = delete;

    // Set dst = a mod this
    inline void mod(FMPZ &dst, const FMPZ &src) const {
        assert(initialized && "FMPZModCtx not initialized");
        fmpz_mod_set_fmpz(dst.get_fmpz(), src.get_fmpz(), ctx);
    }

    // Set dst = src^{-1} mod this
    inline void invmod(FMPZ &dst, const FMPZ &src) const {
        assert(initialized && "FMPZModCtx not initialized");
        assert(fmpz_mod_is_canonical(src.get_fmpz(), ctx) && "src must be canonical");
        fmpz_mod_inv(dst.get_fmpz(), src.get_fmpz(), ctx);
    }

    // Set dst = a^c mod this
    inline void powm(FMPZ &dst, const FMPZ &src, const FMPZ &exp) const {
        assert(initialized && "FMPZModCtx not initialized");
        assert(fmpz_mod_is_canonical(src.get_fmpz(), ctx) && "src must be canonical");
        fmpz_mod_pow_fmpz(dst.get_fmpz(), src.get_fmpz(), exp.get_fmpz(), ctx);
    }

    inline fmpz_mod_ctx_t &get_ctx() {
        assert(initialized && "FMPZModCtx not initialized");
        return ctx;
    }

    inline const fmpz_mod_ctx_t &get_ctx() const {
        assert(initialized && "FMPZModCtx not initialized");
        return ctx;
    }
};

#endif // FMPZ_MOD_CTX_HH
