#pragma once
#include <vector>
#include <optional>

#include "precompute_io.hpp"
#include "vector/multiplication.hpp"
#include "reduction/montgomery.hpp"

#include "fmpz_class.hh"
#include "fmpz_mod.hh"
#include "multi_crt.hh"
#include "util.hh"

template <int LIMBS>
struct MontInt {
    AVXVector<LIMBS> a1;
    AVXVector<LIMBS> a2;
    bool dummy_one = false; // representation of 1 (dummy)

    inline MontInt() = default;

    inline MontInt(bool dummy_one) : dummy_one(dummy_one) {
        // This constructor is used to represent 1 in Montgomery form
    }

    inline MontInt(const AVXVector<LIMBS> &a1, const AVXVector<LIMBS> &a2) : a1(a1), a2(a2) {}

    void print() const {
        if (dummy_one) {
            std::cout << "MontInt(1)" << std::endl;
            return;
        }
        std::cout << "MontInt(";
        a1.print("  .a1");
        a2.print("  .a2");
        std::cout << ")" << std::endl;
    }
};

inline double us_to_rns_total = 0.0;
inline double us_from_rns_total = 0.0;

template <int LIMBS>
class RNSCtx {
    using RNSMultiplier = IntRNS2<MontgomeryReduce, LIMBS, LIMBS>;
    public:
    RNSMultiplier multiplier;
    // std::vector<uint64_t> moduli1;
    FMPZ modulus1;
    std::vector<uint64_t> moduli2;
    FMPZ modulus2;
    FMPZModCtx modulus2_ctx;
    FMPZ inv_m1_mod_m2;
    FMPZ target;
    FMPZModCtx target_ctx;
    std::optional<MultiCRTCtx> m1_crt_ctx;
    // std::optional<MultiCRTCtx> m2_crt_ctx;

    // from_mont precomputed stuff
    FMPZ modulus1_half;
    FMPZModCtx modulus1_ctx;
    // std::vector<FMPZ> icrt_factors;
    FMPZ inv_2w_mod_m1;
    FMPZ inv_m1_mod_target;

    static constexpr int w = 52;
    static constexpr int buffer_size = AVXVector<LIMBS>::VEC_LIMBS * AVXVector<LIMBS>::LIMBS_PER_VEC;

    // Delete copy constructor
    RNSCtx(const RNSCtx&) = delete;
    // Delete move constructor
    RNSCtx(RNSCtx&&) = delete;
    // Delete copy assignment operator
    RNSCtx& operator=(const RNSCtx&) = delete;
    // Delete move assignment operator
    RNSCtx& operator=(RNSCtx&&) = delete;

    inline RNSCtx(PrecomputeIO &io, const FMPZ &target_modulus) : multiplier(io) {
        auto reducer_m1 = multiplier.m1;
        auto reducer_m2 = multiplier.m2;

        std::vector<uint64_t> moduli1;
        moduli1.resize(buffer_size);
        reducer_m1.moduli.store(moduli1.data());
        moduli1.resize(LIMBS);
        moduli2.resize(buffer_size);
        reducer_m2.moduli.store(moduli2.data());
        moduli2.resize(LIMBS);
        // modulus1 = product of moduli1
        modulus1.set(1);
        for (const auto &mod : moduli1) {
            modulus1.mul(modulus1, mod);
        }
        fmpz_tdiv_q_ui(modulus1_half.get_fmpz(), modulus1.get_fmpz(), 2);
        modulus1_ctx.set_modulus(modulus1);
        m1_crt_ctx.emplace(moduli1.data(), moduli1.size());

        // modulus2 = product of moduli2
        modulus2.set(1);
        for (const auto &mod : moduli2) {
            modulus2.mul(modulus2, mod);
        }
        modulus2_ctx.set_modulus(modulus2);
        // m2_crt_ctx.emplace(moduli2.data(), moduli2.size());

        target.set(target_modulus);
        target_ctx.set_modulus(target);

        inv_m1_mod_m2.invmod(modulus1, modulus2);

        // std::cout << "moduli1: ";
        // for (const auto &mod : moduli1) {
        //     std::cout << mod << " ";
        // }
        // std::cout << std::endl;
        // std::cout << "moduli2: ";
        // for (const auto &mod : moduli2) {
        //     std::cout << mod << " ";
        // }
        // std::cout << std::endl;
        // std::cout << "modulus1: " << modulus1 << std::endl;
        // std::cout << "modulus2: " << modulus2 << std::endl;
        // std::cout << "inv_m1_mod_m2: " << inv_m1_mod_m2 << std::endl;
        // std::cout << "target: " << target << std::endl;

        // icrt_factors.reserve(LIMBS);
        // for (const auto &mi : moduli1) {
        //     FMPZ m1_over_mi;
        //     fmpz_fdiv_q_ui(m1_over_mi.get_fmpz(), modulus1.get_fmpz(), mi);
        //     // std::cout << "m1_over_mi: " << m1_over_mi << std::endl;
        //     auto &icrt_factor = icrt_factors.emplace_back();
        //     icrt_factor.invmod(m1_over_mi, mi);
        //     // std::cout << "icrt_factor: " << icrt_factor << std::endl;
        //     icrt_factor.mul(icrt_factor, m1_over_mi);
        //     modulus1_ctx.mod(icrt_factor, icrt_factor);
        // }

        fmpz_ui_pow_ui(inv_2w_mod_m1.get_fmpz(), 2, w);
        inv_2w_mod_m1.invmod(inv_2w_mod_m1, modulus1);
        inv_m1_mod_target.invmod(modulus1, target);

        // std::cout << "icrt_factors: ";
        // for (const auto &factor : icrt_factors) {
        //     std::cout << factor << " ";
        // }
        // std::cout << std::endl;
        // std::cout << "inv_2w_mod_m1: " << inv_2w_mod_m1 << std::endl;
        // std::cout << "inv_m1_mod_target: " << inv_m1_mod_target << std::endl;
    }

    // SLOW! Use to_mont_avx instead
    inline void to_mont(AVXVector<LIMBS> &a1, AVXVector<LIMBS> &a2, const FMPZ &a) {
        std::cerr << "RNSCtx::to_mont: This function is slow and should only be used for debugging!";
        
        FMPZ tmp;
        tmp.mul(a, modulus1);
        target_ctx.mod(tmp, tmp);
        tmp.mul_2exp(tmp, w);

        uint64_t v1[buffer_size];
        // us_to_rns_total += measure_time_seconds([&]() {
            // for (int i = 0; i < LIMBS; i++) {
            //     FMPZ v1_i;
            //     v1_i.mod(tmp, moduli1[i]);
            //     v1[i] = fmpz_get_ui(v1_i.get_fmpz());
            // }
            m1_crt_ctx->to_rns(v1, tmp);
        // }) * 1e6; // Convert to us

        a1.load(v1);

        tmp.mul(inv_m1_mod_m2, tmp);
        modulus2_ctx.mod(tmp, tmp);

        uint64_t v2[buffer_size];
        // us_to_rns_total += measure_time_seconds([&]() {
            for (int i = 0; i < LIMBS; i++) {
                FMPZ v2_i;
                v2_i.mod(tmp, moduli2[i]);
                v2[i] = fmpz_get_ui(v2_i.get_fmpz());
            }
            // m2_crt_ctx->to_rns(v2, tmp);
        // }) * 1e6; // Convert to us
        
        a2.load(v2);
    }

    // SLOW! Use to_mont_avx instead
    inline void to_mont(MontInt<LIMBS> &mont, const FMPZ &a) {
        to_mont(mont.a1, mont.a2, a);
        mont.dummy_one = false;
    }

    inline void to_mont_avx(MontInt<LIMBS> &mont, const FMPZ &a) {
        multiplier.to_mont_avx(a, mont.a1, mont.a2);
        mont.dummy_one = false;
    }

    // SLOW! Use from_mont_avx instead
    inline void from_mont(FMPZ &a, const AVXVector<LIMBS> &a1/*, const AVXVector<LIMBS> &a2*/) {
        std::cerr << "RNSCtx::from_mont: This function is slow and should only be used for debugging!";

        FMPZ acc;
        FMPZ tmp;

        uint64_t v1[buffer_size];
        a1.store(v1);

        // us_from_rns_total += measure_time_seconds([&]() {
            m1_crt_ctx->from_rns(acc, v1);
            // for (int i = 0; i < LIMBS; i++) {
            //     fmpz_mul_ui(tmp.get_fmpz(), icrt_factors[i].get_fmpz(), v1[i]);
            //     acc.add(acc, tmp);
            // }
            // modulus1_ctx.mod(acc, acc);
        // }) * 1e6; // Convert to us
        a.mul(acc, inv_2w_mod_m1);
        modulus1_ctx.mod(a, a);
        if (a.cmp(modulus1_half) > 0) { // A very rare corner case where a is "negative" mod modulus1 due to redundancy
            a.sub(a, modulus1);
        }
        a.mul(a, inv_m1_mod_target);
        target_ctx.mod(a, a);
    }

    inline void from_mont(FMPZ &a, const MontInt<LIMBS> &mont) {
        if (mont.dummy_one) {
            a.set(1); // If it's the dummy representation of 1, just set a to 1
            return;
        }
        from_mont(a, mont.a1/*, mont.a2*/);
    }

    inline void from_mont_avx(FMPZ &a, const MontInt<LIMBS> &mont) {
        if (mont.dummy_one) {
            a.set(1); // If it's the dummy representation of 1, just set a to 1
            return;
        }
        multiplier.from_mont_avx(a, mont.a2);
        target_ctx.mod(a, a);
    }

    inline void mul_mont(MontInt<LIMBS> &c, const MontInt<LIMBS> &a, const MontInt<LIMBS> &b) {
        if (a.dummy_one) {
            c = b; // If a is the dummy representation of 1, just set c to b
            return;
        }
        if (b.dummy_one) {
            c = a; // If b is the dummy representation of 1, just set c to a
            return;
        }
        multiplier.mul_reduce(a.a1, a.a2, b.a1, b.a2, c.a1, c.a2);
        c.dummy_one = false;
    }
};
