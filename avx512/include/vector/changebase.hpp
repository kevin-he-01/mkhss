#pragma once
#include "avx512ifma.hpp"
#include "../precompute_io.hpp"

template <int in_limbs1, int out_limbs, bool accumulate>
class RNSMatrix {

    private:
        AVXVector<out_limbs> rns_mat[in_limbs1];
        uint64_t shifted_quotient_estimations[in_limbs1];
        AVXVector<out_limbs> correction;

        const int precision = 64;
        static constexpr int mulbits = 52;

    public:

    inline RNSMatrix(PrecomputeIO &io) {
        // printf("Before %d\n", io.read_counter);
        for (int i = 0; i < in_limbs1; i++) {
            io.read(rns_mat[i]);
        }
        // printf("RNS mat %d %d %d\n", io.read_counter, in_limbs1 + in_limbs2, out_limbs);
        io.vector(shifted_quotient_estimations);
        io.read(correction);
        // io.vector(target_shifted_quotient_estimations);
        // io.scalar(tqc);
        // io.read(tp);
    }

    template<int in_limbs, int offset>
    inline void accumulate_loop(AVXVector<out_limbs> &out_hi, AVXVector<out_limbs> &out_lo, const AVXVector<in_limbs> &residues, uint128_t &k_raw) const {
        uint64_t scalars[nonzero(residues.VEC_LIMBS * residues.LIMBS_PER_VEC)];
        residues.store((uint64_t*)&scalars);
        for (int i = 0; i < in_limbs; i++) {
            uint64_t s = scalars[i];
            AVXScalar sv = Scalar(s);
            out_hi = out_hi.mulhi_scalar(rns_mat[i + offset], sv);
            out_lo = out_lo.mullo_scalar(rns_mat[i + offset], sv);
            k_raw += (uint128_t)(s) * (uint128_t)(shifted_quotient_estimations[i + offset]);
        }
    }

    inline uint64_t add_correction(AVXVector<out_limbs> &out_hi, AVXVector<out_limbs> &out_lo, uint128_t raw, const AVXVector<out_limbs> &c_value) const {
        uint64_t k = raw >> precision;
        AVXScalar sk = Scalar(k);
        out_hi = out_hi.mulhi_scalar(c_value, sk);
        out_lo = out_lo.mullo_scalar(c_value, sk);
        sk = scalar_shift(sk, mulbits);
        out_hi = out_hi.mullo_scalar(c_value, sk);
        AVXVector<out_limbs> lp = AVXVector<out_limbs>(0).mulhi_scalar(c_value, sk);
        out_hi = out_hi.add(lp.slli(mulbits));
        return k;
    }

    template<class Reducer>
    inline AVXVector<out_limbs> rns_reduce(const AVXVector<in_limbs1> &residues1, const AVXVector<out_limbs> &acc1, const AVXVector<out_limbs> &acc2, const Reducer& reducer) const {
        AVXVector<out_limbs> out_hi = AVXVector<out_limbs>(0);
        AVXVector<out_limbs> out_lo = AVXVector<out_limbs>(0);

        if (accumulate) {
            out_hi = out_hi.mulhi(acc1, acc2);
            out_lo = out_lo.mullo(acc1, acc2);
        }
        uint128_t k_raw = 0;
        accumulate_loop<in_limbs1, 0>(out_hi, out_lo, residues1, k_raw);
        uint64_t k = add_correction(out_hi, out_lo, k_raw, correction);
        return reducer.reduce_full(out_hi, out_lo);
    }

    inline void rns_reduce_raw(AVXVector<out_limbs> &out_hi, AVXVector<out_limbs> &out_lo, const AVXVector<in_limbs1> &residues1) const {
        // useful for convert from RMRNS to integer
        out_hi = AVXVector<out_limbs>(0);
        out_lo = AVXVector<out_limbs>(0);

        assert(!accumulate && "accumulate must be false for rns_reduce_raw");
        uint128_t k_raw = 0;
        accumulate_loop<in_limbs1, 0>(out_hi, out_lo, residues1, k_raw);
        uint64_t k = add_correction(out_hi, out_lo, k_raw, correction);
    }
};