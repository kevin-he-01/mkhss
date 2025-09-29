#pragma once
#include "../vector/avx512ifma.hpp"

template <int limbs>
class MontgomeryReduce {
    public:

    AVXVector<limbs> mont_factor;
    AVXVector<limbs> moduli;
    AVXVector<limbs> t_squared;

    static constexpr int mulbits = 52;
    static constexpr int bits = 52;
    static constexpr int hi_bit_shift = bits + bits - mulbits;

    const AVXScalar hi_mask = get_mask(hi_bit_shift);
    const AVXScalar lo_mask = get_mask(mulbits);

    inline MontgomeryReduce(PrecomputeIO& io) {
        io.read(this->moduli);
        io.read(this->mont_factor);
        io.read(this->t_squared);
    }

    inline AVXVector<limbs> reduce_small(const AVXVector<limbs> &hi, const AVXVector<limbs> &low) const {
        auto q = AVXVector<limbs>(0).mullo(low, mont_factor);
        auto h = AVXVector<limbs>(0).mulhi(q, moduli);
        auto c = hi.sub(h);
        return c.cmp_add(h, hi, moduli);
    }

    // inline AVXVector<limbs> reduce_full(const AVXVector<limbs> &hi, const AVXVector<limbs> &low) const {
        
    //     auto hi_hi = hi.srli(hi_bit_shift);
    //     auto hi_rem = hi.and_scalar(hi_mask);
    //     auto lo_acc = low.mullo(hi_hi, t_squared);

    //     auto lo_hi = lo_acc.srli(mulbits);
    //     auto lo_rem = lo_acc.and_scalar(lo_mask);
    //     auto hi_acc = hi_rem.add(lo_hi);
    //     auto hi_cor = hi_acc.cmp_sub(moduli);
    //     return reduce_small(hi_cor, lo_rem);
        
    // }

     inline AVXVector<limbs> reduce_full(const AVXVector<limbs> &hi, const AVXVector<limbs> &low) const {
        auto lo_hi = low.srli(mulbits);
        auto lo_rem = low.and_scalar(lo_mask);
        auto hi_t = hi.add(lo_hi);

        auto hi_hi = hi_t.srli(hi_bit_shift);
        auto hi_rem = hi_t.and_scalar(hi_mask);
        auto hi_acc = hi_rem.mullo(hi_hi, t_squared);

        auto hi_cor = hi_acc.cmp_sub(moduli);
        return reduce_small(hi_cor, lo_rem);
        
    }

};