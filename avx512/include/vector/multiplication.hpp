#pragma once
#include "changebase.hpp"
#include "avx512ifma.hpp"
#include "../precompute_io.hpp"
#include "fmpz_class.hh"

template <template <int> class Reduction, int limbs>
inline AVXVector<limbs> ele_mult(const AVXVector<limbs> &a, const AVXVector<limbs> & b, const Reduction<limbs> &r) {
    auto ab_lo = AVXVector<limbs>(0).mullo(a, b);
    auto ab_hi = AVXVector<limbs>(0).mulhi(a, b);
    return r.reduce_small(ab_hi, ab_lo);
}

inline void split_52bit_limbs(uint64_t *limbs, size_t num_limbs, const FMPZ &n_in) {
    int sign = n_in.cmp(0);
    assert(sign >= 0 && "Negative numbers not supported");
    if (sign == 0) {
        memset(limbs, 0, num_limbs * sizeof(uint64_t));
        return;
    }

    nn_ptr n_limbs;
    FMPZ n(n_in);
    int num_64_limbs = fmpz_get_mpn(&n_limbs, n.get_fmpz());
    // std::cout << "num_64_limbs: " << num_64_limbs << std::endl;
    assert(num_limbs * 52 >= num_64_limbs * 64 && "Not enough limbs for 52-bit representation");
    size_t n_data_size = ((num_limbs + 1) / 2) * 13 + 1; // +1 to allow efficient memcpy of size 8 bytes at a time (rather than 7, which is inefficient)
    uint8_t *n_data = new uint8_t[n_data_size];
    assert(n_data_size >= num_64_limbs * 8);
    memset(n_data, 0, n_data_size);
    memcpy(n_data, n_limbs, num_64_limbs * 8);
    for (size_t i = 0; i < num_limbs; i++) {
        if (i % 2 == 0) {
            memcpy(limbs + i, n_data + (i / 2) * 13, 8);
        } else {
            memcpy(limbs + i, n_data + (i / 2) * 13 + 6, 8);
            limbs[i] >>= 4;
        }
        limbs[i] &= (1ULL << 52) - 1; // Mask to 52 bits
    }
    delete[] n_data;
    flint_free(n_limbs);
}

// inline void reconstruct_from_52bit_limbs(FMPZ &n_out, const uint64_t *limbs, size_t num_limbs) {
//     // ASSUME: upper 12 bits of limbs[i] are all zero. Undefined behavior otherwise. Please add an assertion
//     size_t n_data_size_u64 = ((num_limbs + 1) / 2) * 13 / 8 + 1; // Number of 64-bit limbs
//     size_t n_data_size = n_data_size_u64 * 8;
//     uint8_t *n_data = new uint8_t[n_data_size];
//     memset(n_data, 0, n_data_size);
//     for (size_t i = 0; i < num_limbs; i += 2) {
//         memcpy(n_data + (i / 2) * 13, limbs + i, 8);
//     }
//     for (size_t i = 1; i < num_limbs; i += 2) {
//         uint64_t *target = (uint64_t *)(n_data + (i / 2) * 13 + 6);
//         *target |= limbs[i] << 4;
//     }
//     uint64_t *n_limbs = (uint64_t *)n_data;
//     // Print n_limbs for debugging
//     // for (size_t i = 0; i < n_data_size_u64; i++) {
//     //     std::cout << "n_limbs[" << i << "]: " << std::hex << n_limbs[i] << std::dec << std::endl;
//     // }
//     slong last_nonzero = n_data_size_u64 - 1;
//     while (last_nonzero >= 0) {
//         if (n_limbs[last_nonzero] != 0) {
//             break;
//         }
//         --last_nonzero;
//     }
//     slong limb_count = last_nonzero + 1;
//     // std::cout << "limb_count: " << limb_count << std::endl;
//     if (limb_count <= 1) {
//         fmpz_set_ui(n_out.get_fmpz(), n_limbs[0]);
//     } else {
//         fmpz_set_mpn_large(n_out.get_fmpz(), n_limbs, limb_count, 0);
//     }
//     delete[] n_data;
// }

inline void fast_import_fmpz(FMPZ &n_out, uint8_t *n_data, size_t n_data_size_u64) {
    uint64_t *n_limbs = (uint64_t *)n_data;
    // Print n_limbs for debugging
    // for (size_t i = 0; i < n_data_size_u64; i++) {
    //     std::cout << "n_limbs[" << i << "]: " << std::hex << n_limbs[i] << std::dec << std::endl;
    // }
    slong last_nonzero = n_data_size_u64 - 1;
    while (last_nonzero >= 0) {
        if (n_limbs[last_nonzero] != 0) {
            break;
        }
        --last_nonzero;
    }
    slong limb_count = last_nonzero + 1;
    // std::cout << "limb_count: " << limb_count << std::endl;
    if (limb_count <= 1) {
        fmpz_set_ui(n_out.get_fmpz(), n_limbs[0]);
    } else {
        fmpz_set_mpn_large(n_out.get_fmpz(), n_limbs, limb_count, 0);
    }
}

inline void reconstruct_from_ovf_52bit_limbs(FMPZ &n_out, const uint64_t *limbs, size_t num_limbs) {
    // Upper 12 bits of limbs[i] need not be zero. This will be incorporated correctly.
    size_t n_data_size_u64 = ((num_limbs + 1) / 2) * 13 / 8 + 1; // Number of 64-bit limbs
    size_t n_data_size = n_data_size_u64 * 8;
    uint8_t *n_data_even = new uint8_t[n_data_size]; // The 1st, 3rd, 5th... limbs (even in 0-based indexing)
    uint8_t *n_data_odd = new uint8_t[n_data_size]; // The 2nd, 4th, 6th... limbs
    memset(n_data_even, 0, n_data_size);
    memset(n_data_odd, 0, n_data_size);
    for (size_t i = 0; i < num_limbs; i += 2) {
        memcpy(n_data_even + (i / 2) * 13, limbs + i, 8);
    }
    for (size_t i = 1; i < num_limbs; i += 2) {
        memcpy(n_data_odd + (i / 2) * 13, limbs + i, 8);
    }
    fast_import_fmpz(n_out, n_data_even, n_data_size_u64);
    FMPZ tmp;
    fast_import_fmpz(tmp, n_data_odd, n_data_size_u64);
    tmp.mul_2exp(tmp, 52);
    n_out.add(n_out, tmp);
    delete[] n_data_odd;
    delete[] n_data_even;
}

template<template <int> class Reduction, int limbs1, int limbs2>
class IntRNS2 {
    static constexpr int w = 52;
    static constexpr int buffer_size = AVXVector<limbs1>::VEC_LIMBS * AVXVector<limbs1>::LIMBS_PER_VEC;

    const RNSMatrix<limbs1, limbs2, true> r1;
    const RNSMatrix<limbs2, limbs1, false> r2;
    
    public:
    
    const Reduction<limbs1> m1;
    const Reduction<limbs2> m2;

    private:
    // const RNSMatrix<limbs1, limbs1, false> to_rns_m1; // makes things slightly faster (300 ns vs 1600 ns for 6144-bit integers, but cost more in terms of precomputation). Suspect some data dependency issues
    const RNSMatrix<limbs1, limbs2, false> to_rns_m2;
    const RNSMatrix<limbs2, limbs1, false> from_rns_m2;

    public:

    inline IntRNS2(PrecomputeIO &io) : r1(io), r2(io), m1(io), m2(io),
                                        to_rns_m2(io), from_rns_m2(io)
    {
    }

    inline void mul_reduce(const AVXVector<limbs1> &a1, const AVXVector<limbs2> &a2, const AVXVector<limbs1> &b1, const AVXVector<limbs2> &b2, AVXVector<limbs1> &out1, AVXVector<limbs2> &out2) {
        auto ab_m1 = ele_mult(a1, b1, m1);
        AVXVector<limbs1> unused1;
        out2 = r1.rns_reduce(ab_m1, a2, b2, m2);
        out1 = r2.rns_reduce(out2, unused1, unused1, m1);
    }

    inline void to_mont_avx(const FMPZ &input, AVXVector<limbs1> &out1, AVXVector<limbs2> &out2) {
        AVXVector<limbs1> in;
        AVXVector<limbs1> unused1;
        AVXVector<limbs2> unused2;
        uint64_t in_limbs[buffer_size];
        //// The slow O(n^2) way to get w-bit limbs from input
        // FMPZ tmp;
        // FMPZ input_copy(input);
        // for (int i = 0; i < limbs1; i++) {
        //     fmpz_tdiv_r_2exp(tmp.get_fmpz(), input_copy.get_fmpz(), w);
        //     in_limbs[i] = fmpz_get_ui(tmp.get_fmpz());
        //     fmpz_fdiv_q_2exp(input_copy.get_fmpz(), input_copy.get_fmpz(), w);
        // }
        //// The faster way
        split_52bit_limbs(in_limbs, buffer_size, input);
        in.load(in_limbs);
        //// Print in_limbs for debugging
        // in.print("in");
        // out1 = to_rns_m1.rns_reduce(in, unused1, unused1, m1);
        out2 = to_rns_m2.rns_reduce(in, unused2, unused2, m2);
        out1 = r2.rns_reduce(out2, unused1, unused1, m1);
    }

    // WARNING: output is NOT reduced mod target prime. This is the responsibility of the caller
    inline void from_mont_avx(FMPZ &output, const AVXVector<limbs2> &in2) {
        AVXVector<limbs1> out_hi;
        AVXVector<limbs1> out_lo;
        from_rns_m2.rns_reduce_raw(out_hi, out_lo, in2);
        // out_hi.print("out_hi");
        // out_lo.print("out_lo");

        static_assert(limbs1 <= buffer_size);
        uint64_t out_limbs_low[buffer_size];
        uint64_t out_limbs_high[buffer_size + 1];
        out_limbs_high[0] = 0;
        out_hi.store(out_limbs_high + 1);
        out_lo.store(out_limbs_low);
        // Print out_limbs_high and out_limbs_low
        // for (size_t i = 0; i < limbs1 + 1; i++) {
        //     std::cout << "out_limbs_high[" << i << "]: " << out_limbs_high[i] << std::endl;
        // }
        // for (size_t i = 0; i < limbs1; i++) {
        //     std::cout << "out_limbs_low[" << i << "]: " << out_limbs_low[i] << std::endl;
        // }
        reconstruct_from_ovf_52bit_limbs(output, out_limbs_low, limbs1);
        FMPZ tmp;
        reconstruct_from_ovf_52bit_limbs(tmp, out_limbs_high, limbs1 + 1);
        output.add(output, tmp);
        // std::cout << "output: " << output << std::endl;
    }
};
