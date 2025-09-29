#ifndef FUZZY_ANIKE_HH
#define FUZZY_ANIKE_HH

// Attribute-based non-interactive key exchange (ANIKE) protocol for fuzzy-word matching
#include <iostream>
#include <vector>
#include "mkhss.hh"
#include "util.hh"
#include "random_oracle.hh"
#include <cassert>

namespace fuzzy_anike {
    constexpr size_t KEY_LENGTH = 16; // Length of a single key in bytes

    typedef mkhss::party_id party_id;

    struct parameters {
        int stat_sec_param; // Tau. Statistical security parameter (for correctness)
        int sec_param; // Lambda. Computational security parameter (against offline attacks)
        int n_width; // When lambda >= 128, need n_width >= 3072
        int l;
        int total_length; // L: Number of bits in the binary string (2048 for an Iris code, for example)
        int word_length; // Width of a word in bits
        int bitcnt_l; // The number of bits needed to represent L
    };

    inline void create_parameters(parameters &params,
                                   int stat_sec_param,
                                   int sec_param,
                                   int min_n_width,
                                   int l,
                                   int b) {
        assert(stat_sec_param > 0);
        assert(sec_param > 0);
        assert(min_n_width > 0);
        assert(l > 0);
        assert(b > 0);
        params.stat_sec_param = stat_sec_param;
        params.sec_param = sec_param;
        params.n_width = min_n_width;
        params.l = l;
        params.total_length = l * b;
        params.word_length = b;
        params.bitcnt_l = bitcnt_n(l);
    }

    inline void print_parameters(std::ostream &out, const parameters &params) {
        out << "# ===== ANIKE Parameters =====" << std::endl;
        out << "# stat_sec_param: " << params.stat_sec_param << std::endl;
        out << "# sec_param: " << params.sec_param << std::endl;
        out << "# n_width: " << params.n_width << std::endl;
        out << "# l: " << params.l << std::endl;
        out << "# b: " << params.word_length << std::endl;
        out << "# ============================" << std::endl << std::endl;
    }

    struct crs {
        parameters params;
        mkhss::crs mkhss_crs; // MKHSS CRS
    };

    struct public_encoding { // pe_sigma
        mkhss::public_key mkhss_pk; // MKHSS public key
        std::vector<mkhss::public_share> x_i; // Shares of the bits of the binary string
    };

    struct private_state { // st_sigma
        mkhss::private_key mkhss_sk; // MKHSS private key
        std::vector<mkhss::private_share> x_i; // Shares of the bits of the binary string
    };

    inline void setup(crs &crs, const parameters &params, mpz_srcptr n = nullptr) {
        mkhss::parameters mkhss_params;
        crs.params = params;
        mkhss::create_parameters(mkhss_params,
                                   true, // short_exponent_assumption
                                   params.stat_sec_param,
                                   params.sec_param,
                                   params.n_width,
                                   params.bitcnt_l);// Max integer size : 2 * L
        mkhss::setup(crs.mkhss_crs, mkhss_params, n);
    }

    // x \in {0,1}^{l*b} is the binary string
    inline void attr_key_gen(const crs &crs, const std::vector<bool> &x, public_encoding &pe, private_state &st) {
        assert(x.size() == crs.params.total_length);
        FMPZ x_i;
        mkhss::keygen(crs.mkhss_crs, pe.mkhss_pk, st.mkhss_sk);
        // Generate the shares of the binary string
        pe.x_i.clear();
        pe.x_i.reserve(crs.params.total_length);
        st.x_i.clear();
        st.x_i.reserve(crs.params.total_length);
        for (int i = 0; i < crs.params.total_length; i++) {
            x_i.set(x[i] ? 1 : 0);
            pe.x_i.emplace_back();
            st.x_i.emplace_back();
            mkhss::public_share &x_i_pk = pe.x_i.back();
            mkhss::private_share &x_i_sk = st.x_i.back();
            mkhss::share(crs.mkhss_crs, pe.mkhss_pk, x_i, x_i_sk, x_i_pk);
        }
    }

    inline std::vector<uint8_t[KEY_LENGTH]> attr_key_der(const crs &crs, party_id sigma, const private_state &st_self, const public_encoding &pe_other, int threshold) {
        assert(0 <= threshold);
        assert(threshold <= crs.params.l);
        std::vector<uint8_t[KEY_LENGTH]> shared_keys(threshold);

        mkhss::rms_context ctx;
        mkhss::rms_bootstrap(crs.mkhss_crs, sigma, st_self.mkhss_sk, pe_other.mkhss_pk, ctx);
        
        mkhss::memory_share accumulator;
        mkhss::memory_share product;
        // Accumulator <- L
        mkhss::rms_mcmult(ctx.get_one(), crs.params.l, accumulator);

        for (int i = 0; i < crs.params.l; i++) {
            // Initialize product to <<1>>
            mkhss::rms_mset(product, ctx.get_one());
    
            for (int j = 0; j < crs.params.word_length; j++) {
                mkhss::input_share x_self, x_other;
                mkhss::sync_share_self(crs.mkhss_crs, ctx, st_self.mkhss_sk, pe_other.mkhss_pk, st_self.x_i[i * crs.params.word_length + j], x_self);
                mkhss::sync_share_other(crs.mkhss_crs, st_self.mkhss_sk, pe_other.mkhss_pk, pe_other.x_i[i * crs.params.word_length + j], x_other);
    
                // Programs MUST compute syntactically the same way for RMS to work, including the order (x*y vs y*x)
                mkhss::input_share &x_A = sigma == party_id::ALICE ? x_self : x_other;
                mkhss::input_share &x_B = sigma == party_id::ALICE ? x_other : x_self;
                mkhss::input_share xA_minus_xB;
                mkhss::memory_share tmp;
                mkhss::rms_isub(crs.mkhss_crs, x_A, x_B, xA_minus_xB);
                ctx.mult_serial(crs.mkhss_crs, xA_minus_xB, product, tmp, 1);
                ctx.mult_serial(crs.mkhss_crs, xA_minus_xB, tmp, tmp, 1);
                mkhss::rms_msub(product, tmp, product);
                // mkhss::rms_madd(accumulator, tmp, accumulator);
            }
    
            mkhss::rms_msub(accumulator, product, accumulator);
        }

        //// ^^^ MKHSS.Eval Above ^^^

        //// vvv ANIKE.AttrKeyDer Below vvv

        uint8_t shared_key[KEY_LENGTH]; // Instead of concatenating the identities id_A || id_B, use the NIKE shared key to simplify implementation
        PRF &shared = ctx.get_prf_shared();
        shared.evaluate_raw(0, shared_key, KEY_LENGTH);

        // auto start = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < threshold; i++) {
            FMPZ k_i;
            // ctx.output_serial(crs.mkhss_crs, key_accumulator, k_i);
            k_i.set(accumulator.y_s); // Non-black-box use of MKHSS
            size_t k_i_size = k_i.get_size_in_bytes();
            uint8_t *buffer = new uint8_t[k_i_size + sizeof(i) + KEY_LENGTH + sizeof(threshold) + 1];
            k_i.export_(buffer, k_i_size);
            serialize_integer(buffer + k_i_size, i);
            memcpy(buffer + k_i_size + sizeof(i), shared_key, KEY_LENGTH);
            serialize_integer(buffer + k_i_size + sizeof(i) + KEY_LENGTH, threshold);
            buffer[k_i_size + sizeof(i) + KEY_LENGTH + sizeof(threshold)] = random_oracle::MAGIC_ANIKE; // Magic number for ANIKE
            uint8_t hash[random_oracle::OUTPUT_LENGTH];
            random_oracle::call(buffer, k_i_size + sizeof(i) + KEY_LENGTH + sizeof(threshold) + 1, hash); // K_i = H(k_i || i || id_A || id_B || C), Description of C is just `threshold`
            static_assert(KEY_LENGTH <= random_oracle::OUTPUT_LENGTH);
            memcpy(shared_keys[i], hash, KEY_LENGTH);
            // std::cout << "Key " << i << ": " << k_i << std::endl;
            mkhss::rms_msub(accumulator, ctx.get_one(), accumulator);
            delete[] buffer;
        }

        // auto end = std::chrono::high_resolution_clock::now();
        // std::chrono::duration<double> duration = (end - start) * 1000.0; // Convert to milliseconds
        // std::cout << "Time taken for hashing: " << duration.count() << " ms" << std::endl;

        return shared_keys;
    }

    inline void attr_key_der_debug(const crs &crs, party_id sigma, const private_state &st_self, const public_encoding &pe_other, FMPZ &output) {
        mkhss::rms_context ctx;
        mkhss::rms_bootstrap(crs.mkhss_crs, sigma, st_self.mkhss_sk, pe_other.mkhss_pk, ctx);
        
        mkhss::memory_share accumulator;
        mkhss::memory_share product;
        // Accumulator <- L
        mkhss::rms_mcmult(ctx.get_one(), crs.params.l, accumulator);

        for (int i = 0; i < crs.params.l; i++) {
            // Initialize product to <<1>>
            mkhss::rms_mset(product, ctx.get_one());
    
            for (int j = 0; j < crs.params.word_length; j++) {
                mkhss::input_share x_self, x_other;
                mkhss::sync_share_self(crs.mkhss_crs, ctx, st_self.mkhss_sk, pe_other.mkhss_pk, st_self.x_i[i * crs.params.word_length + j], x_self);
                mkhss::sync_share_other(crs.mkhss_crs, st_self.mkhss_sk, pe_other.mkhss_pk, pe_other.x_i[i * crs.params.word_length + j], x_other);
    
                // Programs MUST compute syntactically the same way for RMS to work, including the order (x*y vs y*x)
                mkhss::input_share &x_A = sigma == party_id::ALICE ? x_self : x_other;
                mkhss::input_share &x_B = sigma == party_id::ALICE ? x_other : x_self;
                mkhss::input_share xA_minus_xB;
                mkhss::memory_share tmp;
                mkhss::rms_isub(crs.mkhss_crs, x_A, x_B, xA_minus_xB);
                ctx.mult_serial(crs.mkhss_crs, xA_minus_xB, product, tmp, 1);
                ctx.mult_serial(crs.mkhss_crs, xA_minus_xB, tmp, tmp, 1);
                mkhss::rms_msub(product, tmp, product);
                // mkhss::rms_madd(accumulator, tmp, accumulator);
            }
    
            mkhss::rms_msub(accumulator, product, accumulator);
        }
        
        ctx.output_serial(crs.mkhss_crs, accumulator, output);
    }
}

#endif // FUZZY_ANIKE_HH
