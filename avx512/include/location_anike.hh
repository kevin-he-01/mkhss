#ifndef GEOLOCATION_ANIKE_HH
#define GEOLOCATION_ANIKE_HH

// Attribute-based non-interactive key exchange (ANIKE) protocol for geolocaiton proximity
#include <iostream>
#include <vector>
#include <cassert>
#include <iomanip>
#include <optional>
#include <nlohmann/json.hpp>

#include "util.hh"
#include "random_oracle.hh"
#include "anike.hh"
#include "nike.hh"

#ifdef BASELINE_MKHSS
#include "baseline_mkhss.hh"
#define mkhss baseline_mkhss
#else
#include "mkhss.hh"
#endif

namespace anike::geolocation {
    using namespace anike;

    typedef mkhss::party_id party_id;
    typedef uint64_t coord_t; // Type for a single coordinate, e.g., latitude or longitude

    struct parameters {
        int stat_sec_param; // Tau. Statistical security parameter (for correctness)
        int sec_param; // Lambda. Computational security parameter (against offline attacks)
        int n_width; // When lambda >= 128, need n_width >= 3072
        int D; // Dimensionality of the coordinate (e.g., 2 for latitude and longitude)
        int l;
    };

    inline void create_parameters(parameters &params,
                                   int stat_sec_param,
                                   int sec_param,
                                   int min_n_width,
                                   int l,
                                   int D = 2) {
        assert(stat_sec_param > 0);
        assert(sec_param > 0);
        assert(min_n_width > 0);
        assert(l > 0);
        assert(l <= sizeof(coord_t) * 8); // Ensure l is reasonable for coord_t
        params.stat_sec_param = stat_sec_param;
        params.sec_param = sec_param;
        params.n_width = min_n_width;
        params.l = l;
        params.D = D; // Dimensionality of the coordinate
    }

    inline void print_parameters(std::ostream &out, const parameters &params) {
        out << "# ===== ANIKE Parameters =====" << std::endl;
        out << "# stat_sec_param: " << params.stat_sec_param << std::endl;
        out << "# sec_param: " << params.sec_param << std::endl;
        out << "# n_width: " << params.n_width << std::endl;
        out << "# D (dimensionality of coordinate): " << params.D << std::endl;
        out << "# L (number of bits to represent a coordinate): " << params.l << std::endl;
        out << "# ============================" << std::endl << std::endl;
    }

    struct crs {
        parameters params;
        mkhss::crs mkhss_crs; // MKHSS CRS
        nike::crs nike_crs; // NIKE CRS
    };

    inline void setup(crs &crs, const parameters &params, mpz_srcptr n = nullptr) {
        mkhss::parameters mkhss_params;
        crs.params = params;
        mkhss::create_parameters(mkhss_params,
                                   true, // short_exponent_assumption
                                   params.stat_sec_param,
                                   params.sec_param,
                                   params.n_width,
                                   1); // B = 2 to support bits
        mkhss::setup(crs.mkhss_crs, mkhss_params, n);
        nike::setup(crs.nike_crs);
    }

    inline void check_integer(const crs &crs, coord_t x) {
        // Check if x is a valid coordinate
        // std::cout << "Checking coordinate: " << x << std::endl;
        // std::cout << "VS L: " << (((coord_t) 1) << crs.params.l) << std::endl;
        assert(crs.params.l == sizeof(coord_t) * 8 || x < (((coord_t) 1) << ((coord_t) crs.params.l))); // Ensure x fits in l bits
    }

    // Appends the binary decomposition of x into the public and private shares
    inline void append_bits(const crs &crs, const mkhss::public_key &pk,
                            std::vector<mkhss::public_share> &pe, 
                            std::vector<mkhss::private_share> &st, coord_t x) {
        // assert(crs.params.l == sizeof(coord_t) * 8 || x < (((coord_t) 1) << crs.params.l)); // Ensure x fits in l bits
        check_integer(crs, x);
        FMPZ x_i;
        for (int i = 0; i < crs.params.l; i++) {
            x_i.set((x >> (crs.params.l - 1 - i)) & 1); // Extract the (L-i)-th bit of x (MSB first)
            pe.emplace_back();
            st.emplace_back();
            mkhss::public_share &x_i_pk = pe.back();
            mkhss::private_share &x_i_sk = st.back();
            mkhss::share(crs.mkhss_crs, pk, x_i, x_i_sk, x_i_pk);
        }
    }

    struct public_encoding { // pe_sigma
        mkhss::public_key mkhss_pk; // MKHSS public key
        std::vector<mkhss::public_share> public_shares_alice;
        std::vector<mkhss::public_share> public_shares_bob;

        nike::public_key nike_pk; // NIKE public key
    };

    struct private_state { // st_sigma
        mkhss::private_key mkhss_sk; // MKHSS private key
        std::vector<mkhss::private_share> private_shares_alice;
        std::vector<mkhss::private_share> private_shares_bob;

        nike::private_key nike_sk; // NIKE private key
    };

    struct synced_shares_alice {
        // std::vector<mkhss::input_share> u;
        // std::vector<mkhss::input_share> v;
        // Refactor to allow arbitrary number of coordinates
        std::vector<std::vector<mkhss::input_share>> x; // coords
    };

    struct synced_shares_bob {
        // std::vector<mkhss::input_share> u_minus_d; // u - d
        // std::vector<mkhss::input_share> u_plus_d; // u + d
        // std::vector<mkhss::input_share> v_minus_d; // v - d
        // std::vector<mkhss::input_share> v_plus_d; // v + d
        std::vector<std::pair<std::vector<mkhss::input_share>, std::vector<mkhss::input_share>>> x; // (x_i - d, x_i + d) for each coordinate
    };

    inline void get_pointers_alice(const crs &crs, synced_shares_alice &shares, std::vector<mkhss::input_share *> &pointers) {
        size_t expected_size = crs.params.D * crs.params.l; // D shares for Alice

        pointers.clear();
        pointers.reserve(expected_size);

        // shares.u.resize(crs.params.l);
        // for (auto &share : shares.u) {
        //     pointers.push_back(&share);
        // }
        // shares.v.resize(crs.params.l);
        // for (auto &share : shares.v) {
        //     pointers.push_back(&share);
        // }
        for (int i = 0; i < crs.params.D; i++) {
            shares.x.emplace_back();
            shares.x.back().resize(crs.params.l);
            for (auto &share : shares.x.back()) {
                pointers.push_back(&share);
            }
        }

        assert(pointers.size() == expected_size);
    }

    inline void get_pointers_bob(const crs &crs, synced_shares_bob &shares, std::vector<mkhss::input_share *> &pointers) {
        size_t expected_size = 2 * crs.params.D * crs.params.l; // 2 * D shares for Bob

        pointers.clear();
        pointers.reserve(expected_size);

        // shares.u_minus_d.resize(crs.params.l);
        // for (auto &share : shares.u_minus_d) {
        //     pointers.push_back(&share);
        // }
        // shares.u_plus_d.resize(crs.params.l);
        // for (auto &share : shares.u_plus_d) {
        //     pointers.push_back(&share);
        // }
        // shares.v_minus_d.resize(crs.params.l);
        // for (auto &share : shares.v_minus_d) {
        //     pointers.push_back(&share);
        // }
        // shares.v_plus_d.resize(crs.params.l);
        // for (auto &share : shares.v_plus_d) {
        //     pointers.push_back(&share);
        // }
        for (int i = 0; i < crs.params.D; i++) {
            shares.x.emplace_back();
            shares.x.back().first.resize(crs.params.l); // x_i - d
            shares.x.back().second.resize(crs.params.l); // x_i + d
            for (auto &share : shares.x.back().first) {
                pointers.push_back(&share);
            }
            for (auto &share : shares.x.back().second) {
                pointers.push_back(&share);
            }
        }

        assert(pointers.size() == expected_size);
    }

    struct keygen_profile {
        double keygen_time; // Time taken for MKHSS::KeyGen
        double precompute_time; // Time taken for precomputation
        double share_time; // Time taken for calling MKHSS::Share
        double nike_time = 0.0; // Time taken for NIKE key generation

        inline nlohmann::json to_json() const {
            nlohmann::json j;
            j["keygen_time"] = keygen_time;
            j["precompute_time"] = precompute_time;
            j["share_time"] = share_time;
            j["nike_time"] = nike_time;
            return j;
        }
    };

    // d only needs to be specified by one party, here it will be Bob
    inline void attr_key_gen(const crs &crs, std::optional<party_id> sigma, const std::vector<coord_t> &x, const std::vector<coord_t> &d, public_encoding &pe, private_state &st, 
                             keygen_profile &profile) {
        profile.keygen_time = measure_time_seconds([&]() {
            mkhss::keygen(crs.mkhss_crs, pe.mkhss_pk, st.mkhss_sk);
        });

        // sigma = \bot means this party need to act as both Alice and Bob
        bool play_alice = !sigma.has_value() || *sigma == party_id::ALICE;
        bool play_bob = !sigma.has_value() || *sigma == party_id::BOB;

        int num_shares = play_alice * (crs.params.l * 2) + play_bob * (crs.params.l * 4);
        profile.precompute_time = measure_time_seconds([&]() {
            mkhss::precompute_share(crs.mkhss_crs, pe.mkhss_pk, num_shares); // Precompute shares for Alice and/or Bob
        });

        profile.share_time = measure_time_seconds([&]() {
            if (play_alice) {
                // Alice only needs u and v
                // append_bits(crs, pe.mkhss_pk, pe.public_shares_alice, st.private_shares_alice, u);
                // append_bits(crs, pe.mkhss_pk, pe.public_shares_alice, st.private_shares_alice, v);
                for (int i = 0; i < crs.params.D; i++) {
                    append_bits(crs, pe.mkhss_pk, pe.public_shares_alice, st.private_shares_alice, x[i]);
                }
            }
            if (play_bob) {
                // Bob needs u - d, u + d, v - d, v + d
                assert(d.size() == crs.params.D && "d is required for Bob");
                for (int i = 0; i < crs.params.D; i++) {
                    assert(x[i] >= d[i] && "x[i] must be greater than or equal to d[i]");
                    append_bits(crs, pe.mkhss_pk, pe.public_shares_bob, st.private_shares_bob, x[i] - d[i]);
                    // assert(x[i] + d[i] < (1u << crs.params.l) && "x[i] + d[i] must fit in l bits");
                    check_integer(crs, x[i] + d[i]);
                    append_bits(crs, pe.mkhss_pk, pe.public_shares_bob, st.private_shares_bob, x[i] + d[i]);
                }
                // assert(u >= d.value() && "u must be greater than or equal to d");
                // append_bits(crs, pe.mkhss_pk, pe.public_shares_bob, st.private_shares_bob, u - d.value());
                // append_bits(crs, pe.mkhss_pk, pe.public_shares_bob, st.private_shares_bob, u + d.value());
                // assert(v >= d.value() && "v must be greater than or equal to d");
                // append_bits(crs, pe.mkhss_pk, pe.public_shares_bob, st.private_shares_bob, v - d.value());
                // append_bits(crs, pe.mkhss_pk, pe.public_shares_bob, st.private_shares_bob, v + d.value());
            }
        });
        // NIKE
        profile.nike_time = measure_time_seconds([&]() {
            nike::keygen(crs.nike_crs, st.nike_sk, pe.nike_pk);
        });
    }

    // inline void attr_key_gen(const crs &crs, std::optional<party_id> sigma, coord_t u, coord_t v, std::optional<coord_t> d, public_encoding &pe, private_state &st) {
    //     assert(crs.params.D == 2 && "This function only supports 2D coordinates (e.g., latitude and longitude)");
    //     std::vector<coord_t> x(crs.params.D);
    //     x[0] = u; // Latitude
    //     x[1] = v; // Longitude
    //     std::vector<coord_t> d_vec(crs.params.D, 0);
    //     if (d.has_value()) {
    //         d_vec[0] = d.value(); // Latitude deviation
    //         d_vec[1] = d.value(); // Longitude deviation
    //     }
    //     attr_key_gen(crs, sigma, x, d_vec, pe, st);
    // }

    inline void less_than(const crs &crs, mkhss::rms_context &ek,
                            const std::vector<mkhss::input_share> &x, const std::vector<mkhss::input_share> &y,
                            const mkhss::memory_share &z_0, mkhss::memory_share &result) {
        mkhss::memory_share z;
        // Initialize z with the initial value
        mkhss::rms_mset(z, z_0);
        mkhss::memory_share lessthan_accumulator;
        mkhss::rms_mzero(lessthan_accumulator);
        for (int i = 0; i < crs.params.l; i++) {
            mkhss::memory_share x_z, xbar_z, xbar_y_z;
            // <<x z>> <- {{x}} * <<z>>
            ek.mult_serial(crs.mkhss_crs, x[i], z, x_z);
            // <<x_bar z>> <- <<z>> - {{x}} * <<z>>
            mkhss::rms_msub(z, x_z, xbar_z);
            // <<x_bar y z>> <- {{y}} * <<x_bar z>>
            ek.mult_serial(crs.mkhss_crs, y[i], xbar_z, xbar_y_z);

            mkhss::rms_madd(lessthan_accumulator, xbar_y_z, lessthan_accumulator);
            if (i != crs.params.l - 1) { // slight optimization: don't compute the equality value for the last bit
                mkhss::memory_share x_y_z, xbar_ybar_z;
                // <<x_bar y_bar z>> <- <<x_bar z>> - <<x_bar y z>>
                mkhss::rms_msub(xbar_z, xbar_y_z, xbar_ybar_z);
                // <<x y z>> <- {{y}} * <<x z>>
                ek.mult_serial(crs.mkhss_crs, y[i], x_z, x_y_z);
                // <<z>> <- <<x_bar y_bar z>> + <<x y z>>
                mkhss::rms_madd(xbar_ybar_z, x_y_z, z);
            }
        }
        mkhss::rms_mset(result, lessthan_accumulator);
    }

    inline void less_than_or_equal(const crs &crs, mkhss::rms_context &ek,
                            const std::vector<mkhss::input_share> &x, const std::vector<mkhss::input_share> &y,
                            const mkhss::memory_share &z_0, mkhss::memory_share &result) {
        mkhss::memory_share z;
        // Initialize z with the initial value
        mkhss::rms_mset(z, z_0);
        mkhss::memory_share lessthan_accumulator;
        mkhss::rms_mzero(lessthan_accumulator);
        for (int i = 0; i < crs.params.l; i++) {
            mkhss::memory_share x_z, xbar_z, xbar_y_z;
            // <<x z>> <- {{x}} * <<z>>
            ek.mult_serial(crs.mkhss_crs, x[i], z, x_z);
            // <<x_bar z>> <- <<z>> - {{x}} * <<z>>
            mkhss::rms_msub(z, x_z, xbar_z);
            // <<x_bar y z>> <- {{y}} * <<x_bar z>>
            ek.mult_serial(crs.mkhss_crs, y[i], xbar_z, xbar_y_z);

            mkhss::rms_madd(lessthan_accumulator, xbar_y_z, lessthan_accumulator);

            mkhss::memory_share x_y_z, xbar_ybar_z;
            // <<x_bar y_bar z>> <- <<x_bar z>> - <<x_bar y z>>
            mkhss::rms_msub(xbar_z, xbar_y_z, xbar_ybar_z);
            // <<x y z>> <- {{y}} * <<x z>>
            ek.mult_serial(crs.mkhss_crs, y[i], x_z, x_y_z);
            // <<z>> <- <<x_bar y_bar z>> + <<x y z>>
            mkhss::rms_madd(xbar_ybar_z, x_y_z, z);
        }
        mkhss::rms_madd(z, lessthan_accumulator, result); // result = \mathbb{1}(x = y) + \mathbb{1}(x < y)
    }

    inline void greater_than(const crs &crs, mkhss::rms_context &ek,
                            const std::vector<mkhss::input_share> &x, const std::vector<mkhss::input_share> &y,
                            const mkhss::memory_share &z_0, mkhss::memory_share &result) {
        mkhss::memory_share z;
        less_than_or_equal(crs, ek, x, y, z_0, z); // z = \mathbb{1}(x <= y) * z_0
        mkhss::rms_msub(z_0, z, result); // result = z_0 - z = \mathbb{1}(x > y) * z_0
    }

    struct keyder_profile {
        double rms_init_time; // Time taken for RMSInit
        mkhss::sync_shares_profile sync_shares_profile; // Time taken for MKHSS::SyncShares
        double rms_mult; // Time taken for RMS multiplications
        double rms_precompute; // Time taken for RMS precomputation
        double memory_management; // Time taken for memory management (e.g., clearing precomputed data)
        double key_derivation_time = 0.0; // Time taken for key derivation + hashing
        uint64_t num_rms_multiplications; // Number of RMS multiplications performed

        inline nlohmann::json to_json() const {
            nlohmann::json j;
            j["rms_init_time"] = rms_init_time;
            j["sync_shares_profile"] = sync_shares_profile.to_json();
            j["rms_mult"] = rms_mult;
            j["rms_precompute"] = rms_precompute;
            j["memory_management"] = memory_management;
            j["key_derivation_time"] = key_derivation_time;

            j["num_rms_muls"] = num_rms_multiplications;
            return j;
        }
    };

    inline void attr_key_der(const crs &crs, party_id sigma,
                                    const private_state &st_self, const public_encoding &pe_other,
                                    uint8_t key[KEY_LENGTH], keyder_profile& profile) {
        mkhss::rms_context ctx;
        profile.rms_init_time = measure_time_seconds([&]() {
            mkhss::rms_bootstrap(crs.mkhss_crs, sigma, st_self.mkhss_sk, pe_other.mkhss_pk, ctx);
        });
        
        synced_shares_alice inputs_A;
        synced_shares_bob inputs_B;
        std::vector<mkhss::input_share *> inputs_A_pointers, inputs_B_pointers;
        get_pointers_alice(crs, inputs_A, inputs_A_pointers);
        get_pointers_bob(crs, inputs_B, inputs_B_pointers);
        profile.sync_shares_profile = mkhss::sync_shares(crs.mkhss_crs, sigma, ctx, st_self.mkhss_sk, pe_other.mkhss_pk,
                           sigma == party_id::ALICE ? st_self.private_shares_alice : st_self.private_shares_bob,
                           inputs_A_pointers,
                           sigma == party_id::ALICE ? pe_other.public_shares_bob : pe_other.public_shares_alice,
                           inputs_B_pointers);
        
        mkhss::memory_share z;
        mkhss::rms_mset(z, ctx.get_one());

        // for (int i = 0; i < crs.params.l; i++) {
        //     for (int j = 0; j < crs.params.D; j++) {
        //         // Precompute shares for Alice and Bob
        //         inputs_A.x[j][i].precompute(crs.mkhss_crs, 4, bitcnt_n(crs.params.l) + 2);
        //     }
        //     // inputs_A.u[i].precompute(crs.mkhss_crs, 4, bitcnt_n(crs.params.l * 6) + 1);
        //     // inputs_A.v[i].precompute(crs.mkhss_crs, 4, bitcnt_n(crs.params.l * 6) + 1);
        // }

        double rms_mult_total = 0.0;
        double rms_precompute_total = 0.0;
        double memory_magagement_total = 0.0;

        // // u_B - d < u_A < u_B + d
        // less_than(crs, ctx, inputs_B.u_minus_d, inputs_A.u, z, z);
        // // less_than(crs, ctx, inputs_A.u, inputs_B.u_plus_d, z, z);
        // greater_than(crs, ctx, inputs_B.u_plus_d, inputs_A.u, z, z); // Consistently using u as the right hand size (y) helps with precomputation
        // // v_B - d < v_A < v_B + d
        // less_than(crs, ctx, inputs_B.v_minus_d, inputs_A.v, z, z);
        // // less_than(crs, ctx, inputs_A.v, inputs_B.v_plus_d, z, z);
        // greater_than(crs, ctx, inputs_B.v_plus_d, inputs_A.v, z, z);
        for (int j = 0; j < crs.params.D; j++) {
            rms_precompute_total += measure_time_seconds([&] {
                for (int i = 0; i < crs.params.l; i++) {
                    // Precompute shares for Alice and Bob
                    inputs_A.x[j][i].precompute(crs.mkhss_crs, 4, bitcnt_n(crs.params.l) + 2);
                }
            });
            rms_mult_total += measure_time_seconds([&] {
                // x_j_B - d < x_j_A < x_j_B + d
                less_than(crs, ctx, inputs_B.x[j].first, inputs_A.x[j], z, z); // x_j_B - d < x_j_A
                greater_than(crs, ctx, inputs_B.x[j].second, inputs_A.x[j], z, z); // x_j_B + d > x_j_A
            });
            memory_magagement_total += measure_time_seconds([&] {
                // Uncompute to save memory
                for (int i = 0; i < crs.params.l; i++) {
                    inputs_A.x[j][i].clear_precomputed_data();
                }
            });
        }

        // z is now a logical conjunction of all the conditions

        // Invert z (z = 1: polic satisfied, z = 0: not satisfied)
        mkhss::rms_msub(ctx.get_one(), z, z);

        // for (int i = 0; i < crs.params.l; i++) {
        //     std::cout << "usage count u: " << inputs_A.u[i].usage_count++;
        //     std::cout << ", v: " << inputs_A.v[i].usage_count++ << std::endl;
        // }

        profile.key_derivation_time = measure_time_seconds([&]() {
            uint8_t k_nike[nike::SHARED_KEY_LEN];
            nike::keyder(crs.nike_crs, st_self.nike_sk, pe_other.nike_pk, k_nike);
            derive_key(k_nike, z, key);
        });
        
        profile.rms_mult = rms_mult_total;
        profile.rms_precompute = rms_precompute_total;
        profile.num_rms_multiplications = ctx.get_next_inst_id();
        profile.memory_management = memory_magagement_total;
    }
}

#undef mkhss

#endif
