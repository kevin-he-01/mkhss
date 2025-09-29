#ifndef BASELINE_MKHSS_HH
#define BASELINE_MKHSS_HH

#include <gmpxx.h>
#include <optional>
#include <vector>
#include <nlohmann/json.hpp>

#include "util.hh"
#include "prf.hh"
#include "fmpz_class.hh"
#include "fmpz_mod.hh"
#include "nike.hh"

#include "fixed_base_exp.hh"
using namespace fixed_base_exp;

namespace baseline_mkhss {
    constexpr bool IS_BASELINE = true; // Set to true if using the baseline MKHSS scheme

    enum party_id {
        ALICE, BOB
    };

    constexpr int NO_MOD_REDUCE = -1;

    struct parameters {
        bool short_exponent_assumption;
        int stat_sec_param; // Tau. Statistical security parameter (for correctness)
        int sec_param; // Lambda. Computational security parameter (against offline attacks)
        int n_width; // When lambda >= 128, need n_width >= 3072
        int log2_b; // log_2(B), where B is the bound on the absolute value of all integers in memory shares
        int w; // Ciphertexts modulus is N^{w+1}. Plaintext modulus is N^w
        // depenent variables below
        int prf_key_length; // AES key length
        int secret_bits; // The width of t. 2 * lambda is needed to be secure against Pollard's rho.
        int bprime_bits; // log_2(B'), B' = B * 2^{tau}
        int exponent_bits; // The width of the exponent (-1 if no mod reduction should be used)
    };

    inline void print_parameters(const parameters &params, std::ostream &out = std::cout) {
        out << "# ===== Baseline MKHSS Parameters =====" << std::endl;
        out << "# short_exponent_assumption: " << (params.short_exponent_assumption ? "true" : "false") << std::endl;
        out << "# stat_sec_param: " << params.stat_sec_param << std::endl;
        out << "# sec_param: " << params.sec_param << std::endl;
        out << "# n_width: " << params.n_width << std::endl;
        out << "# log2_b: " << params.log2_b << std::endl;
        out << "# w: " << params.w << std::endl;
        out << "# ----- Dependent variables" << std::endl;
        out << "# prf_key_length: " << params.prf_key_length << std::endl;
        out << "# secret_bits: " << params.secret_bits << std::endl;
        out << "# bprime_bits: " << params.bprime_bits << std::endl;
        out << "# exponent_bits: " << params.exponent_bits << std::endl;
        out << "# ===== End Parameters =====" << std::endl;
    }

    inline void create_parameters(parameters &params,
                                  bool short_exponent_assumption,
                                  int stat_sec_param,
                                  int sec_param,
                                  int n_width,
                                  int log2_b,
                                  int w = 3) {
        assert(short_exponent_assumption); // Non SEI not supported

        params.short_exponent_assumption = short_exponent_assumption;
        params.stat_sec_param = stat_sec_param;
        params.sec_param = sec_param;
        params.n_width = n_width;
        params.log2_b = log2_b;
        params.w = w; // TODO: choose w to optimize performance while fitting the exponent

        // Set dependent variables
        params.prf_key_length = get_prf_key_length(sec_param); // AES key length
        params.secret_bits = short_exponent_assumption ? (2 * sec_param) : n_width;

        // assert(log2_b + stat_sec_param <= params.n_width);
        if (log2_b + stat_sec_param > params.n_width) {
            std::cerr << "log2_b + stat_sec_param (" << log2_b + stat_sec_param
                      << ") exceeds n_width (" << params.n_width << "). "
                      << "This will lead to exponent overflow." << std::endl;
            assert(false);
        }
        params.bprime_bits = params.n_width; // B' = N
        // Secret: s_sigma = t_sigma * N + 1, joint secret: s = s_A * s_B
        params.exponent_bits = 2 * (params.bprime_bits + params.secret_bits) + log2_b + stat_sec_param;
        // assert(params.exponent_bits <= (params.n_width * params.w)); // exponent <= plaintext modulus
        if (params.exponent_bits > (params.n_width * params.w)) {
            std::cerr << "Exponent length " << params.exponent_bits
                      << " exceeds plaintext modulus, which is " << (params.n_width * params.w) << " bits" << std::endl;
            assert(false);
        }
    }

    struct crs {
        parameters params;
        FMPZ n;
        FMPZ g;
        FMPZ h;
        // nike::crs nike_crs; // NIKE CRS
        // Things to speed up computation
        FMPZ n_2; // n^2
        FMPZ n_w; // n^w, plaintext modulus
        FMPZ n_w_1; // n^(w+1), ciphertext modulus
        FMPZModCtx n_2_ctx; // Context for n^2
        FMPZModCtx n_w_1_ctx; // Context for n^(w+1)
    };

    struct public_key {
        // DJEG public key
        FMPZ f;
        // NIM public encoding
        FMPZ nim_pe_A;
        FMPZ nim_pe0_B;
        FMPZ nim_pe1_B;
        // NIKE public key
        // nike::public_key nike_pk;
        // Since precomputation is not an algorithmic optimization (and known in prior work), we apply it to the baseline
        FixedBaseExp g_fbe_n_w_1; // Precomputed FBE for g (mod N^{w+1})
        FixedBaseExp g_fbe_n_2; // Precomputed FBE for g (mod N^2)
        FixedBaseExp f_fbe_n_w_1; // Precomputed FBE for f (mod N^{w+1})
        FixedBaseExp f_fbe_n_2; // Precomputed FBE for f (mod N^2)
        bool precomputed = false; // Whether the FBE is precomputed
    };

    // num_shares is the number of invocations to MKHSS.Share expected.
    inline void precompute_share(const crs &crs, public_key &pk, int num_shares) {
        // Precompute FBE for g and f
        size_t window_size = get_optimal_window_size(crs.params.sec_param * 2, num_shares);
        if (window_size == 0) {
            // If window size is 0, it means we should not precompute
            return;
        }
        pk.g_fbe_n_w_1.precompute(crs.g, crs.n_w_1, crs.params.sec_param * 2, window_size);
        pk.f_fbe_n_w_1.precompute(pk.f, crs.n_w_1, crs.params.sec_param * 2, window_size);
        // pk.g_fbe_n_2.precompute(crs.g, crs.n_2, crs.params.sec_param * 2, window_size);
        // pk.f_fbe_n_2.precompute(pk.f, crs.n_2, crs.params.sec_param * 2, window_size);
        pk.g_fbe_n_2.downcast(pk.g_fbe_n_w_1, crs.n_2);
        pk.f_fbe_n_2.downcast(pk.f_fbe_n_w_1, crs.n_2);
        pk.precomputed = true;
    }

    inline size_t get_public_key_size(const crs &crs, const public_key &pk) {
        assert(pk.f.cmp(0) >= 0 && pk.f.cmp(crs.n_w_1) < 0);
        assert(pk.nim_pe_A.cmp(0) >= 0 && pk.nim_pe_A.cmp(crs.n_w_1) < 0);
        assert(pk.nim_pe0_B.cmp(0) >= 0 && pk.nim_pe0_B.cmp(crs.n_w_1) < 0);
        assert(pk.nim_pe1_B.cmp(0) >= 0 && pk.nim_pe1_B.cmp(crs.n_w_1) < 0);

        // Leading zeros might shrink the empirical size
        size_t empirical_size = pk.f.get_size_in_bytes() + pk.nim_pe_A.get_size_in_bytes()
                                + pk.nim_pe0_B.get_size_in_bytes() + pk.nim_pe1_B.get_size_in_bytes();
                                // + pk.nike_pk.point.get_serialization_size(crs.nike_crs.group);

        // Serialization must have well-marked boundaries between fields, so need to use values from CRS
        size_t serialized_size = 4 * crs.n_w_1.get_size_in_bytes(); //+ pk.nike_pk.point.get_serialization_size(crs.nike_crs.group);
        assert(empirical_size <= serialized_size);
        // std::cerr << "DEBUG: empirical size: " << empirical_size << ", serialized size: " << serialized_size << std::endl;

        return serialized_size;
    }

    struct private_key {
        // DJEG private key (also used in NIM.Decode for Alice)
        FMPZ s;
        // NIM private state
        FMPZ r_A;
        FMPZ r_B;
        // NIKE private key
        // nike::private_key nike_sk;
    };

    // Secret share intended for the sharing party itself
    struct private_share {
        FMPZ x;
        FMPZ r;
        FMPZ r_prime;
        FMPZ c0;
        FMPZ c0_prime;
    };

    // Secret share intended for another party
    struct public_share {
        FMPZ c0;
        FMPZ c1;
        FMPZ c0_prime;
        FMPZ c1_prime;
    };

    inline size_t get_public_share_size(const crs &crs) {
        return 2 * crs.n_w_1.get_size_in_bytes() + 2 * crs.n_2.get_size_in_bytes();
    }

    inline size_t get_public_share_size(const crs &crs, const public_share &x) {
        assert(x.c0.cmp(0) >= 0 && x.c0.cmp(crs.n_w_1) < 0);
        assert(x.c1.cmp(0) >= 0 && x.c1.cmp(crs.n_w_1) < 0);
        assert(x.c0_prime.cmp(0) >= 0 && x.c0_prime.cmp(crs.n_2) < 0);
        assert(x.c1_prime.cmp(0) >= 0 && x.c1_prime.cmp(crs.n_2) < 0);

        // Leading zeros might shrink the empirical size
        size_t empirical_size = x.c0.get_size_in_bytes() + x.c1.get_size_in_bytes() + 
                                x.c0_prime.get_size_in_bytes() + x.c1_prime.get_size_in_bytes();
        
        // Serialization must have well-marked boundaries between fields, so need to use values from CRS
        size_t serialized_size = 2 * crs.n_w_1.get_size_in_bytes() + 2 * crs.n_2.get_size_in_bytes();
        assert(empirical_size <= serialized_size);
        // std::cerr << "DEBUG: empirical size: " << empirical_size << ", serialized size: " << serialized_size << std::endl;
        return serialized_size;
    }

    // Synchronized input share. Essentially a DJEG flipped ciphertext encrypting x * s
    struct input_share {
        FMPZ c0;
        FMPZ c1;
        FMPZ c0_prime;
        FMPZ c1_prime;
        mutable int usage_count = 0;
        FixedBaseExp c0_fbe_n_w_1; // Precomputed FBE for c0 (mod N^{w+1})
        FixedBaseExp c1_fbe_n_w_1; // Precomputed FBE for c1 (mod N^{w+1})
        FixedBaseExp c0_prime_fbe_n_2; // Precomputed FBE for c0 (mod N^2)
        FixedBaseExp c1_prime_fbe_n_2; // Precomputed FBE for c1 (mod N^2)
        bool precomputed = false; // Whether the FBE is precomputed

        inline void precompute(const crs &crs, int expected_usage_count, int leeway = 0) {
            // leeway is to account for adding multiple memory shares leading to larger exponents
            // It should be roughly log_2(# of max memory share additions/subtractions)
            if (precomputed) {
                std::cerr << "WARNING: Input share already precomputed (possible bug)." << std::endl;
            }
            if (expected_usage_count <= 1) {
                // No need to precompute
                return;
            }
            size_t size_n_w_bits = fmpz_sizeinbase(crs.n_w.get_fmpz(), 2);
            size_t size_n_bits = fmpz_sizeinbase(crs.n.get_fmpz(), 2);
            // std::cout << "DEBUG: size_n_w_bits: " << size_n_w_bits << ", size_n_bits: " << size_n_bits << std::endl;
            c0_fbe_n_w_1.precompute(c0, crs.n_w_1, size_n_w_bits + leeway, get_optimal_window_size(size_n_bits * crs.params.w, expected_usage_count)); // <<y*s>> has size log_2(N^w) bits
            c0_prime_fbe_n_2.precompute(c0, crs.n_2, size_n_w_bits + leeway, get_optimal_window_size(size_n_bits * crs.params.w, expected_usage_count)); // <<y*s>> has size log_2(N^w) bits
            c1_fbe_n_w_1.precompute(c1, crs.n_w_1, size_n_bits + leeway, get_optimal_window_size(size_n_bits, expected_usage_count)); // <<y>> has size log_2(N) bits
            c1_prime_fbe_n_2.precompute(c1, crs.n_2, size_n_bits + leeway, get_optimal_window_size(size_n_bits, expected_usage_count)); // <<y>> has size log_2(N) bits
            precomputed = true;
        }

        inline void clear_precomputed_data() {
            // Clear precomputed data
            c0_fbe_n_w_1.clear();
            c1_fbe_n_w_1.clear();
            c0_prime_fbe_n_2.clear();
            c1_prime_fbe_n_2.clear();
            // Reset precomputed flag
            precomputed = false;
        }
    };

    struct memory_share {
        FMPZ y_s; // <y*s>_sigma
        FMPZ y; // <y>_sigma
    };

    // [[deprecated]]
    // void rms_output(const memory_share &mem_share, FMPZ &output);
    void rms_mult(const crs &crs, PRF &prf, uint64_t unique_id, const input_share &b, const memory_share &a, memory_share &result);

    class rms_context { // This is the "evaluation key" for MKHSS, i.e., ek_sigma in the paper
    private:
        // const crs &crs;
        PRF prf_rms;
        PRF prf_nim;
        PRF prf_shared; // A general purpose shared key between two parties. Can use this instead of performing another NIKE in parallel. Used by ANIKE to simulate id_A || id_B
        uint64_t next_inst_id;
        bool initialized = false;

        memory_share one;
        bool one_set = false;

    public:
        FMPZ f; // g^{-s}, precomputed to speed up ExpLinEncS
        FixedBaseExp f_fbe_n_w_1; // Precomputed FBE for f (mod N^{w+1})
        FixedBaseExp f_fbe_n_2; // Precomputed FBE for f (mod N^2)
        bool precomputed = false; // Whether the FBE is precomputed

        inline void init(const crs &crs, party_id sigma, const private_key &sk_self, const public_key &pk_other, uint64_t program_description = 0) {
            // PRF prf(crs.prf_key, crs.params.prf_key_length);
            uint8_t shared_secret[nike::SHARED_KEY_LEN];
            // nike::keyder(crs.nike_crs, sk_self.nike_sk, pk_other.nike_pk, shared_secret);
            assert(crs.params.prf_key_length <= nike::SHARED_KEY_LEN);

            f.powm(pk_other.f, sk_self.s, crs.n_w_1);

            // Hash f into shared_secret
            size_t f_size = f.get_size_in_bytes();
            uint8_t *buf = new uint8_t[f_size + 1];
            f.export_(buf, f_size);
            buf[f_size] = random_oracle::MAGIC_MKHSS;
            random_oracle::call(buf, f_size + 1, shared_secret);

            PRF prf(shared_secret, crs.params.prf_key_length);
            PRF prf_program;

            prf.derive_subkey(program_description, prf_program);

            prf_program.derive_subkey(0, prf_rms);
            prf_program.derive_subkey(1, prf_nim);
            prf_program.derive_subkey(2, prf_shared);

            next_inst_id = 0;

            initialized = true;

            delete[] buf;
        }

        inline void precompute(const crs &crs, int explinencs_count) {
            assert(initialized);
            // Precompute FBE for f
            size_t window_size = get_optimal_window_size(crs.params.sec_param * 2, explinencs_count);
            if (window_size == 0) {
                // If window size is 0, it means we should not precompute
                return;
            }
            f_fbe_n_w_1.precompute(f, crs.n_w_1, crs.params.sec_param * 2, window_size);
            // f_fbe_n_2.precompute(f, crs.n_2, crs.params.sec_param * 2, window_size);
            f_fbe_n_2.downcast(f_fbe_n_w_1, crs.n_2); // roughly as slow as recomputation from scratch
            precomputed = true;
        }

        inline PRF &get_prf_nim() {
            assert(initialized);
            return prf_nim;
        }

        inline PRF &get_prf_shared() {
            assert(initialized);
            return prf_shared;
        }

        inline void set_one(const memory_share &one) {
            assert(initialized);
            this->one.y.set(one.y);
            this->one.y_s.set(one.y_s);
            one_set = true;
        }

        inline void mult(const crs &crs, uint64_t inst_id, const input_share &b, const memory_share &a, memory_share &result, std::optional<int> max_log_y = std::nullopt) {
            rms_mult(crs, prf_rms, inst_id, b, a, result);
        }

        // Not thread safe
        inline void mult_serial(const crs &crs, const input_share &b, const memory_share &a, memory_share &result, std::optional<int> max_log_y = std::nullopt) {
            mult(crs, next_inst_id++, b, a, result);
        }

        inline const memory_share &get_one() const {
            assert(one_set);
            return one;
        }

        inline void convert(const crs &crs, uint64_t inst_id, const input_share &b, memory_share &result, std::optional<int> max_log_y = std::nullopt) {
            rms_mult(crs, prf_rms, inst_id, b, get_one(), result);
        }

        // Not thread safe
        inline void convert_serial(const crs &crs, const input_share &b, memory_share &result, std::optional<int> max_log_y = std::nullopt) {
            convert(crs, next_inst_id++, b, result);
        }

        // Achieves external security by making each output uniformly random and independent
        inline void output(const crs &crs, uint64_t inst_id, const memory_share &mem_share, FMPZ &randomized_output) {
            FMPZ output;
            output.set(mem_share.y);
            // Randomize output
            FMPZ r;
            prf_rms.evaluate_mod_2exp(inst_id, crs.params.bprime_bits, r);
            randomized_output.add(output, r);
            randomized_output.mod_2exp(randomized_output, crs.params.bprime_bits);
        }

        inline void output_serial(const crs &crs, const memory_share &mem_share, FMPZ &randomized_output) {
            output(crs, next_inst_id++, mem_share, randomized_output);
        }

        inline uint64_t get_next_inst_id() {
            assert(initialized);
            return next_inst_id;
        }
    };

    void generate_generator(const crs &crs, FMPZ & g);
    void setup(crs &crs, const parameters &params, mpz_srcptr n = nullptr);
    void keygen(const crs &crs, public_key &pk, private_key &sk);
    void share(const crs &crs, const public_key &pk, const FMPZ & x, private_share &x_self, public_share &x_other);
    void sync_share_self(const crs &crs, const private_key &sk_self, const public_key &pk_other, const private_share &x_self, input_share &x_synced);
    void sync_share_self(const crs &crs, const rms_context &ctx, const private_key &sk_self, const public_key &pk_other, const private_share &x_self, input_share &x_synced);

    void sync_share_other(const crs &crs, const private_key &sk_self, const public_key &pk_other, const public_share &x_other, input_share &x_synced);
    // void nim_decode(const crs &crs, const private_key &sk_self, const public_key &pk_other, FMPZ &product_ss);
    void rms_bootstrap(const crs &crs, party_id sigma, const private_key &sk_self, const public_key &pk_other, rms_context &ctx, uint64_t program_description = 0);

    void rms_mzero(memory_share &zero);
    inline void rms_mset(memory_share &dst, const memory_share &src) {
        dst.y.set(src.y);
        dst.y_s.set(src.y_s);
    }
    void rms_madd(const memory_share &a, const memory_share &b, memory_share &result);
    void rms_msub(const memory_share &a, const memory_share &b, memory_share &result);
    void rms_mcmult(const memory_share &a, const FMPZ &b, memory_share &result);
    void rms_iadd(const crs &crs, const input_share &a, const input_share &b, input_share &result);
    void rms_isub(const crs &crs, const input_share &a, const input_share &b, input_share &result);
    void rms_icmult(const crs &crs, const input_share &a, const FMPZ &b, input_share &result);

    struct sync_shares_profile {
        double sync_self_precompute = 0.0; // Time spent on precomputing FBE for self shares
        double sync_self_explinencs = 0.0; // Time spent on syncing self shares (Excluding precompuption)
        double sync_other_time = 0.0; // Time spent on syncing other shares

        inline nlohmann::json to_json() const {
            nlohmann::json j;
            j["sync_self_precompute"] = sync_self_precompute;
            j["sync_self_explinencs"] = sync_self_explinencs;
            j["sync_other_time"] = sync_other_time;
            return j;
        }
    };

    inline sync_shares_profile sync_shares(const crs &crs, party_id sigma, rms_context &ctx, const private_key &sk_self, const public_key &pk_other,
                             const std::vector<private_share> &x_self, std::vector<input_share *> &x_A_synced,
                             const std::vector<public_share> &x_other, std::vector<input_share *> &x_B_synced) {
        std::vector<input_share *> &x_self_synced = sigma == party_id::ALICE ? x_A_synced : x_B_synced;
        std::vector<input_share *> &x_other_synced = sigma == party_id::ALICE ? x_B_synced : x_A_synced;
        sync_shares_profile profile;
        // std::cout << "DEBUG: x_self.size() = " << x_self.size() << ", x_self_synced.size() = " << x_self_synced.size() << std::endl;
        // std::cout << "DEBUG: x_other.size() = " << x_other.size() << ", x_other_synced.size() = " << x_other_synced.size() << std::endl;
        assert(x_self.size() == x_self_synced.size());
        assert(x_other.size() == x_other_synced.size());
        profile.sync_self_precompute = measure_time_seconds([&]() {
            ctx.precompute(crs, x_self_synced.size()); // Precompute FBE for f
        });
        profile.sync_self_explinencs = measure_time_seconds([&]() {
            for (size_t i = 0; i < x_self_synced.size(); i++) {
                assert(x_self_synced[i] != nullptr);
                sync_share_self(crs, ctx, sk_self, pk_other, x_self[i], *x_self_synced[i]);
            }
        });
        profile.sync_other_time = measure_time_seconds([&]() {
            for (size_t i = 0; i < x_other_synced.size(); i++) {
                assert(x_other_synced[i] != nullptr);
                sync_share_other(crs, sk_self, pk_other, x_other[i], *x_other_synced[i]);
            }
        });
        return profile;
    }
}

#endif
