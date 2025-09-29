#include <iostream>
#include <cassert>

#include "mkhss.hh"
#include "util.hh"


// This is the optimized construction of MKHSS introduced by this new paper
namespace mkhss {
    bool deprecate_warned = false;
    // This also allows hashing the transcript into a PRF key for malicious security

    party_id other_party(party_id sigma) {
        return (sigma == ALICE) ? BOB : ALICE;
    }

    inline void gen_random_exponent(const crs &crs, FMPZ &rnd) {
        // Generate a random exponent r, such that g^r is indistinguishable from a random group element
        // When using the short exponent assumption, sample from [0, 2^{2\lambda}-1], otherwise sample from [0, N-1]
        if (crs.params.short_exponent_assumption) {
            // Generate a random exponent in the range [0, 2^{2*\lambda}-1]
            gen_random_fmpz(rnd, 2 * crs.params.sec_param);
        } else {
            // Generate a random exponent in the range [0, N-1]
            gen_random_under_fmpz(rnd, crs.n);
        }
    }

    // Generates a generator of order phi(N) / 4 with high probability.
    // Writes generator to g
    void generate_generator(const crs &crs, FMPZ & g) {
        FMPZ g_0, two_n_w;
        gen_random_under_fmpz(g_0, crs.n_w_1);
        two_n_w.mul(crs.n_w, 2);
        g.powm(g_0, two_n_w, crs.n_w_1);
    }

    // Writes to the NIM private state and public encoding parts in sk and pk
    void nim_encode(const crs &crs, public_key &pk, private_key &sk) {
        FMPZ tmp;
        gen_random_exponent(crs, sk.r_A);
        gen_random_exponent(crs, sk.r_B);

        // pk.nim_pe_A = g^r_A * h^{t_A} mod N^{w+1}
        pk.nim_pe_A.powm(crs.g, sk.r_A, crs.n_w_1);
        tmp.powm(crs.h, sk.t, crs.n_w_1);
        pk.nim_pe_A.mul(pk.nim_pe_A, tmp);
        pk.nim_pe_A.mod(pk.nim_pe_A, crs.n_w_1);

        // pk.nim_pe0_B = g^{r_B} mod N^{w+1}
        // pk.nim_pe1_B = (N+1)^{t_B} * h^{r_B} mod N^{w+1}
        pk.nim_pe0_B.powm(crs.g, sk.r_B, crs.n_w_1);
        pow_base_n_1_mod_n_w_1(pk.nim_pe1_B, sk.t, crs.n, crs.params.w, crs.n_w_1);
        tmp.powm(crs.h, sk.r_B, crs.n_w_1);
        pk.nim_pe1_B.mul(pk.nim_pe1_B, tmp);
        pk.nim_pe1_B.mod(pk.nim_pe1_B, crs.n_w_1);
    }

    void setup(crs &crs, const parameters &params, mpz_srcptr n) {
        mpz_class p, q;

        crs.params = params;
        if (n == nullptr) {
            generate_safe_prime(p.get_mpz_t(), params.n_width / 2 + (params.n_width % 2));
            generate_safe_prime(q.get_mpz_t(), params.n_width / 2);
            crs.n.mul(p, q);
        } else {
            crs.n.set_mpz(n);
        }
        auto n_bits = fmpz_sizeinbase(crs.n.get_fmpz(), 2);
        assert(n_bits == params.n_width || n_bits == params.n_width - 1 && "N is of the wrong size");
        // Set N^w and N^{w+1}
        fmpz_pow_ui(crs.n_w.get_fmpz(), crs.n.get_fmpz(), params.w);
        fmpz_pow_ui(crs.n_w_1.get_fmpz(), crs.n.get_fmpz(), params.w + 1);
        crs.n_w_1_ctx.set_modulus(crs.n_w_1);
        // Generate random generators g and h
        generate_generator(crs, crs.g);
        generate_generator(crs, crs.h);
        
        // Setup NIKE (for external security)
        // nike::setup(crs.nike_crs);
    }

    void keygen(const crs &crs, public_key &pk, private_key &sk) {
        // Generate random s = t * B' + 1, where t is chosen depending on the assumption
        // First generate t, then multiply by B' and add 1
        if (crs.params.short_exponent_assumption) {
            gen_random_fmpz(sk.t, crs.params.secret_bits);
        } else {
            gen_random_under_fmpz(sk.t, crs.n);
        }

        // s = t * B' + 1
        fmpz_mul_2exp(sk.s.get_fmpz(), sk.t.get_fmpz(), crs.params.bprime_bits);
        sk.s.add(sk.s, 1);
        // f = g^{-s} mod N^{w+1}
        pk.f.powm(crs.g, sk.s, crs.n_w_1);
        // mpz_invert(pk.f, pk.f, crs.n_w_1);
        fmpz_invmod(pk.f.get_fmpz(), pk.f.get_fmpz(), crs.n_w_1.get_fmpz());
        // nim_encode with x = sk.s
        nim_encode(crs, pk, sk);

        // Generate NIKE keypair
        // nike::keygen(crs.nike_crs, sk.nike_sk, pk.nike_pk);
    }

    void keygen(const crs &crs, party_id sigma, public_key &pk, private_key &sk) {
        keygen(crs, pk, sk);
    }

    void share(const crs &crs, const public_key &pk, const FMPZ & x, private_share &x_self, public_share &x_other) {
        FMPZ tmp;
        // r <-$ Exponent distribution
        gen_random_exponent(crs, x_self.r);
        // c0 = (N+1)^x * g^r mod N^{w+1}
        // c1 = f^r mod N^{w+1}
        pow_base_n_1_mod_n_w_1(x_other.c0, x, crs.n, crs.params.w, crs.n_w_1);
        if (pk.precomputed) {
            // If the public key is precomputed, use the precomputed FBE
            pk.g_fbe.compute(x_self.r, tmp);
        } else {
            tmp.powm(crs.g, x_self.r, crs.n_w_1);
        }
        x_other.c0.mul(x_other.c0, tmp);
        x_other.c0.mod(x_other.c0, crs.n_w_1);
        if (pk.precomputed) {
            // If the public key is precomputed, use the precomputed FBE
            pk.f_fbe.compute(x_self.r, x_other.c1);
        } else {
            x_other.c1.powm(pk.f, x_self.r, crs.n_w_1);
        }
        x_self.c0.set(x_other.c0);
        // x_self.c1 = x_other.c1;
    }

    void sync_share_self(const crs &crs, const private_key &sk_self, const public_key &pk_other, const private_share &x_self, input_share &x_synced) {
        if (!deprecate_warned) {
            std::cerr << "[!] This version of ExpLinEncS is slow and does not take advantage of the precomputed shared secret g^{-s}!" << std::endl;
            deprecate_warned = true;
        }
        // ExpLinEncS (sender)
        x_synced.c0.set(x_self.c0);
        // c1 = f_{1-sigma}^(r*s_sigma), where f_{1-sigma} = g^{s_{1-sigma}} comes from the other party's public key
        x_synced.c1.powm(pk_other.f, x_self.r, crs.n_w_1);
        x_synced.c1.powm(x_synced.c1, sk_self.s, crs.n_w_1);
        // Below is an alternative way to compute c1. Roughly the same speed.
        // FMPZ r_s;
        // r_s.mul(x_self.r, sk_self.s);
        // x_synced.c1.powm(pk_other.f, r_s, crs.n_w_1);
    }

    void sync_share_self(const crs &crs, const rms_context &ctx, const private_key &sk_self, const public_key &pk_other, const private_share &x_self, input_share &x_synced) {
        // ExpLinEncS (sender)
        // Optimized version by precomputing shared secret and put it in the evaluation key, for example
        x_synced.c0.set(x_self.c0);
        // c1 = f^r, where f = g^{-s} = g^{-s_A s_B}
        if (ctx.precomputed) {
            // If the context is precomputed, use the precomputed FBE
            ctx.f_fbe.compute(x_self.r, x_synced.c1);
        } else {
            // Otherwise, compute it directly
            x_synced.c1.powm(ctx.f, x_self.r, crs.n_w_1);
        }
    }

    void sync_share_other(const crs &crs, const private_key &sk_self, const public_key &pk_other, const public_share &x_other, input_share &x_synced) {
        // ExpLinEncR (receiver)
        x_synced.c0.set(x_other.c0);
        // c1 = (c1_other)^s mod N^{w+1}
        x_synced.c1.powm(x_other.c1, sk_self.s, crs.n_w_1);
    }

    void nim_decode(const crs &crs, party_id sigma, const private_key &sk_self, const public_key &pk_other, FMPZ &product_ss) {
        FMPZ tmp;
        FMPZ z_sigma;
        switch(sigma) {
            case party_id::ALICE:
                // Alice: $z_A = \pe_{B,0}^{r_{\nim,A}}\pe_{B,1}^{t_A}$
                z_sigma.powm(pk_other.nim_pe0_B, sk_self.r_A, crs.n_w_1);
                tmp.powm(pk_other.nim_pe1_B, sk_self.t, crs.n_w_1);
                z_sigma.mul(z_sigma, tmp);
                z_sigma.mod(z_sigma, crs.n_w_1); 
                break;
            case party_id::BOB:
                // Bob: $z_B = \pe_A^{r_{\nim,B}$
                z_sigma.powm(pk_other.nim_pe_A, sk_self.r_B, crs.n_w_1);
                break;
        }
        // <product>_sigma = DDLog(z_sigma)
        ddlog(z_sigma, product_ss, crs.n, crs.params.w);
    }

    // InitRMS in the paper
    void rms_bootstrap(const crs &crs, party_id sigma, const private_key &sk_self, const public_key &pk_other, rms_context &ctx, uint64_t program_description) {
        FMPZ nim_rand_shift;
        
        ctx.init(crs, sigma, sk_self, pk_other, program_description);
        ctx.get_prf_nim().evaluate_mod_t(0, crs.n_w, crs.params.stat_sec_param, nim_rand_shift);
        
        FMPZ tmp;
        memory_share one;
        nim_decode(crs, sigma, sk_self, pk_other, tmp);
        // Add NIM random shift and mod N^{w}
        tmp.add(tmp, nim_rand_shift);
        if (tmp.cmp(crs.n_w) >= 0) {
            tmp.sub(tmp, crs.n_w);
        }
        reduce_exponent(tmp, tmp, crs.params.nim_product_bits);

        // tmp = z_\sigma = <t_A * t_B>_\sigma
        // y = 1
        // y * s = t_A * t_B * B'^2 + (t_A + t_B) * B' + 1
        tmp.mul_2exp(tmp, crs.params.bprime_bits * 2);
        
        switch(sigma) {
            case party_id::ALICE:
                // Alice:
                // <y>_A = 1
                // <y*s>_A = z_A * B'^2 + t_A * B' + 1
                one.y_s.mul_2exp(sk_self.t, crs.params.bprime_bits);
                one.y_s.add(one.y_s, tmp);
                one.y_s.add(one.y_s, 1);
                one.y.set(1);
                break;
            case party_id::BOB:
                // Bob: <y>_B = 0
                // <y*s>_B = z_B * B'^2 - t_B * B'
                one.y_s.set(tmp);
                tmp.mul_2exp(sk_self.t, crs.params.bprime_bits);
                one.y_s.sub(one.y_s, tmp);
                one.y.set(0);
                break;
        }

        ctx.set_one(one);
    }

    // Outputs an integer share
    void rms_output(const memory_share &mem_share, FMPZ &output) {
        output.set(mem_share.y);
    }

    void rms_mzero(memory_share &zero) {
        zero.y_s.set(0);
        zero.y.set(0);
    }

    // Add 2 memory shares. No modulo reduction is performed
    void rms_madd(const memory_share &a, const memory_share &b, memory_share &result) {
        result.y_s.add(a.y_s, b.y_s);
        result.y.add(a.y, b.y);
    }

    void rms_msub(const memory_share &a, const memory_share &b, memory_share &result) {
        result.y_s.sub(a.y_s, b.y_s);
        result.y.sub(a.y, b.y);
    }

    // Multiply memory shares by a constant. No modulo reduction is performed, so only use small b
    void rms_mcmult(const memory_share &a, const FMPZ &b, memory_share &result) {
        result.y_s.mul(a.y_s, b);
        result.y.mul(a.y, b);
    }

    // Add 2 input shares
    void rms_iadd(const crs &crs, const input_share &a, const input_share &b, input_share &result) {
        result.c0.mul(a.c0, b.c0);
        crs.n_w_1_ctx.mod(result.c0, result.c0); // Precomputed inverse to speed up mod reduction

        result.c1.mul(a.c1, b.c1);
        crs.n_w_1_ctx.mod(result.c1, result.c1);
    }

    // Subtract 2 input shares
    void rms_isub(const crs &crs, const input_share &a, const input_share &b, input_share &result) {
        crs.n_w_1_ctx.invmod(result.c0, b.c0);
        result.c0.mul(a.c0, result.c0);
        crs.n_w_1_ctx.mod(result.c0, result.c0);

        crs.n_w_1_ctx.invmod(result.c1, b.c1);
        result.c1.mul(a.c1, result.c1);
        crs.n_w_1_ctx.mod(result.c1, result.c1);
    }

    // Multiply input share by a constant.
    void rms_icmult(const crs &crs, const input_share &a, const FMPZ &b, input_share &result) {
        crs.n_w_1_ctx.powm(result.c0, a.c0, b);
        crs.n_w_1_ctx.powm(result.c1, a.c1, b);
    }

    // Multiplies a memory share by an input share
    // unique_id must be different for every RMS multiplication, but must agree between 2 parties for the same multiplication
    // In sequential code, this can be implemented with an incrementing counter
    void rms_mult(const crs &crs, PRF &prf, uint64_t unique_id, const input_share &b, const memory_share &a, memory_share &result, std::optional<int> max_log_y) {
        FMPZ d, tmp, r;

        // std::cout << "a.y_s size: " << fmpz_sizeinbase(a.y_s.get_fmpz(), 2) << std::endl;
        // std::cout << "a.y size: " << fmpz_sizeinbase(a.y.get_fmpz(), 2) << std::endl;

        // d = c0^a.ys * c1^a.y mod N^{w+1}
        b.usage_count++;
        if (b.precomputed) {
            // Use precomputed FBE
            b.c0_fbe.compute(a.y_s, tmp);
            b.c1_fbe.compute(a.y, d);
        } else {
            // Compute c0^a.y_s and c1^a.y
            tmp.powm(b.c0, a.y_s, crs.n_w_1);
            d.powm(b.c1, a.y, crs.n_w_1);
        }
        // tmp.powm(b.c0, a.y_s, crs.n_w_1);
        // d.powm(b.c1, a.y, crs.n_w_1);
        d.mul(d, tmp);
        d.mod(d, crs.n_w_1);

        ddlog(d, result.y_s, crs.n, crs.params.w);
        prf.evaluate_mod_t(unique_id, crs.n_w, crs.params.stat_sec_param, r);
        // Add random shift and mod N^{w}
        result.y_s.add(result.y_s, r);
        if (result.y_s.cmp(crs.n_w) >= 0) {
            result.y_s.sub(result.y_s, crs.n_w);
        }

        int log_y_bound = std::clamp(max_log_y.value_or(crs.params.log2_b), 0, crs.params.log2_b);
        reduce_exponent(result.y_s, result.y_s, crs.params.exponent_bits - crs.params.log2_b + log_y_bound);
        fmpz_fdiv_r_2exp(result.y.get_fmpz(), result.y_s.get_fmpz(), crs.params.stat_sec_param + log_y_bound);
    }
}
