#include "baseline_mkhss.hh"
#include "mkhss.hh"
#include "util.hh"
#include "flint/fmpz.h"
#include "flint/flint.h"
#include <iostream>
#include <chrono>
#include <functional>
#include <cassert>

#ifdef BASELINE_MKHSS
#define mkhss baseline_mkhss
#endif

// p = 1846676733488631055471255842241675918524634969943310544305024546510335111754942467855175334391304886738089012161892806684010679369213995855689802335628273984374651322457100666803755783382146601342956837891522725255623204154321951439681667011710360653069057516378690718722316644878184359035522184553662331270365968446822008203609426500960578829207809935403657564047100423746644192595412017524710813488670533965603293464156021031202360140589859051220120924452326967
// q = 1725876987217986047624500396971859681137265214517892022764955133522145650537095430314045088828089319223467973204866689257287084472431639057265828454653668217896505768149095595788843000748135974986150770180016555436749641844269146519701313369089352948184435258565377439750312707561252924698656404987527139949369182619563623089969086836486319123231640554043654596504290380661637096995077501815992908334347140447815056438846057257331379512924564768728797003234542679
mpz_class test_n("0x8c70e73ca01499a327a6733c4ebc9a6b23c1724eb2a2ba6023764ba3f7c09206f4936209cc83925c27dc15e1ee542050fe5f75d388c2e938ba76b7b56ac3131d09538dbd157e267b888d17d95053bff688744f76589d5e9cc8e48d7fb82079b2d6127e12ff00bbb89a7df722f950fb54c1f5856bfd5ac348af1847bdaeebd1398a2a34c2de9ff94a7101098714eae611eebed32f37dc1c33de4a0975b16c293973f37caa3368eeae049d0c70a15226198f07aecf051dd6e7ea303379726cf560cb444d0889714928501b83accbef555425218f28120c0758e97a857d18559556339161e185b343680052a3e7ef01e90b0c1f0339f06f01003958f3c8d5e3af0de618a6620bed48a489006a2042f206926aa12bd4b611fcf9bd6cf46c7cb9ad481caa7ac1bf6ff578967fd16bc4a0e47cef10fd3d4dc278922b6f4284959e3149f0325f06a8e531563326df1dcc34941775e75e28468e1c98fd0c12c2bb4dfd4aee3f55ec404aef48e8f5de5273645f26cd5d12726c08ef3f84845b38e54730b1");

// Parameters
constexpr bool SHORT_EXPONENT_ASSUMPTION = true;
constexpr int SEC_PARAM = 128; // Lambda. Computational security parameter (against offline attacks)
constexpr int STAT_SEC_PARAM = TAU; // Tau. Statistical security parameter (for correctness)

constexpr int OURS_W = 1;
constexpr int BASELINE_W = 3;
constexpr int LOG2_B = 14; // log_2(B), where B is the bound on the absolute value of all integers in memory shares
constexpr int OURS_N_WIDTH = 3072; // Needed for 128 bits of security against factoring
constexpr int BASELINE_N_WIDTH = 3072; // Needed for 128 bits of security against factoring

// constexpr int OURS_W = 1;
// constexpr int BASELINE_W = 3;
// constexpr int LOG2_B = 1024; // log_2(B), where B is the bound on the absolute value of all integers in memory shares
// constexpr int OURS_N_WIDTH = 3968;
// constexpr int BASELINE_N_WIDTH = 3072; // Needed for 128 bits of security against factoring

// constexpr int OURS_W = 2;
// constexpr int BASELINE_W = 3;
// constexpr int LOG2_B = 2048; // log_2(B), where B is the bound on the absolute value of all integers in memory shares
// constexpr int OURS_N_WIDTH = 3520;
// constexpr int BASELINE_N_WIDTH = 3072; // Needed for 128 bits of security against factoring

// constexpr int OURS_W = 4;
// constexpr int BASELINE_W = 4;
// constexpr int LOG2_B = 4096; // log_2(B), where B is the bound on the absolute value of all integers in memory shares
// constexpr int OURS_N_WIDTH = 3296;
// constexpr int BASELINE_N_WIDTH = 4224;

constexpr int W = mkhss::IS_BASELINE ? BASELINE_W : OURS_W;
constexpr int N_WIDTH = mkhss::IS_BASELINE ? BASELINE_N_WIDTH : OURS_N_WIDTH; // Use the same N_WIDTH for both baseline and optimized MKHSS

void eval_rms(const mkhss::crs &crs, mkhss::rms_context &ctx, const mkhss::input_share &x_synced, const mkhss::input_share &y_synced, FMPZ &output) {
    // uint64_t inst_id = 0;
    mkhss::memory_share x_mem, x_1_mem, x_minus_y_mem;
    mkhss::input_share two_x, x_minus_y;

    mkhss::rms_iadd(crs, x_synced, x_synced, two_x);
    // std::cout << "2x.c0 bits: " << fmpz_bits(two_x.c0.get_fmpz()) << std::endl;
    // std::cout << "2x.c1 bits: " << fmpz_bits(two_x.c1.get_fmpz()) << std::endl;
    mkhss::rms_isub(crs, x_synced, y_synced, x_minus_y);
    // std::cout << "{x-y}.c0 bits: " << fmpz_bits(x_minus_y.c0.get_fmpz()) << std::endl;
    // std::cout << "{x-y}.c1 bits: " << fmpz_bits(x_minus_y.c1.get_fmpz()) << std::endl;
    mkhss::rms_icmult(crs, x_synced, 20, two_x);

    ctx.convert_serial(crs, two_x, x_mem);
    ctx.mult_serial(crs, x_minus_y, x_mem, x_minus_y_mem);
    // mkhss::rms_mcmult(x_mem, 10, x_mem);
    // ctx.mult_serial(zero, x_mem, x_mem);
    // mkhss::rms_madd(x_mem, one, x_1_mem);
    mkhss::rms_msub(x_minus_y_mem, ctx.get_one(), x_1_mem);
    // mkhss::rms_output(x_1_mem, output);
    // mkhss::rms_output(x_1_mem, output);
    ctx.output_serial(crs, x_1_mem, output);
}

void eval_rms_simple(const mkhss::crs &crs, mkhss::rms_context &ctx, const mkhss::input_share &x_synced, FMPZ &output) {
    uint64_t inst_id = 0;
    mkhss::memory_share x_mem;

    ctx.convert_serial(crs, x_synced, x_mem);

    ctx.output_serial(crs, x_mem, output);
}

int main() {
    mkhss::crs crs;
    mkhss::parameters params;

    if constexpr (mkhss::IS_BASELINE) {
        std::cout << "Testing correctness of baseline MKHSS" << std::endl;
    } else {
        std::cout << "Testing correctness of this work (optimized MKHSS)" << std::endl;
    }
    
    mkhss::create_parameters(params, SHORT_EXPONENT_ASSUMPTION, STAT_SEC_PARAM, SEC_PARAM, N_WIDTH, LOG2_B, W);
    mkhss::print_parameters(params);

    std::cout << "Generating CRS..." << std::endl;
    // mkhss::setup(crs, params, test_n.get_mpz_t()); // Use pregenerated n to speed up edit-debug cycle
    mkhss::setup(crs, params); // Use pregenerated n to speed up edit-debug cycle
    std::cout << "CRS generated" << std::endl;

    // return 0;

    // crs = test_crs;
    mkhss::public_key pk_0, pk_1;
    mkhss::private_key sk_0, sk_1;
    mkhss::keygen(crs, pk_0, sk_0);
    mkhss::keygen(crs, pk_1, sk_1);
    std::cout << "Public key size (Alice): " << mkhss::get_public_key_size(crs, pk_0) << std::endl;
    std::cout << "Public key size (Bob): " << mkhss::get_public_key_size(crs, pk_1) << std::endl;
    mkhss::private_share x_0;
    mkhss::public_share x_1;
    mkhss::public_share y_0;
    mkhss::private_share y_1;
    // FMPZ x = 2;
    // fmpz_ui_pow_ui(x.get_fmpz(), 2, 900);
    FMPZ x = 42, y = 30;
    mkhss::share(crs, pk_0, x, x_0, x_1);
    mkhss::input_share x_synced_0, x_synced_1;
    mkhss::sync_share_self(crs, sk_0, pk_1, x_0, x_synced_0);
    mkhss::sync_share_other(crs, sk_1, pk_0, x_1, x_synced_1);

    mkhss::share(crs, pk_1, y, y_1, y_0);
    std::cout << "Share size (Alice): " << mkhss::get_public_share_size(crs, x_1) << std::endl;
    std::cout << "Share size (Bob): " << mkhss::get_public_share_size(crs, y_0) << std::endl;
    mkhss::input_share y_synced_0, y_synced_1;
    mkhss::sync_share_other(crs, sk_0, pk_1, y_0, y_synced_1);
    mkhss::sync_share_self(crs, sk_1, pk_0, y_1, y_synced_0);

    assert(x_synced_0.c0.cmp(x_synced_1.c0) == 0);
    assert(x_synced_0.c1.cmp(x_synced_1.c1) == 0);
    assert(y_synced_0.c0.cmp(y_synced_1.c0) == 0);
    assert(y_synced_0.c1.cmp(y_synced_1.c1) == 0);

    #ifdef BASELINE_MKHSS
    assert(x_synced_0.c0_prime.cmp(x_synced_1.c0_prime) == 0);
    assert(x_synced_0.c1_prime.cmp(x_synced_1.c1_prime) == 0);
    assert(y_synced_0.c0_prime.cmp(y_synced_1.c0_prime) == 0);
    assert(y_synced_0.c1_prime.cmp(y_synced_1.c1_prime) == 0);
    #endif
    // std::cout << "Sync share self: " << x_synced_0.c0 << ", " << x_synced_0.c1 << std::endl;
    // std::cout << "Sync share other: " << x_synced_1.c0 << ", " << x_synced_1.c1 << std::endl;
    std::cout << "Sync share OK" << std::endl;

    std::vector<std::pair<const mkhss::private_key &, const mkhss::public_key &>> keys = {
        {sk_0, pk_0},
        {sk_1, pk_1}
    };

    FMPZ z_true;
    for (int i = 0; i <= 1; i++) {
        // Test reusability by making each party alternate between being Alice and Bob
        const mkhss::private_key &sk_alice = keys[i].first;
        const mkhss::public_key &pk_alice = keys[i].second;
        const mkhss::private_key &sk_bob = keys[1-i].first;
        const mkhss::public_key &pk_bob = keys[1-i].second;
        FMPZ s, s_expected, one;
        mkhss::rms_context ctx_a, ctx_b;
        mkhss::rms_bootstrap(crs, mkhss::party_id::ALICE, sk_alice, pk_bob, ctx_a);
        mkhss::rms_bootstrap(crs, mkhss::party_id::BOB, sk_bob, pk_alice, ctx_b);
        const mkhss::memory_share &m_a = ctx_a.get_one(), &m_b = ctx_b.get_one();

        s.sub(m_a.y_s, m_b.y_s);
        one.sub(m_a.y, m_b.y);
        s_expected.mul(sk_0.s, sk_1.s);
        assert(s_expected.cmp(crs.n_w) < 0);
        // std::cout << "s         : " << s << std::endl;
        // std::cout << "s_expected: " << s_expected << std::endl;
        assert(s.cmp(s_expected) == 0);
        assert(one.cmp(1) == 0);
        std::cout << "NIM decode OK" << std::endl;
    
        FMPZ z_A, z_B, z;
        eval_rms(crs, ctx_a, x_synced_0, y_synced_0, z_A);
        eval_rms(crs, ctx_b, x_synced_0, y_synced_0, z_B);
        // eval_rms_simple(crs, ctx_a, x_synced_0, z_A);
        // eval_rms_simple(crs, ctx_b, x_synced_0, z_B);
        z.sub(z_A, z_B);
        if (i == 0) {
            z_true.set(z);
            std::cout << "z: " << z << std::endl;
        } else {
            assert(z.cmp(z_true) == 0);
        }
    }

    return 0;
}
