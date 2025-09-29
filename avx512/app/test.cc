#include "prf.hh"
#include "mkhss.hh"
#include "util.hh"
#include "flint/fmpz.h"
#include "flint/flint.h"
#include <iostream>
#include <chrono>
#include <functional>
#include <cassert>
#include "ec.hh"
#include "nike.hh"
#include "random_oracle.hh"
#include "fixed_base_exp.hh"

#include "rns_helper.hh"
#include "rns_fbe.hh"
#include "rns_modexp.hh"

/////// DEBUGGING functions below ///////
// This function runs the code cold, would be better to use Google Benchmark or similar
template <typename Func>
double time_function(Func&& func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = (end - start) * 1000.0; // Convert to milliseconds
    std::cout << "Time taken: " << duration.count() << " ms" << std::endl;
    return duration.count();
}

// inline void escape(void *p) {
//     asm volatile("" : : "g"(p) : "memory");
// }
/////// DEBUGGING functions above ///////

void next_safe_prime(mpz_ptr result, const mpz_class& start_base);

void next_safe_base_prime(mpz_class& result, const mpz_class& start) {
    next_safe_prime(Z(result), start);
    result = (result - 1) / 2;
}

void generate_divisor_table();

void test_prf_int() {
    unsigned char key[16];
    get_random_bytes(key, sizeof(key));
    PRF prf(key, 16);
    FMPZ output;
    FMPZ t("100000000000000000000000000", 16);
    std::cout << "     t: 0x" << std::hex << t << std::endl;
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 1024; i += 256) {
            prf.evaluate_mod_t(i, t, 30, output);
            std::cout << "Input: " << std::dec << i << std::endl;
            std::cout << std::endl;
            std::cout << "Output: 0x";
            std::cout << std::hex << output;
            std::cout << std::endl << std::endl;
        }
    }
    std::cout << std::dec;
}

void test_prf() {
    unsigned char key[16];
    get_random_bytes(key, sizeof(key));
    PRF prf(key, 16);
    unsigned char output[32];
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 1024; i += 256) {
            prf.evaluate_raw(i, output, 32);
            std::cout << "Input: " << std::dec << i << std::hex << std::endl;
            std::cout << std::endl;
            std::cout << "Output: ";
            for (int i = 0; i < 32; i++) {
                std::cout << std::hex << (int)output[i];
                if (i != 31) {
                    std::cout << ", ";
                }
            }
            std::cout << std::endl << std::endl;
        }
    }
    std::cout << std::dec;
}

void test_primegen() {
    time_function([&]() {
        generate_divisor_table();
    });

    mpz_class start;
    mpz_ui_pow_ui(Z(start), 2, 1598);
    // The next prime after 2^1598 is 2823 numbers away
    // The next safe prime i 1693069 numbers away
    mpz_class result = -1, gap;
    std::cout << "Start: " << start << std::endl;
    // for (int i = 0; i < 10; i++) { // Warmup
    //     next_safe_base_prime(result, start);
    //     if (i == 0) {
    //         std::cout << "Result: " << result << std::endl;
    //         result -= start;
    //         std::cout << "Gap: " << result << std::endl;
    //     }
    // }
    for (int i = 0; i < 10; i++) {
        time_function([&]() {
            next_safe_base_prime(result, start);
        });
    }

    for (int i = 0; i < 10; i++) {
        mpz_class start;
        gen_random(Z(start), 1535);
        std::cout << "Start: " << start << std::endl;
        mpz_class result = -1;
        double time_us = time_function([&]() {
            next_safe_base_prime(result, start);
        }) * 1000.0; // Convert to microseconds
        result -= start;
        uint64_t gap = mpz_get_ui(Z(result));
        std::cout << "Gap: " << gap << std::endl;
        std::cout << "us per integer: " << time_us / gap << std::endl;
    }
}

void test_dlog() {
    int n_int = 13 * 17;
    int w = 5;
    FMPZ n = n_int;
    FMPZ n_w_1 = n_int;
    for (int i = 0; i < w; i++) {
        n_w_1.mul(n_w_1, n);
    }
    for (int x = 0; x < n_int * n_int; x++) {
        FMPZ x_expected = x;
        FMPZ result;
        FMPZ x_actual;
        pow_base_n_1_mod_n_w_1(result, x_expected, n, 1, n_w_1);
        dlog_z_sp1(n, w, result, x_actual);
        std::cout << "x_expected: " << x_expected << std::endl;
        std::cout << "x_actual: " << x_actual << std::endl;
        assert(x_expected.cmp(x_actual) == 0);
    }
}

void test_import_export_fmpz() {
    // FMPZ a("1234567890123456789012345678901234567890", 10);
    FMPZ a("1fffefd", 16);
    std::cout << "a: " << a << std::endl;
    size_t byte_count = a.get_size_in_bytes();
    std::cout << "byte_count: " << byte_count << std::endl;
    uint8_t *data = (uint8_t *)malloc(byte_count + 1);
    data[byte_count] = 0x42; // Just to check that the buffer is not overrun
    a.export_(data, byte_count);
    assert(data[byte_count] == 0x42);
    for (size_t i = 0; i < byte_count; i++) {
        std::cout << std::hex << (int)data[i];
        if (i != byte_count - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
    FMPZ b;
    b.import(data, byte_count);
    std::cout << "b: " << b << std::endl;
    assert(a.cmp(b) == 0);
    free(data);
}

#include <openssl/ec.h>
#include <openssl/bn.h>
#include <openssl/obj_mac.h>
#include <openssl/err.h>

void handleErrors(void)
{
    ERR_print_errors_fp(stderr);
    abort();
}

// #define NID_CURVE NID_X9_62_prime256v1
// // #define NID_CURVE NID_ED25519

// unsigned char *ecdh_low(size_t *secret_len)
// {
// 	EC_KEY *key, *peerkey;
// 	int field_size;
// 	unsigned char *secret;

// 	/* Create an Elliptic Curve Key object and set it up to use the ANSI X9.62 Prime 256v1 curve */
// 	if(NULL == (key = EC_KEY_new_by_curve_name(NID_CURVE))) handleErrors();

// 	/* Generate the private and public key */
// 	if(1 != EC_KEY_generate_key(key)) handleErrors();

// 	/* Get the peer's public key, and provide the peer with our public key -
// 	 * how this is done will be specific to your circumstances */
// 	// peerkey = get_peerkey_low(key);
//     peerkey = EC_KEY_new_by_curve_name(NID_CURVE);
//     EC_KEY_copy(peerkey, key);

// 	/* Calculate the size of the buffer for the shared secret */
// 	field_size = EC_GROUP_get_degree(EC_KEY_get0_group(key));
// 	*secret_len = (field_size+7)/8;

// 	/* Allocate the memory for the shared secret */
// 	if(NULL == (secret = (unsigned char *) OPENSSL_malloc(*secret_len))) handleErrors();

// 	/* Derive the shared secret */
// 	*secret_len = ECDH_compute_key(secret, *secret_len, EC_KEY_get0_public_key(peerkey),
// 						key, NULL);

// 	/* Clean up */
// 	EC_KEY_free(key);
// 	EC_KEY_free(peerkey);

// 	if(*secret_len <= 0)
// 	{
// 		OPENSSL_free(secret);
// 		return NULL;
// 	}

// 	return secret;
// }

void test_ec() {
    int ret = 1;
    BN_CTX *ctx = NULL;
    EC_GROUP *group = NULL;
    BIGNUM *order = BN_new();
    // BIGNUM *rand_bn = BN_new();
    EC_POINT *pub_key = NULL;
    EC_POINT *shared_key = NULL;
    BIGNUM *priv_key = NULL;
    char *hex = NULL;

    // 1. Create a new BN_CTX (context structure)
    ctx = BN_CTX_new();
    if (ctx == NULL) goto err;

    // 2. Create a new EC_GROUP object with a named curve (e.g., secp256k1)
    group = EC_GROUP_new_by_curve_name(NID_X9_62_prime256v1);
    if (group == NULL) goto err;

    if (EC_GROUP_get_order(group, order, ctx) != 1) {
        goto err;
    }

    // Print the order of the group
    printf("Order of the group: ");
    BN_print_fp(stdout, order);

    // 3. Create a new BIGNUM for the private key
    priv_key = BN_new();
    if (priv_key == NULL) goto err;

    // Set the private key to a specific value (e.g., 1)
    // if (!BN_set_word(priv_key, 2)) goto err;
    // BN_copy(priv_key, order);
    // if (!BN_add_word(priv_key, 42)) goto err;

    if (BN_rand_range(priv_key, order) != 1) goto err;

    // 4. Create a new EC_POINT for the public key
    pub_key = EC_POINT_new(group);
    if (pub_key == NULL) goto err;
    shared_key = EC_POINT_new(group);
    if (shared_key == NULL) goto err;

    // 5. Perform scalar multiplication: pub_key = priv_key * G
    if (!EC_POINT_mul(group, pub_key, priv_key, NULL, NULL, ctx)) goto err;
    // 5. Perform scalar multiplication: shared_key = priv_key * pub_key
    if (!EC_POINT_mul(group, shared_key, NULL, pub_key, priv_key, ctx)) goto err;

    // 6. Print the public key in hexadecimal form
    hex = EC_POINT_point2hex(group, shared_key, POINT_CONVERSION_UNCOMPRESSED, ctx);
    if (hex == NULL) goto err;
    // printf("Shared Key: %s\n", hex);
    OPENSSL_free(hex);

    ret = 0; // Success

err:
    if (ret) {
        handleErrors();
    }

    // Clean up
    EC_POINT_free(pub_key);
    BN_free(priv_key);
    BN_CTX_free(ctx);
    // BN_free(rand_bn);
    BN_free(order);
    EC_GROUP_free(group);
}

void test_nike() {
    // Initialize the CRS
    nike::crs crs;
    nike::setup(crs);

    // Generate keys
    nike::private_key sk_A, sk_B;
    nike::public_key pk_A, pk_B;
    uint8_t shared_key_A[nike::SHARED_KEY_LEN];
    uint8_t shared_key_B[nike::SHARED_KEY_LEN];
    nike::keygen(crs, sk_A, pk_A);
    nike::keygen(crs, sk_B, pk_B);
    nike::keyder(crs, sk_A, pk_B, shared_key_A);
    nike::keyder(crs, sk_B, pk_A, shared_key_B);
    std::cout << "Shared Key A: ";
    for (size_t i = 0; i < nike::SHARED_KEY_LEN; i++) {
        std::cout << std::hex << (int)shared_key_A[i];
        if (i != nike::SHARED_KEY_LEN - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
    std::cout << "Shared Key B: ";
    for (size_t i = 0; i < nike::SHARED_KEY_LEN; i++) {
        std::cout << std::hex << (int)shared_key_B[i];
        if (i != nike::SHARED_KEY_LEN - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
    // Check if the shared keys are equal
    assert(memcmp(shared_key_A, shared_key_B, nike::SHARED_KEY_LEN) == 0);
}

void test_random_oracle() {
    uint8_t input[32];
    uint8_t output[33];
    output[32] = 0x42; // Just to check that the buffer is not overrun
    for (int i = 0; i < 32; i++) {
        input[i] = i + 1;
    }
    random_oracle::call(input, sizeof(input), output);
    std::cout << "Output: ";
    assert(output[32] == 0x42);
    for (int i = 0; i < 32; i++) {
        std::cout << std::hex << (int)output[i];
        if (i != 31) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
}

mpz_class test_n("0x8c70e73ca01499a327a6733c4ebc9a6b23c1724eb2a2ba6023764ba3f7c09206f4936209cc83925c27dc15e1ee542050fe5f75d388c2e938ba76b7b56ac3131d09538dbd157e267b888d17d95053bff688744f76589d5e9cc8e48d7fb82079b2d6127e12ff00bbb89a7df722f950fb54c1f5856bfd5ac348af1847bdaeebd1398a2a34c2de9ff94a7101098714eae611eebed32f37dc1c33de4a0975b16c293973f37caa3368eeae049d0c70a15226198f07aecf051dd6e7ea303379726cf560cb444d0889714928501b83accbef555425218f28120c0758e97a857d18559556339161e185b343680052a3e7ef01e90b0c1f0339f06f01003958f3c8d5e3af0de618a6620bed48a489006a2042f206926aa12bd4b611fcf9bd6cf46c7cb9ad481caa7ac1bf6ff578967fd16bc4a0e47cef10fd3d4dc278922b6f4284959e3149f0325f06a8e531563326df1dcc34941775e75e28468e1c98fd0c12c2bb4dfd4aee3f55ec404aef48e8f5de5273645f26cd5d12726c08ef3f84845b38e54730b1");

void test_fixed_base_exp() {
    PrecomputeIO io("precompute_test_nsquared.txt");
    FMPZ target(test_n);
    target.mul(target, target); // target = n^2
    std::optional<RNSCtx<119>> rns_ctx;
    rns_ctx.emplace(io, target);

    rns_fbe::FixedBaseExp rns_fbe;
    FMPZ base("2", 10); // TODO: larger bases
    int exponent_bits = 900;
    rns_fbe.precompute(base, target, exponent_bits, 1, rns_ctx);
    fixed_base_exp::FixedBaseExp fbe;
    fbe.precompute(base, target, exponent_bits, 1, rns_ctx);
    FMPZ result_expected, result_actual;

    for (int i = 0; i < 1000; i++) {
        // TODO: larger (900-bit exponents)
        // FMPZ exponent("62761dcdde2d39227e4ff53e14852ec15c2f7f0234ec8cd10139237c82f8a257f236626c69de9199a8b0c5b6131636aa5798bb0b936f5640b4eca3f1fe92bb104c48d4a9ef6fa543a4a8abe6f1619e30d2c675f9a7338245de30ff91ff29a4ce3b2198dc", 16);
        FMPZ exponent;
        if (rand() % 10 == 0) {
            exponent.set(0); // Test with zero exponent
        } else {
            gen_random_fmpz(exponent, exponent_bits);
        }
        fbe.compute(exponent, result_expected);
        // rns_fbe.compute(exponent, result_actual);
        // if (result_expected.cmp(result_actual) != 0) {
        //     std::cout << "Mismatch!" << std::endl;
        //     std::cout << "Exponent: 0x" << std::hex << exponent << std::oct << std::endl;
        //     std::cout << "Expected: 0x" << std::hex << result_expected << std::endl;
        //     std::cout << "Actual:   0x" << std::hex << result_actual << std::endl;
        // } 
    }
}

void test_fmpz() {
    FMPZ a, b, c;
    gen_random_fmpz(a, 1000);
    gen_random_fmpz(b, 1000);
    for (int i = 0; i < 10000000; i++) {
        c.mul(a, b);
    }
}

void test_fmpz_mod() {
    FMPZ a, m, c;
    gen_random_fmpz(a, 6144 * 2);
    gen_random_fmpz(m, 6144);
    FMPZModCtx mod_ctx;
    mod_ctx.set_modulus(m);
    for (int i = 0; i < 100000; i++) {
        mod_ctx.mod(c, a);
    }
}

inline void escape(void *p) {
    asm volatile("" : : "g"(p) : "memory");
}

void test_avx_vector() {
    // AVXVector<16> vec1;
    // uint64_t array[16] = {0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7, 0x8,
    //                                 0x9, 0xa, 0xb, 0xc, 0xd, 0xe, 0xf, 0x10};
    // vec1.load(array);
    // vec1.print("array");
    // escape(&vec1);

    constexpr int limbs = 2;
    PrecomputeIO io("/home/ec2-user/mkhss/precompute64.txt");
    FMPZ target("14686114091472174187", 10);

    // constexpr int limbs = 119;
    // PrecomputeIO io("precompute6144.txt");
    // FMPZ target("33534903879232140606314897416440635292778318121545011141787771098890211829116661635200361482790933174257658124765219232207796161958513410038576831349012154285962466325981943093928272305040478823942257124600098538261590473230470630343814668918960377308544749673289050029738457557535118909216817114026816652800744544959306352937535594677428531682581085575329903587005994455326620655513620734869615628995995566431410596276888506543079337189269377879164671187214065648846124686857441046702728286088458348713573043031867255811634380465091679593790892563833949843367899763624556368664648397435728982024787219592465304158314633292702814458308775739986343398068455410911842890459605585047108209859868366846083112089271125878780347798285065647753430443309991527682055167321527826742476367540645077197882800191107364200994428561422791080067454119527748138139576042951722928409119895635143301420004057224629163127771159182975790872133383595595552320710718108043900788739720208199835836842041824100015456316610465017907116324049832283818746240800594413811649422853852466142108052540345081430915291415825734915500675549023806586118714138437155907960471832686016154023256880299891650750579801894468390599172614029282694584371483359228714408523349405308181309144761604801290478162777736009930907654334157137407685960469753880549676039749795347473023354330051239719186915909528715406192606683304578475158159357892837920654789984648888014915040634381071847035775247384311627654435808082053881324411131118153139534240626902431361044229830496277999553029888317569016481642339797126217445356296924528167807806487702804833404830369912680515930001804569866855104159054037369267562960186530271662351372477373680857393000958389906330189389957376397618954610287457223683528293569955158824893658874288052702081827223318373790943134817245035290788780423769037307653604596894347", 10);

    // N^2
    // constexpr int limbs = 119;
    // PrecomputeIO io("precompute_test_nsquared.txt");
    // FMPZ target(test_n);
    // target.mul(target, target); // target = n^2

    RNSCtx<limbs> rns_ctx(io, target);
    // Testing a single mul
    {
        MontInt<limbs> am, bm, cm, dm;
        FMPZ a;
        FMPZ b;
        // gen_random_under_fmpz(a, target);
        // gen_random_under_fmpz(b, target);
        a.set(1);
        b.set(1);
        std::cout << "a: " << a << std::endl;
        std::cout << "b: " << b << std::endl;
        std::cout << "mod target: " << target << std::endl;
        rns_ctx.to_mont_avx(dm, a);
        dm.print();
        am = dm;
        rns_ctx.to_mont(bm, b);
        rns_ctx.mul_mont(cm, am, bm);
        FMPZ c_actual, c_expected;
        c_expected.mul(a, b);
        c_expected.mod(c_expected, target);
        // rns_ctx.from_mont(c_actual, cm);
        rns_ctx.from_mont_avx(c_actual, cm);
        std::cout << "c_expected: " << c_expected << std::endl;
        std::cout << "c_actual: " << c_actual << std::endl;
        assert(c_expected.cmp(c_actual) == 0);
        std::cout << "test passed (single multiplication)" << std::endl;
    }

    // Chain multiplications
    int iters = 10000;
    {
        // AVXVector<limbs> a1, a2, c1, c2;
        MontInt<limbs> am, cm;
        FMPZ a;
        // FMPZ b;
        gen_random_under_fmpz(a, target);
        // gen_random_under_fmpz(b, target);
        // std::cout << "a: " << a << std::endl;
        // std::cout << "b: " << b << std::endl;
        // std::cout << "mod target: " << target << std::endl;
        rns_ctx.to_mont_avx(am, a);
        rns_ctx.to_mont_avx(cm, a);
        // rns_ctx.to_mont(b1, b2, b);
        // am.print();
        // cm.print();
        for (int i = 0; i < iters; i++) {
            rns_ctx.mul_mont(cm, am, cm);
        }
        FMPZ c_actual, c_expected(a);
        for (int i = 0; i < iters; i++) {
            c_expected.mul(c_expected, a);
            c_expected.mod(c_expected, target);
        }
        // rns_ctx.from_mont(c_actual, cm);
        rns_ctx.from_mont_avx(c_actual, cm);
        std::cout << "c_expected: " << c_expected << std::endl;
        std::cout << "c_actual: " << c_actual << std::endl;
        assert(c_expected.cmp(c_actual) == 0);
        std::cout << "test passed (multiple multiplications)" << std::endl;
    }

    // Test to_rns_m1 and to_rns_m2
    // AVXVector<limbs> input;
    // uint64_t input_array[8] = {184012938129ULL};
    // input.load(input_array);
    // AVXVector<limbs> unused; // limbs1
    // AVXVector<limbs> output1, output2;
    // output1 = rns_ctx.multiplier.to_rns_m1.rns_reduce(input, unused, unused, rns_ctx.multiplier.m1);
    // output1.print("output1");
    // output2 = rns_ctx.multiplier.to_rns_m2.rns_reduce(input, unused, unused, rns_ctx.multiplier.m2);
    // output2.print("output2");
}

// void test_exp() {
//     FMPZ a, b, c, c_expected;
//     PrecomputeIO io("precompute_test_nsquared.txt");
//     FMPZ mod(test_n);
//     mod.mul(mod, mod); // target = n^2
//     std::optional<RNSCtx<119>> rns_ctx;
//     rns_ctx.emplace(io, mod);

//     gen_random_under_fmpz(a, mod);
//     for (int iter = 1; iter < 10000; iter++) {
//         gen_random_fmpz(b, 899);
//         rns_modexp::sliding_window(c, a, b, mod, rns_ctx);
//         c_expected.powm(a, b, mod);
//         // assert(c.cmp(c_expected) == 0);
//         if (!(c.cmp(c_expected) == 0)) {
//             std::cout << "Mismatch on iteration " << iter << std::endl;
//             std::cout << "Exponent: 0x" << std::hex << b << std::oct << std::endl;
//             std::cout << "Expected: 0x" << std::hex << c_expected << std::endl;
//             std::cout << "Actual:   0x" << std::hex << c << std::endl;
//         }
//         // std::cout << a << "^" << b << " (mod " << mod << ") = " << c << std::endl;
//         gen_random_fmpz(b, 131);
//         rns_modexp::sliding_window(c, a, b, mod, rns_ctx);
//         c_expected.powm(a, b, mod);
//         // assert(c.cmp(c_expected) == 0);
//         if (!(c.cmp(c_expected) == 0)) {
//             std::cout << "Mismatch on iteration " << iter << std::endl;
//             std::cout << "Exponent: 0x" << std::hex << b << std::oct << std::endl;
//             std::cout << "Expected: 0x" << std::hex << c_expected << std::endl;
//             std::cout << "Actual:   0x" << std::hex << c << std::endl;
//         }
//         // std::cout << a << "^" << b << " (mod " << mod << ") = " << c << std::endl;
//     }
// }

// void test_to_mont() {
//     FMPZ a;
//     MontInt<119> a_mont;
//     PrecomputeIO io("precompute_test_nsquared.txt");
//     FMPZ target(test_n);
//     target.mul(target, target); // target = n^2
//     RNSCtx<119> rns_ctx(io, target);

//     gen_random_under_fmpz(a, target);
//     int num_iters = 10000;
//     double total_time_us = measure_time_seconds([&]() {
//         for (int i = 0; i < num_iters; i++) {
//             rns_ctx.to_mont(a_mont, a);
//         }
//     }) * 1000000.0; // Convert to microseconds
//     std::cout << "Time taken per to_mont call: " << total_time_us / num_iters << " us" << std::endl;
//     std::cout << "Time taken in to_rns: " << us_to_rns_total/ num_iters << " us" << std::endl;

//     // from_mont
//     FMPZ c;
//     double total_from_mont_time_us = measure_time_seconds([&]() {
//         for (int i = 0; i < num_iters; i++) {
//             rns_ctx.from_mont(c, a_mont);
//         }
//     }) * 1000000.0; // Convert to microseconds
//     std::cout << "Time taken per from_mont call: " << total_from_mont_time_us / num_iters << " us" << std::endl;
//     std::cout << "Time taken in from_rns: " << us_from_rns_total / num_iters << " us" << std::endl;
// }

int main() {
    // test_prf();
    // test_prf_int();
    // test_primegen();
    // test_dlog();
    // test_import_export_fmpz();
    // for (int i = 0; i < 1; i++) {
        // test_ec();
    // }
    // size_t secret_len;
    // unsigned char *secret;
    // secret = ecdh_low(&secret_len);
    // for (size_t i = 0; i < secret_len; i++) {
    //     std::cout << std::hex << (int)secret[i];
    //     if (i != secret_len - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // OPENSSL_free(secret);
    // test_random_oracle();
    // test_nike();
    test_fixed_base_exp();
    // test_fmpz();
    // test_fmpz_mod();
    // test_avx_vector();
    // test_exp();
    // test_to_mont();
    return 0;
}
