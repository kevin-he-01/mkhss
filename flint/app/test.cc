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
#include "modexp.hh"

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

void test_fixed_base_exp() {
    FixedBaseExp fbe, fbe2;
    // fbe2.precompute(2, 10000, 8, 1);
    // fbe.downcast(fbe2, 1000); // Downcast to a smaller modulus
    fbe.precompute(2, 1000, 8, 4); // Precompute powers of 2 mod 1000
    FMPZ result;
    // for (int i = -10; i < 10; i++) {
    //     FMPZ exponent;
    //     exponent.set(i);
    //     fbe.compute(exponent, result);
    //     std::cout << "2^" << i << " mod 1001 = " << result << std::endl;
    // }
    for (int i = 0; i < 100; i++) {
        FMPZ exponent;
        exponent.set(i);
        fbe.compute(exponent, result);
        std::cout << "2^" << i << " mod 1000 = " << result << std::endl;
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

void test_modexp() {
    FMPZ base, exponent, modulus, result;
    gen_random_fmpz(base, 2048);
    gen_random_fmpz(exponent, 2048);
    gen_random_fmpz(modulus, 2048);
    // FMPZ base("2", 10), exponent("2", 10), modulus("1000", 10), result;
    // modulus.add(modulus, 1); // Make sure modulus is not zero
    std::cout << "Base: " << base << std::endl;
    std::cout << "Exponent: " << exponent << std::endl;
    std::cout << "Modulus: " << modulus << std::endl;

    modexp::sliding_window(result, base, exponent, modulus);
    std::cout << "Result: " << result << std::endl;

    FMPZ expected;
    expected.powm(base, exponent, modulus);
    std::cout << "Result (expected): " << expected << std::endl;
    assert(result.cmp(expected) == 0);
}

void test_straus() {
    FMPZ base1, base2, exponent1, exponent2, modulus(
        "5103107934809809939552244207115506393310956973793356350273061938637513624233398735886878273640175010327914486504003513299938314617211182700029517511845263419434737062340722923261868297648541688741941971318797127173031366926596836787437926876276970899985914650252830605497052526704505178416480357751254601345373648084449815066346641899706590033987943151860925374744679452680044717753967707317322707976876274854788495657109053368187367373277148444153011564450396825553314693768811337112289124302339669838098050348458397195220135848861092623082643504914608564584988068331278698691040405941939394670877291579544599169117",
        10
    ), result;
    gen_random_fmpz(base1, 2048);
    gen_random_fmpz(base2, 2048);
    gen_random_fmpz(exponent1, 2048);
    fmpz_neg(exponent1.get_fmpz(), exponent1.get_fmpz()); // Make exponent1 negative
    gen_random_fmpz(exponent2, 100);
    // fmpz_neg(exponent2.get_fmpz(), exponent2.get_fmpz()); // Make exponent2 negative
    // generate_safe_prime(modulus, 2048);
    // FMPZ base("2", 10), exponent("2", 10), modulus("1000", 10), result;
    // modulus.add(modulus, 1); // Make sure modulus is not zero
    std::cout << "Base 1: " << base1 << std::endl;
    std::cout << "Base 2: " << base2 << std::endl;
    std::cout << "Exponent 1: " << exponent1 << std::endl;
    std::cout << "Exponent 2: " << exponent2 << std::endl;
    std::cout << "Modulus: " << modulus << std::endl;

    std::vector<FMPZ> bases = {base1, base2};
    std::vector<FMPZ> exponents = {exponent1, exponent2};
    modexp::straus(result, bases, exponents, modulus);
    // modexp::sliding_window(result, base, exponent, modulus);
    std::cout << "Result: " << result << std::endl;

    FMPZ expected, tmp;
    expected.powm(base1, exponent1, modulus);
    tmp.powm(base2, exponent2, modulus);
    expected.mul(expected, tmp);
    expected.mod(expected, modulus);
    std::cout << "Result (expected): " << expected << std::endl;
    assert(result.cmp(expected) == 0);
}

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
    // test_fixed_base_exp();
    // test_fmpz();
    // test_fmpz_mod();
    // test_modexp();
    test_straus();
    flint_cleanup_master();
    return 0;
}
