#ifndef PRF_HH
#define PRF_HH

#include <openssl/evp.h>
#include <cassert>

#include "util.hh"
#include "fmpz_class.hh"

constexpr int MAX_PRF_KEY_LENGTH = 32; // Maximum key length for AES

inline int get_prf_key_length(int lambda) {
    assert(lambda <= 256);
    if (lambda <= 128) {
        return 16;
    } else if (lambda <= 192) {
        return 24;
    } else {
        return 32;
    }
}

class PRF {
private:
    EVP_CIPHER_CTX *ctx;
    mpz_class tmp;
    bool initialized = false;
    int key_length;
public:
    inline PRF() {}

    inline void init(const uint8_t *key, int key_length) {
        if (initialized) {
            EVP_CIPHER_CTX_free(ctx);
        }
        EVP_CIPHER_CTX *ctx = EVP_CIPHER_CTX_new();
        assert(key_length == 16 || key_length == 24 || key_length == 32);
        EVP_EncryptInit_ex(ctx,
                           key_length == 16 ? EVP_aes_128_ecb() :
                           key_length == 24 ? EVP_aes_192_ecb() :
                           key_length == 32 ? EVP_aes_256_ecb() : NULL, NULL, key, NULL);
        EVP_CIPHER_CTX_set_padding(ctx, 0);
        this->ctx = ctx;
        this->key_length = key_length;
        initialized = true;
    }

    inline PRF(const uint8_t *key, int key_length) {
        init(key, key_length);
    }

    // Delete rule of five functions to prevent double free
    PRF(const PRF&) = delete; // Copy constructor
    PRF& operator=(const PRF&) = delete; // Copy assignment operator
    PRF(PRF&&) = delete; // Move constructor
    PRF& operator=(PRF&&) = delete; // Move assignment operator

    inline ~PRF() {
        if (initialized) {
            EVP_CIPHER_CTX_free(ctx);
        }
    }

    inline void evaluate_raw(uint64_t input, uint8_t *output, size_t output_length) {
        uint8_t aes_input[16];
        int len;
        assert(initialized);
        assert(output_length % 16 == 0);
        
        // Set the first 8 bytes of aes_input to input
        for (int i = 0; i < 8; i++) {
            aes_input[i] = (uint8_t) (input >> (i * 8));
        }

        for (size_t i = 0; i < output_length; i += 16) {
            // Set the last 8 bytes of aes_input to i / 16
            for (int j = 0; j < 8; j++) {
                aes_input[8 + j] = (uint8_t) ((i / 16) >> (j * 8));
            }
            // Run the AES pseudorandom permutation
            EVP_EncryptUpdate(ctx, output + i, &len, aes_input, 16);
        }
    }

    inline void derive_subkey(uint64_t id, PRF &subkey) {
        assert(initialized);
        // The number of bytes needed to store the subkey
        uint8_t output[(MAX_PRF_KEY_LENGTH + 15) / 16 * 16];
        // Evaluate the PRF
        evaluate_raw(id, output, (key_length + 15) / 16 * 16);
        // Set the subkey
        subkey.init(output, key_length);
    }

    // Evaluate PRF output as a uniformly random integer mod t
    // tau is the statistical security parameter. Final result is 2^{-tau} indistinguishable from uniform
    inline void evaluate_mod_t(uint64_t input, const FMPZ &t, int tau, FMPZ &out) {
        assert(initialized);
        // The number of bytes needed to store 2^tau * t
        size_t output_bytes = (fmpz_sizeinbase(t.get_fmpz(), 2) + tau) / 8 + 1;
        output_bytes = (output_bytes + 15) / 16 * 16; // Round up to the nearest multiple of 16
        uint8_t *output = new uint8_t[output_bytes];

        // Evaluate the PRF
        evaluate_raw(input, output, output_bytes);
        // Convert the output to a FMPZ
        mpz_import(tmp.get_mpz_t(), output_bytes, 1, sizeof(uint8_t), 1, 0, output);
        out.set(tmp);
        // std::cout << std::dec << "divisor size: " << fmpz_sizeinbase(out.get_fmpz(), 2) << std::endl;
        // std::cout << std::dec << "quotient size: " << fmpz_sizeinbase(t.get_fmpz(), 2) << std::endl;
        out.mod(out, t);

        delete[] output;
    }

    // Evaluate PRF output as a uniformly random integer mod 2^t
    inline void evaluate_mod_2exp(uint64_t input, ulong t, FMPZ &out) {
        assert(initialized);
        // The number of bytes needed to store 2^t
        size_t output_bytes = t / 8 + 1;
        output_bytes = (output_bytes + 15) / 16 * 16; // Round up to the nearest multiple of 16
        uint8_t *output = new uint8_t[output_bytes];

        // Evaluate the PRF
        evaluate_raw(input, output, output_bytes);
        // Convert the output to a FMPZ
        mpz_import(tmp.get_mpz_t(), output_bytes, 1, sizeof(uint8_t), 1, 0, output);
        out.set(tmp);
        // std::cout << std::dec << "divisor size: " << fmpz_sizeinbase(out.get_fmpz(), 2) << std::endl;
        // std::cout << std::dec << "quotient size: " << fmpz_sizeinbase(t.get_fmpz(), 2) << std::endl;
        out.mod_2exp(out, t);

        delete[] output;
    }
};

#endif // PRF_HH
