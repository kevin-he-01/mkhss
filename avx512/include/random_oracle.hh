#ifndef RANDOM_ORACLE_HH
#define RANDOM_ORACLE_HH

#include <openssl/sha.h>
#include <openssl/evp.h>

namespace random_oracle {
    constexpr size_t OUTPUT_LENGTH = 32; // SHA-3-256: 32 bytes
    // Magic byte to append to hash inputs to ensure independent randomness between different parts of the program that uses it
    constexpr uint8_t MAGIC_NIKE = 0x00; // Magic number for NIKE
    constexpr uint8_t MAGIC_ANIKE = 0x01; // Magic number for ANIKE
    constexpr uint8_t MAGIC_MKHSS = 0x02; // Magic number for MKHSS

    inline void call(const uint8_t *input, size_t input_length, uint8_t output[OUTPUT_LENGTH]) {
        EVP_MD_CTX *ctx = EVP_MD_CTX_new();
        // SHA-3 is a good model for a random oracle (no length extension attacks unlike SHA-2)
        EVP_DigestInit_ex(ctx, EVP_sha3_256(), NULL);
        EVP_DigestUpdate(ctx, input, input_length);
        EVP_DigestFinal_ex(ctx, output, NULL);
        EVP_MD_CTX_free(ctx);
    }
}

#endif // RANDOM_ORACLE_HH
