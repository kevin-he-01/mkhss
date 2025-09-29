#ifndef ANIKE_HH
#define ANIKE_HH

#include "random_oracle.hh"

#ifdef BASELINE_MKHSS
#include "baseline_mkhss.hh"
#define mkhss baseline_mkhss
#else
#include "mkhss.hh"
#endif

namespace anike {
    constexpr size_t KEY_LENGTH = 16; // Length of a single key in bytes
    
    inline void derive_key(uint8_t nike_shared_key[nike::SHARED_KEY_LEN], const mkhss::memory_share &z, uint8_t key[KEY_LENGTH]) {
        // uint8_t shared_key[KEY_LENGTH]; // Instead of concatenating the identities id_A || id_B, use the NIKE shared key f to simplify implementation
        // PRF &shared = ctx.get_prf_shared();
        // shared.evaluate_raw(0, shared_key, KEY_LENGTH);

        FMPZ k_i = z.y_s;
        // std::cout << "K_i: " << k_i << std::endl;
        size_t k_i_size = z.y_s.get_size_in_bytes();
        uint8_t *buffer = new uint8_t[k_i_size + nike::SHARED_KEY_LEN + 1];
        k_i.export_(buffer, k_i_size);
        memcpy(buffer + k_i_size, nike_shared_key, nike::SHARED_KEY_LEN);
        buffer[k_i_size + nike::SHARED_KEY_LEN] = random_oracle::MAGIC_ANIKE; // Magic number for ANIKE
        // Print buffer
        // std::cout << "Buffer: ";
        // for (size_t i = 0; i < k_i_size + KEY_LENGTH + 1; i++) {
        //     std::cout << std::hex << std::setfill('0') << std::setw(2) << (int) buffer[i];
        // }
        // std::cout << std::dec << std::endl;
        uint8_t hash[random_oracle::OUTPUT_LENGTH];
        random_oracle::call(buffer, k_i_size + nike::SHARED_KEY_LEN + 1, hash); // K_i = H(k_i || k_{nike})
        static_assert(KEY_LENGTH <= random_oracle::OUTPUT_LENGTH);
        memcpy(key, hash, KEY_LENGTH);

        delete[] buffer;
    }
}

#undef mkhss

#endif // ANIKE_HH
