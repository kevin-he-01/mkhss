#ifndef UTIL_HH
#define UTIL_HH

#include <gmpxx.h>
#include "fmpz_class.hh"
#include <cassert>
#include <chrono>

#define Z(mpz_class_object) (mpz_class_object.get_mpz_t())

// 20000 seems to be a sweet spot that balances division and primality testing cost, assuming 1536 bit primes
constexpr size_t DIVISOR_TABLE_SIZE = 20000;
// How much memory to reserve for the vector holding a list of primes whose product is less than 2^64
// A smaller number allow the primes to be packed denser in memory
// 3 seems to be optimal for DIVISOR_TABLE_SIZE = 20000
constexpr size_t MEDIAN_NUMPRIMES = 3;
// Statistical security parameter (for correctness)
constexpr int TAU = 128;

void get_random_bytes(uint8_t *buf, size_t byte_count);
void gen_random(mpz_ptr rnd, mp_bitcnt_t len);
void gen_random_fmpz(FMPZ &rnd, mp_bitcnt_t len);
void gen_random_under(mpz_ptr rnd, mpz_srcptr upper_bound);
void gen_random_under_fmpz(FMPZ &rnd, const FMPZ &upper_bound);
void generate_safe_prime(mpz_ptr safe_prime, size_t bits);
void pow_base_n_1_mod_n_w_1(FMPZ &res, const FMPZ &exp, const FMPZ &n, int w, const FMPZ &n_w_1);
void dlog_z_sp1(const FMPZ &n, int w, const FMPZ &a, FMPZ &i);
void ddlog(const FMPZ &g_prime, FMPZ &z, const FMPZ &n, int w);

inline int bitcnt_n(int n) {
    assert(n > 0);
    // Obtain the number of bits needed to represent n
    int bitcnt = 0;
    while (n) {
        n >>= 1;
        bitcnt++;
    }
    return bitcnt;
}

// Machine-endianness-agnostic serialization of integers
template <typename T>
inline void serialize_integer(uint8_t *buf, T integer) {
    static_assert(std::is_integral<T>::value, "T must be an integral type");
    for (size_t i = 0; i < sizeof(T); i++) {
        buf[i] = static_cast<uint8_t>((integer >> (i * 8)) & 0xFF);
    }
}

template <typename Func>
inline double measure_time_seconds(Func&& func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = (end - start);
    return duration.count();
}

#endif // UTIL_HH
