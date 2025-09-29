#include "util.hh"
#include <gmpxx.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <iostream>
#include <cassert>
#include <vector>

// See GMP: findnext
// TODO: asymptotically faster algorithm (instead of trial division): sieve of eratosthenes. Generate sieve up to the average safe prime gap corresponding to the desired bit length (See GMP: calculate_sievelimit)
// Use FLINT's prime iterator to generate small primes, then sieve an array or bitset (more memory efficient), including only values 5 mod 6
// https://cp-algorithms.com/algebra/sieve-of-eratosthenes.html

// A lot of safe prime generation code inside this file can probably be optimized using FLINT
// Especially https://flintlib.org/doc/ulong_extras.html
// and https://flintlib.org/doc/ulong_extras.html#prime-number-generation-and-counting
// Like prime iterator, precomputed modular inverse, etc.

#define BIT2BYTE(a) (a+7)>>3
// #define DEBUG_SAFE_PRIME
// #define DEBUG_TRIAL_DIVISION_INTEGERS 10000000

// TODO: Is this algorithm memory bound? If so, use uint32_t instead of uint64_t in magic_t to speed up memory access

struct magic_t {
    uint64_t inverse_u64; // divisor^{-1} mod 2^64
    uint64_t threshold; // threshold = 2^64 // divisor + 1
    // uint64_t divisor_shr1; // divisor >> 1 = (divisor - 1) / 2, all divisors are odd (in fact, primes), used to check if (2*number + 1) is divisible by divisor
    uint64_t divisor_shr1_times_inverse_u64; // divisor_shr1 * inverse_u64
};

void number_to_magic(magic_t &magic, mpz_srcptr divisor, const mpz_class& u64) {
    // Can be used to quickly test whether m % number == 0 || (2*m + 1) % number == 0
    mpz_class inverse, threshold;
    mpz_invert(Z(inverse), divisor, Z(u64));
    magic.inverse_u64 = mpz_get_ui(Z(inverse));
    mpz_cdiv_q(Z(threshold), Z(u64), divisor); // divisor is odd, so this is equivalent to floor((2^64) / divisor) + 1
    magic.threshold = mpz_get_ui(Z(threshold));
    // magic.divisor_shr1 = mpz_get_ui(divisor) >> 1; // divisor is odd, so this is equivalent to (divisor - 1) / 2
    magic.divisor_shr1_times_inverse_u64 = (mpz_get_ui(divisor) >> 1) * magic.inverse_u64;
    // std::cout << "# magic.divisor = " << mpz_get_ui(divisor) << std::endl;
    // std::cout << "magic.inverse_u64 = " << magic.inverse_u64 << std::endl;
    // std::cout << "magic.threshold = " << magic.threshold << std::endl;
    // std::cout << "magic.divisor_shr1 = " << magic.divisor_shr1 << std::endl;
}

// checks if either number or (2*number+1) is divisible by divisor inside magic
inline bool divisible_by_magic(uint64_t number, const magic_t& magic) {
    // TODO: check, make sure the subtraction does not underflow
    // return (number * magic.inverse_u64) < magic.threshold || ((number - magic.divisor_shr1) * magic.inverse_u64) < magic.threshold;
    uint64_t product = number * magic.inverse_u64;
    return ((product < magic.threshold) || ((product - magic.divisor_shr1_times_inverse_u64) < magic.threshold));
}

// // Naive implementation
// struct magic_t {
//     uint64_t divisor;
// };

// void number_to_magic(magic_t &magic, mpz_srcptr divisor, const mpz_class& u64) {
//     // Can be used to quickly test whether m % number == 0 || (2*m + 1) % number == 0
//     magic.divisor = mpz_get_ui(divisor);
// }

// // checks if either number or (2*number+1) is divisible by divisor inside magic
// inline bool divisible_by_magic(uint64_t number, const magic_t& magic) {
//     // Really naive implementation, but it works
//     // return (number % magic.divisor == 0) || ((2 * number + 1) % magic.divisor == 0);
//     // Slightly faster implementation (reusing the same modulus reduction result)
//     return (number % magic.divisor == 0) || (number % magic.divisor == (magic.divisor >> 1));
// }

std::vector<std::vector<magic_t>> divisor_table;
uint64_t divisor_product_table[DIVISOR_TABLE_SIZE];

#ifdef DEBUG_TRIAL_DIVISION_INTEGERS
uint64_t num_tests = 0;
#endif
// Fermat primality test function base 2. Used to quickly filter out composites before passing to BPSW (mpz_probab_prime_p)
int fermat_test(mpz_srcptr n) {
    #ifdef DEBUG_TRIAL_DIVISION_INTEGERS
    num_tests++;
    return 0;
    #endif
    // Internal function. Assumes n is odd and greater than 2
    mpz_class n_minus_1, a = 2, remainder;
    mpz_sub_ui(Z(n_minus_1), n, 1);  // n_minus_1 = n - 1

    // Compute a^(n-1) mod n
    mpz_powm(Z(remainder), Z(a), Z(n_minus_1), n);

    // If remainder != 1, n is definitely composite
    if (mpz_cmp_ui(Z(remainder), 1) != 0) {
        return 0;  // composite
    }

    return 1;  // probably prime
}

uint64_t primes_with_product_under_u64(std::vector<magic_t> &primes, const mpz_class& u64, mpz_class& p) {
    // MODIFIES p!!!!
    mpz_class product = 1;

    while (true) {
        if (product * p > u64) {
            break;
        }
        number_to_magic(primes.emplace_back(), Z(p), u64);
        product *= p;

        mpz_nextprime(Z(p), Z(p));
    }

    return mpz_get_ui(Z(product));
}

void generate_divisor_table() {
    mpz_class u64;
    mpz_ui_pow_ui(Z(u64), 2, 64);
    // std::cout << "u64 = " << u64 << std::endl;
    mpz_class current_prime = 5; // By considering only 5 mod 6, we can skip dividing by 2 and 3.
    #ifdef DEBUG_SAFE_PRIME
    std::cout << "Generating divisor table..." << std::endl;
    #endif
    for (size_t i = 0; i < DIVISOR_TABLE_SIZE; i++) {
        auto &divisor_table_row = divisor_table.emplace_back();
        divisor_table_row.reserve(MEDIAN_NUMPRIMES);
        uint64_t product = primes_with_product_under_u64(divisor_table_row, u64, current_prime);
        
        // if (i >= DIVISOR_TABLE_SIZE - 10) {
        //     if (!divisor_table_row.empty()) {
        //         ptrdiff_t occupied = (uint8_t *) divisor_table_row.data() - (uint8_t *) divisor_table[divisor_table.size() - 2].data();
        //         size_t theoretical_size = sizeof(magic_t) * divisor_table_row.size();
        //         // std::cout << "divisor_table_row address: " << (void *) divisor_table_row.data() << std::endl;
        //         std::cout << "divisor_table_row address - previous: " << occupied << std::endl;
        //         std::cout << "divisor_table_row size: " << divisor_table_row.size() << std::endl;
        //         std::cout << "magic_t size: " << sizeof(magic_t) << std::endl;
        //         std::cout << "Theoretical size: " << theoretical_size << std::endl;
        //         std::cout << "Overhead: " << ((double) occupied) / theoretical_size << std::endl;
        //     } else {
        //         std::cout << "divisor_table_row is empty, no address to display." << std::endl;
        //     }
        // }

        // std::cout << "product = " << product << std::endl;
        // std::cout << "divisor_table_row = ";
        // for (size_t j = 0; j < divisor_table_row.size(); j++) {
        //     std::cout << divisor_table_row[j] << " ";
        // }
        // std::cout << std::endl;
        divisor_product_table[i] = product;
    }
    #ifdef DEBUG_SAFE_PRIME
    std::cout << "Smallest prime not trial divided: " << current_prime << std::endl;
    #endif
    // Print tables to ensure correctness
    // for (size_t i = 0; i < DIVISOR_TABLE_SIZE; i++) {
    //     std::cout << "divisor_product_table[" << i << "] = " << divisor_product_table[i] << std::endl;
    //     std::cout << "divisor_table[" << i << "] = ";
    //     for (size_t j = 0; j < divisor_table[i].size(); j++) {
    //         std::cout << divisor_table[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
}

void next_safe_prime(mpz_ptr result, const mpz_class& start_base) {
    // start_base = start / 2 - 1, so start_base needs to be a random number with one fewer bit than desired
    assert(divisor_table.size() == DIVISOR_TABLE_SIZE);
    // Bump base to next number that is 5 mod 6
    mpz_class base = start_base;
    mpz_class probable_base_prime, probable_prime;
    base += 5 - (base % 6);
    // std::cout << "next_safe_prime: start = " << start << std::endl << "base = " << base << std::endl;
    uint64_t base_mods[DIVISOR_TABLE_SIZE];
    for (size_t i = 0; i < DIVISOR_TABLE_SIZE; i++) {
        mpz_class result;
        // Modulo reduction for MPZ is a lot slower than for native integers, so only do it once
        mpz_mod_ui(Z(result), Z(base), divisor_product_table[i]);
        // base_mods[i] = result
        base_mods[i] = mpz_get_ui(Z(result));
        // std::cout << "base_mods[" << i << "] = " << base_mods[i] << std::endl;
    }
    #ifdef DEBUG_TRIAL_DIVISION_INTEGERS
    num_tests = 0;
    auto trial_divide_start = std::chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < DEBUG_TRIAL_DIVISION_INTEGERS; i += 6)
    #else
    for (uint64_t i = 0; true; i += 6)
    #endif
    { // All safe primes are 5 mod 6
        mpz_class divisor_mpz, next_divisor;
        for (size_t j = 0; j < DIVISOR_TABLE_SIZE; j++) {
            uint64_t number = base_mods[j] + i; // hopefully won't overflow uint64_t
            for (magic_t &magic : divisor_table[j]) {
                if (divisible_by_magic(number, magic)) {
                    goto out; // not safe prime
                }
            }
        }
        // Probable prime
        probable_base_prime = base + i;
        probable_prime = 2 * probable_base_prime + 1;
        // std::cout << "probable_base_prime = " << probable_base_prime << std::endl;
        if (fermat_test(Z(probable_base_prime)) && fermat_test(Z(probable_prime))) {
            // Do a stronger test to be sure
            // Composite numbers pass this final test with probability <= 4^{-tau/2} = 2^{-tau}
            if (mpz_probab_prime_p(Z(probable_base_prime), TAU / 2) == 0) {
                #ifdef DEBUG_SAFE_PRIME
                std::cerr << "[!]: probable_base_prime is a Fermat liar" << std::endl;
                #endif
                goto out;
            }
            if (mpz_probab_prime_p(Z(probable_prime), TAU / 2) == 0) {
                #ifdef DEBUG_SAFE_PRIME
                std::cerr << "[!]: probable_prime is a Fermat liar" << std::endl;
                #endif
                goto out;
            }
            mpz_set(result, Z(probable_prime));
            return;
        }
        out:;
    }
    #ifdef DEBUG_TRIAL_DIVISION_INTEGERS
    auto trial_divide_end = std::chrono::high_resolution_clock::now();
    auto trial_divide_duration = std::chrono::duration_cast<std::chrono::microseconds>(trial_divide_end - trial_divide_start);
    auto duration_per_int = (double) trial_divide_duration.count() / DEBUG_TRIAL_DIVISION_INTEGERS;
    std::cout << "Trial division took " << duration_per_int << " us per integer unit" << std::endl;
    std::cout << "Number of tests per 1000 integers: " << (double) num_tests / DEBUG_TRIAL_DIVISION_INTEGERS * 1000 << std::endl;
    mpz_set_si(Z(probable_base_prime), -1);
    #endif
}

// Adapted from https://github.com/camillevuillaume/Paillier-GMP/blob/master/src/tools.c#L69
void get_random_bytes(uint8_t *buf, size_t byte_count) {
    FILE *dev_random, *dev_urandom;
    size_t byte_read;

    dev_random = fopen("/dev/random", "r");
    if (dev_random == NULL) {
        fprintf(stderr, "cannot open random number device!\n");
        exit(1);
    }
    dev_urandom = fopen("/dev/urandom", "r");
    if (dev_urandom == NULL) {
        fprintf(stderr, "cannot open random number device!\n");
        exit(1);
    }

    byte_read = 0;
    while (byte_read < 16 && byte_read < byte_count) {
        byte_read += fread(buf, sizeof(uint8_t), byte_count, dev_random);
    }
    fclose(dev_random);
    while (byte_read < byte_count) {
        byte_read += fread(buf, sizeof(uint8_t), byte_count, dev_urandom);
    }
    fclose(dev_urandom);
}

// Adapted from https://github.com/camillevuillaume/Paillier-GMP/blob/master/src/tools.c#L69
void gen_random(mpz_ptr rnd, mp_bitcnt_t len) {
    int byte_count;
    uint8_t *seed;

    byte_count = BIT2BYTE(len);

    seed = (uint8_t *)malloc(sizeof(uint8_t) * byte_count);

    get_random_bytes(seed, byte_count);

    mpz_import(rnd, byte_count, 1, sizeof(uint8_t), 0, 0, seed);
    mpz_mod_2exp(rnd, rnd, len);
    free(seed);
}

void gen_random_fmpz(FMPZ &rnd, mp_bitcnt_t len) {
    mpz_class rnd_mpz;
    gen_random(Z(rnd_mpz), len);
    rnd.set(rnd_mpz);
}

void gen_random_under(mpz_ptr rnd, mpz_srcptr upper_bound) {
    mpz_class upper_bound_minus_one;
    mpz_sub_ui(Z(upper_bound_minus_one), upper_bound, 1);
    do {
        gen_random(rnd, mpz_sizeinbase(Z(upper_bound_minus_one), 2));
    } while (mpz_cmp(rnd, upper_bound) >= 0);
}

void gen_random_under_fmpz(FMPZ &rnd, const FMPZ &upper_bound) {
    assert(upper_bound.cmp(0) > 0); // upper_bound must be positive
    FMPZ upper_bound_minus_one;
    upper_bound_minus_one.sub(upper_bound, 1);
    do {
        gen_random_fmpz(rnd, fmpz_sizeinbase(upper_bound_minus_one.get_fmpz(), 2));
    } while (rnd.cmp(upper_bound) >= 0);
}

// TODO: parallelize this
void generate_safe_prime(mpz_ptr safe_prime, size_t bits) {
    if (divisor_table.empty()) {
        generate_divisor_table();
    }
    mpz_class cursor;
    gen_random(Z(cursor), bits - 1);
    mpz_setbit(Z(cursor), bits - 2); // Ensure prime is exactly `bits`-bit long
    next_safe_prime(safe_prime, cursor);
    #ifdef DEBUG_SAFE_PRIME
    std::cout << "Safe prime generated: " << safe_prime << std::endl;
    #endif
}

void pow_base_n_1_mod_n_w_1(FMPZ &res, const FMPZ &exp, const FMPZ &n, int w, const FMPZ &n_w_1) {
    if (w == 1) {
        // (N+1)^x === x*N + 1 mod N^2
        res.mul(exp, n);
        res.add(res, 1);
        if (exp.cmp(n) >= 0 || exp.cmp(0) < 0) { // Avoid expensive mod operation if not needed
            res.mod(res, n_w_1);
        }
    } else {
        // This is a more optimized implementation using binomial coefficients
        // This will speed up KeyGen since NIM does (N+1)^s where s is a long secret key
        // See https://link.springer.com/content/pdf/10.1007/s10207-010-0119-9.pdf#page=7.15
        // Section 4.2
        constexpr int n_power_offset = 1;
        std::vector<FMPZ> n_powers; // Stores n^1, n^2, ..., n^{w+1}
        n_powers.reserve(w + 1);
        n_powers.emplace_back(n);
        for (int j = 2; j <= w + 1; j++) {
            FMPZ &next = n_powers.emplace_back();
            next.mul(n_powers[j - 2], n);
        }

        FMPZ c, tmp(exp), tmp2;
        FMPZ j_factorial_inverse_mod_n_w_1(1);
        c.mul(exp, n);
        c.add(c, 1);
        for (int j = 2; j <= w; j++) {
            tmp2.set(j);
            tmp2.invmod(tmp2, n_powers[w + 1 - n_power_offset]); // For large enough N, j! is coprime to n^j
            j_factorial_inverse_mod_n_w_1.mul(j_factorial_inverse_mod_n_w_1, tmp2);
            j_factorial_inverse_mod_n_w_1.mod(j_factorial_inverse_mod_n_w_1, n_powers[w + 1 - n_power_offset]);

            tmp2.sub(exp, j - 1);
            tmp.mul(tmp, tmp2);
            tmp.mod(tmp, n_powers[w - j + 1 - n_power_offset]);

            tmp2.set(tmp);
            tmp2.mul(tmp2, j_factorial_inverse_mod_n_w_1);
            tmp2.mod(tmp2, n_powers[w + 1 - n_power_offset]);
            tmp2.mul(tmp2, n_powers[j - n_power_offset]);
            c.add(c, tmp2);
            c.mod(c, n_powers[w + 1 - n_power_offset]);
        }
        res.set(c);
    }
    // Naive implementation:
    // FMPZ base;
    // base.add(n, 1);
    // res.powm(base, exp, n_w_1);
}

// https://brics.dk/RS/00/45/BRICS-RS-00-45.pdf
// See page 4-5 of the paper linked above
void dlog_z_sp1(const FMPZ &n, int w, const FMPZ &a, FMPZ &i) {
    assert(&i != &a); // Does not support in-place computation
    assert(&i != &n); // Does not support in-place computation

    FMPZ a_minus_1;
    FMPZ t1, t2;
    FMPZ k_factorial_inverse_mod_nj;
    FMPZ tmp1, tmp2, tmp3;

    constexpr int n_power_offset = 1;
    std::vector<FMPZ> n_powers; // Stores n^1, n^2, ..., n^{w+1}
    n_powers.reserve(w + 1);
    n_powers.emplace_back(n);
    for (int j = 2; j <= w + 1; j++) {
        FMPZ &next = n_powers.emplace_back();
        next.mul(n_powers[j - 2], n);
    }

    i.set(0);
    a_minus_1.sub(a, 1);

    for (int j = 1; j <= w; j++) {
        tmp3.mod(a_minus_1, n_powers[j + 1 - n_power_offset]);
        fmpz_fdiv_q(t1.get_fmpz(), tmp3.get_fmpz(), n.get_fmpz());
        
        t2.set(i);

        k_factorial_inverse_mod_nj.set(1);
        for (int k = 2; k <= j; k++) {
            tmp2.set(k);
            tmp2.invmod(tmp2, n_powers[j - n_power_offset]); // For large enough N, k! is coprime to n^j
            k_factorial_inverse_mod_nj.mul(k_factorial_inverse_mod_nj, tmp2);
            k_factorial_inverse_mod_nj.mod(k_factorial_inverse_mod_nj, n_powers[j - n_power_offset]);

            i.sub(i, 1);
            t2.mul(t2, i);
            t2.mod(t2, n_powers[j - n_power_offset]);
            
            tmp1.mul(t2, n_powers[k - 1 - n_power_offset]);
            tmp1.mod(tmp1, n_powers[j - n_power_offset]);
            tmp1.mul(tmp1, k_factorial_inverse_mod_nj);

            t1.sub(t1, tmp1);
            t1.mod(t1, n_powers[j - n_power_offset]);
        }
        i.set(t1);
    }
}

void ddlog(const FMPZ &g_prime, FMPZ &z, const FMPZ &n, int w) {
    // z = DDLog(g')
    // When w = 1
    // g' = h + h'N
    // z = h' h^{-1} mod N
    if (w == 1) {
        FMPZ h_prime;
        fmpz_fdiv_qr(h_prime.get_fmpz(), z.get_fmpz(), g_prime.get_fmpz(), n.get_fmpz());
        fmpz_invmod(z.get_fmpz(), z.get_fmpz(), n.get_fmpz());
        z.mul(z, h_prime);
        z.mod(z, n);
    } else {
        // Using the NIDLS framework, Section 4.1 instantiation of DDLog:
        // https://eprint.iacr.org/2022/363.pdf#page=14.64
        FMPZ h;
        FMPZ n_w_1;
        fmpz_pow_ui(n_w_1.get_fmpz(), n.get_fmpz(), w + 1);
        h.mod(g_prime, n); // h = g' mod N, the coset label
        // std::cout << "h: " << h << std::endl;
        h.invmod(h, n_w_1);
        h.mul(h, g_prime);
        h.mod(h, n_w_1);
        dlog_z_sp1(n, w, h, z);
    }
}
