#ifndef FMPZ_CLASS_HH
#define FMPZ_CLASS_HH

#include <gmpxx.h>
#include <flint/fmpz.h>
#include <iostream>
#include <cassert>
#include <cstdint>

// Wrapper around FLINT fmpz_t in a class to simplify memory management
// and provide a more C++-like interface.
class FMPZ {
private:
    fmpz_t value;
    bool is_modulus = false; // RNS Montgomery modulus indicator
public:
    inline FMPZ() {
        fmpz_init(value);
    }

    inline FMPZ(slong val) {
        // std::cout << "Creating FMPZ from slong: " << val << std::endl;
        fmpz_init_set_si(value, val);
    }

    inline FMPZ(const char *str, int b) {
        fmpz_init(value);
        fmpz_set_str(value, str, b);
    }

    inline FMPZ(const mpz_class &mpz_value) {
        fmpz_init(value);
        fmpz_set_mpz(value, mpz_value.get_mpz_t());
    }

    inline FMPZ(const FMPZ &other) {
        fmpz_init_set(value, other.value);
        is_modulus = other.is_modulus;
    }

    // Can technically implement copy assignment. Just use set() for consistent style
    // inline FMPZ &operator=(const FMPZ &other) {
    //     std::cout << "Copying FMPZ from another FMPZ" << std::endl;
    //     fmpz_set(value, other.value);
    //     return *this;
    // }

    FMPZ &operator=(const FMPZ &other) = delete;

    // Don't know how to implement the move constructor properly
    FMPZ(FMPZ &&other) = delete;

    FMPZ &operator=(FMPZ &&other) = delete;

    inline ~FMPZ() {
        fmpz_clear(value);
    }

    inline void set(slong val) {
        fmpz_set_si(value, val);
    }

    inline void set(const mpz_class &other) {
        fmpz_set_mpz(value, other.get_mpz_t());
    }

    inline void set_is_modulus(bool is_mod) {
        is_modulus = is_mod;
    }

    inline size_t get_size_in_bytes() const {
        return (fmpz_sizeinbase(value, 2) + 7) / 8; // Convert bits to bytes
    }

    inline void import(const uint8_t *data, size_t length) {
        mpz_class mpz_value;
        mpz_import(mpz_value.get_mpz_t(), length, 1, sizeof(uint8_t), 1, 0, data);
        fmpz_set_mpz(value, mpz_value.get_mpz_t());
    }

    inline void export_(uint8_t *data, size_t length) const {
        assert(length >= get_size_in_bytes() && "Buffer too small for export");
        mpz_class mpz_value;
        fmpz_get_mpz(mpz_value.get_mpz_t(), value);
        mpz_export(data, nullptr, 1, sizeof(uint8_t), 1, 0, mpz_value.get_mpz_t());
    }

    inline void set(const FMPZ &other) {
        // std::cout << "Setting FMPZ from another FMPZ" << std::endl;
        fmpz_set(value, other.value);
    }

    inline void set_mpz(mpz_srcptr mpz_value) {
        fmpz_set_mpz(value, mpz_value);
    }

    inline void add(const FMPZ &a, slong b) {
        fmpz_add_si(value, a.value, b);
    }

    inline void add(const FMPZ &a, const FMPZ &b) {
        fmpz_add(value, a.value, b.value);
    }

    inline void sub(const FMPZ &a, slong b) {
        fmpz_sub_si(value, a.value, b);
    }

    inline void sub(const FMPZ &a, const FMPZ &b) {
        fmpz_sub(value, a.value, b.value);
    }

    inline void mul(const FMPZ &a, slong b) {
        fmpz_mul_si(value, a.value, b);
    }

    inline void mul(const FMPZ &a, const FMPZ &b) {
        fmpz_mul(value, a.value, b.value);
    }

    inline void mul_2exp(const FMPZ &a, ulong exp) {
        fmpz_mul_2exp(value, a.value, exp);
    }

    inline void mod(const FMPZ &a, const FMPZ &b) {
        fmpz_mod(value, a.value, b.value);
    }

    inline void mod(const FMPZ &a, ulong b) {
        fmpz_mod_ui(value, a.value, b);
    }

    inline void mod_2exp(const FMPZ &a, ulong exp) {
        fmpz_fdiv_r_2exp(value, a.value, exp);
    }

    inline void invmod(const FMPZ &a, const FMPZ &b) {
        fmpz_invmod(value, a.value, b.value);
    }

    inline int cmp(const FMPZ &b) const {
        return fmpz_cmp(value, b.value);
    }

    inline int cmp(slong b) const {
        return fmpz_cmp_si(value, b);
    }

    inline void pow(const FMPZ &base, ulong exp) {
        fmpz_pow_ui(value, base.value, exp);
    }

    // inline void powm(const FMPZ &base, const FMPZ &exp, const FMPZ &mod) {
    //     fmpz_powm(value, base.value, exp.value, mod.value);
    // }
    void powm(const FMPZ &base, const FMPZ &exp, const FMPZ &mod);

    inline ulong get_ui() const {
        return fmpz_get_ui(value);
    }

    inline slong get_si() const {
        return fmpz_get_si(value);
    }

    inline void to_mpz(mpz_ptr mpz) const {
        fmpz_get_mpz(mpz, value);
    }

    inline void to_mpz(mpz_class &mpz_value) const {
        fmpz_get_mpz(mpz_value.get_mpz_t(), value);
    }

    inline const fmpz *get_fmpz() const {
        return value;
    }

    inline fmpz *get_fmpz() {
        return value;
    }
};

inline std::ostream & operator<<
(std::ostream &o, const FMPZ &expr)
{
    mpz_class mpz;
    expr.to_mpz(mpz);
    return o << mpz; // gmpxx.h provides operator<< for mpz_class
}

#endif
