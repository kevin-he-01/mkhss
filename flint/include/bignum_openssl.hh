#ifndef BIGNUM_OPENSSL_HH
#define BIGNUM_OPENSSL_HH

#include <openssl/bn.h>

// A wrapper around OpenSSL BIGNUM to simplify memory management
class BIGNUMClass {
private:
    BIGNUM *bn;
public:
    inline BIGNUMClass() {
        bn = BN_new();
    }

    // Delete copy constructor and assignment operator
    BIGNUMClass(const BIGNUMClass&) = delete;
    BIGNUMClass& operator=(const BIGNUMClass&) = delete;
    // Delete move constructor and assignment operator
    BIGNUMClass(BIGNUMClass&&) = delete;
    BIGNUMClass& operator=(BIGNUMClass&&) = delete;

    inline ~BIGNUMClass() {
        BN_free(bn);
    }

    inline BIGNUM *get_bn() {
        return bn;
    }

    inline const BIGNUM *get_bn() const {
        return bn;
    }
};

#endif // BIGNUM_OPENSSL_HH
