#ifndef EC_HH
#define EC_HH

#include <openssl/ec.h>
#include <openssl/obj_mac.h>
#include <openssl/err.h>
#include "bignum_openssl.hh"

// A wrapper around OpenSSL EC_GROUP to simplify memory management
class ECGroup {
private:
    mutable BN_CTX *ctx;
    EC_GROUP *group;
public:
    inline ECGroup() {
        ctx = BN_CTX_new();
        group = EC_GROUP_new_by_curve_name(NID_X9_62_prime256v1);
    }

    // Delete copy constructor and assignment operator
    ECGroup(const ECGroup&) = delete;
    ECGroup& operator=(const ECGroup&) = delete;
    // Delete move constructor and assignment operator
    ECGroup(ECGroup&&) = delete;
    ECGroup& operator=(ECGroup&&) = delete;

    inline ~ECGroup() {
        EC_GROUP_free(group);
        BN_CTX_free(ctx);
    }

    inline const EC_GROUP *get_group() const {
        return group;
    }

    inline EC_GROUP *get_group() {
        return group;
    }

    inline BN_CTX *get_ctx() const {
        return ctx;
    }
};

// A wrapper around OpenSSL EC_POINT to simplify memory management
// Must free EC_POINT BEFORE freeing EC_GROUP
class ECPoint {
private:
    EC_POINT *point = nullptr;
    bool initialized = false;
public:
    inline ECPoint() {}

    inline void init(const ECGroup &ec_group) {
        if (!initialized) {
            point = EC_POINT_new(ec_group.get_group());
            assert(point != nullptr && "Failed to create EC_POINT");
            initialized = true;
        }
    }

    inline ECPoint(ECGroup &ec_group) {
        init(ec_group);
    }
    
    // Delete copy constructor and assignment operator
    ECPoint(const ECPoint&) = delete;
    ECPoint& operator=(const ECPoint&) = delete;
    // Delete move constructor and assignment operator
    ECPoint(ECPoint&&) = delete;
    ECPoint& operator=(ECPoint&&) = delete;

    inline ~ECPoint() {
        if (initialized) {
            EC_POINT_free(point);
        }
    }

    inline EC_POINT *get_point() {
        return point;
    }

    inline const EC_POINT *get_point() const {
        return point;
    }

    inline size_t get_serialization_size(const ECGroup &ec_group) const {
        const EC_GROUP *group = ec_group.get_group();
        BN_CTX *ctx = ec_group.get_ctx();
        return EC_POINT_point2oct(group, point, POINT_CONVERSION_COMPRESSED, NULL, 0, ctx);
    }

    inline void serialize(const ECGroup &ec_group, uint8_t *buf) const {
        const EC_GROUP *group = ec_group.get_group();
        BN_CTX *ctx = ec_group.get_ctx();
        size_t needed_len = EC_POINT_point2oct(group, point, POINT_CONVERSION_COMPRESSED, NULL, 0, ctx);
        size_t len = EC_POINT_point2oct(group, point, POINT_CONVERSION_COMPRESSED, buf, needed_len, ctx);
        assert(len > 0);
    }

    // inline uint8_t *serialize(const ECGroup &ec_group, size_t &len) const {
    //     const EC_GROUP *group = ec_group.get_group();
    //     BN_CTX *ctx = ec_group.get_ctx();
    //     len = EC_POINT_point2oct(group, point, POINT_CONVERSION_COMPRESSED, NULL, 0, ctx);
    //     uint8_t *buf = new uint8_t[len];
    //     EC_POINT_point2oct(group, point, POINT_CONVERSION_COMPRESSED, buf, len, ctx);
    //     return buf;
    // }
};

#endif // EC_HH
