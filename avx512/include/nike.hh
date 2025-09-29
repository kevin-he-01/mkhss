#ifndef NIKE_HH
#define NIKE_HH

// Implementation of an elliptic-curve-based Non-Interactive Key Exchange (NIKE) protocol
#include "ec.hh"
#include "bignum_openssl.hh"
#include "random_oracle.hh"

namespace nike {
    constexpr size_t SHARED_KEY_LEN = random_oracle::OUTPUT_LENGTH; // Length of the shared key in bytes

    struct crs {
        ECGroup group;
    };
    struct public_key {
        ECPoint point;
    };
    struct private_key {
        BIGNUMClass s;
    };

    inline void setup(crs &crs) {}

    inline void keygen(const crs &crs, private_key &sk, public_key &pk) {
        BIGNUMClass order;
        pk.point.init(crs.group);
        assert(EC_GROUP_get_order(crs.group.get_group(), order.get_bn(), crs.group.get_ctx()) == 1);
        // Generate a random private key
        assert(BN_rand_range(sk.s.get_bn(), order.get_bn()) == 1);
        // Compute the public key point
        assert(EC_POINT_mul(crs.group.get_group(), pk.point.get_point(), sk.s.get_bn(), NULL, NULL, crs.group.get_ctx()) == 1);
    }

    inline size_t get_public_key_size(const crs &crs, public_key &pk) {
        return pk.point.get_serialization_size(crs.group);
    }

    inline void keyder(const crs &crs, const private_key &sk, const public_key &pk, uint8_t shared_key_bytes[SHARED_KEY_LEN]) {
        ECPoint shared_key;
        shared_key.init(crs.group);
        // Compute the shared key point
        assert(EC_POINT_mul(crs.group.get_group(), shared_key.get_point(), NULL, pk.point.get_point(), sk.s.get_bn(), crs.group.get_ctx()) == 1);
        // Print the public key in hex form
        // char *hex = NULL;
        // hex = EC_POINT_point2hex(crs.group.get_group(), shared_key.get_point(), POINT_CONVERSION_COMPRESSED, crs.group.get_ctx());
        // std::cout << "Shared Key: " << hex << std::endl;
        // OPENSSL_free(hex);
        size_t buflen = shared_key.get_serialization_size(crs.group);
        uint8_t *buf = new uint8_t[buflen + 1];
        shared_key.serialize(crs.group, buf);
        buf[buflen] = random_oracle::MAGIC_NIKE; // Append magic byte to ensure independent randomness
        // Hash the shared key point to get the shared key
        random_oracle::call(buf, buflen + 1, shared_key_bytes);
        delete[] buf;
    }
};

#endif // NIKE_HH
