from random import randint, randrange
from common import generate_generator, keygen
from damgard_jurik import dlog_z_sp1
from gmpy2 import powmod # type: ignore

# A Length-Flexible Threshold Cryptosystem with Applications
# https://brics.dk/RS/03/16/BRICS-RS-03-16.pdf
# Section 3. The Basic Cryptosystem
S = 0
# The original cryptosystem sets s = 0, but the MKHSS paper sets it to s = 1
# See paper, Section 3.1. Be sure to check more up-to-date references
SMALL_GROUP_EXP = S + 1

# Section 3.1. Theorem 1. "A simple modification makes it possible to show semantic security based only on [DCR], using a technique from [15].""
DCR_ONLY = False

def djeg_setup():
    """Return crs"""
    (p, q), N = keygen()
    g_0 = randrange(1, N ** SMALL_GROUP_EXP)
    g = powmod(g_0, 2 * N ** (SMALL_GROUP_EXP - 1), N ** SMALL_GROUP_EXP)
    return N, g

def djeg_keygen(crs):
    """Returns (pk, sk)."""
    N, g = crs
    s = randrange(0, N)
    f = powmod(g, s, N ** SMALL_GROUP_EXP)
    return f, s

# In the construction, w = 1 or 5
# Inspired by [DJ03, Section 3]: https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=ffa5976cb7c9c3694b5e2e275f07742afa656dd9
# MKHSS->DJ03 variable names: w->s, f->h, s->alpha
def djeg_encrypt(crs, pk, x, w: int, r: int | None = None):
    """Returns ct"""
    N, g = crs
    modulus = N ** (w + 1)
    f = pk
    if r is None:
        r = randrange(0, N)
    if DCR_ONLY:
        g_prime = powmod(g, N ** (w - S), modulus)
        f_prime = powmod(f, N ** (w - S), modulus)
        c0 = powmod(g_prime, r, modulus)
        # TODO: perhaps can optimize (N+1)^x mod N^w+1 using some binomial expansion (See page 16 and 17 of https://brics.dk/RS/00/45/BRICS-RS-00-45.pdf)
        c1 = (powmod(N + 1, x, modulus) * powmod(f_prime, r, modulus)) % modulus
    else:
        c0 = powmod(g, r, N ** SMALL_GROUP_EXP)
        # TODO: perhaps can optimize (N+1)^x mod N^w+1 using some binomial expansion (See page 16 and 17 of https://brics.dk/RS/00/45/BRICS-RS-00-45.pdf)
        c1 = (powmod(N + 1, x, modulus) * powmod((powmod(f, r, N ** SMALL_GROUP_EXP)), N ** (w - S), modulus)) % modulus
        # c1 = (powmod(N + 1, x, modulus) * powmod((powmod(f, r, N ** SMALL_GROUP_EXP) + 42 * N ** SMALL_GROUP_EXP), N ** (w - S), modulus)) % modulus
        # c1 = (powmod(N + 1, x, modulus) * powmod(powmod(f, r, N ** SMALL_GROUP_EXP), N ** (w - 1), modulus)) % modulus # will cause homomorphic operations to fail
    return c0, c1, w

def djeg_encrypt_raw(crs, pk, x, w: int, r: int | None = None):
    c0, c1, _ = djeg_encrypt(crs, pk, x, w, r)
    return c0, c1

# Not checked, will be good to verify what is the correctness definition for this
# def djeg_flip_encrypt(crs, pk, x, w: int, r: int | None = None):
#     """Returns ct"""
#     N, g = crs
#     modulus = N ** (w + 1)
#     f = pk
#     if r is None:
#         r = randrange(0, N)
#     c0 = (powmod(N + 1, x, modulus) * powmod(g, r, modulus)) % modulus
#     c1 = powmod(f, r, modulus)
#     return c0, c1, w

# def djeg_flip_encrypt_raw(crs, pk, x, w: int, r: int | None = None):
#     c0, c1, _ = djeg_flip_encrypt(crs, pk, x, w, r)
#     return c0, c1

# For now, we assume w = 1
def djeg_decrypt(crs, sk, ct):
    """Returns x"""
    N, _ = crs
    c0, c1, w = ct
    modulus = N ** (w + 1)
    # N_w = N ** w
    # c_prime = (c1 * inverse_mod(powmod(c0, sk, N**2), modulus)) % modulus
    if DCR_ONLY:
        c_prime = (c1 * powmod(c0, -sk, modulus)) % modulus
    else:
        c_prime = (c1 * powmod(powmod(c0, sk, N ** SMALL_GROUP_EXP), -N ** (w - S), modulus)) % modulus
        # c_prime = (c1 * powmod(powmod(c0, sk, N ** SMALL_GROUP_EXP), -N ** (w - 1), modulus)) % modulus # will cause homomorphic operations to fail, despite decrypting fresh ciphertexts correctly
    # print('c_prime =', c_prime)
    return dlog_z_sp1(N, w, c_prime)

if __name__ == "__main__":
    print('Initialized')
    print(f'DCR_ONLY = {DCR_ONLY}')
    print('SMALL_GROUP_EXP =', SMALL_GROUP_EXP)
    crs = djeg_setup()
    # print('crs =', crs)
    pk, sk = djeg_keygen(crs)
    # print('(pk, sk) =', (pk, sk))

    print('Done with keygen')

    N, g = crs
    # x = 1100
    # x2 = 11
    w = 5
    x = randint(0, N ** w - 1)
    x2 = randint(0, N ** w - 1)
    ct = djeg_encrypt(crs, pk, x, w)
    ct2 = djeg_encrypt(crs, pk, x2, w)
    sum_ct = (ct[0] * ct2[0] % N**(w + 1), ct[1] * ct2[1] % N**(w + 1), w)
    # print('ct=', ct)
    actual_x = djeg_decrypt(crs, sk, ct)
    assert x == actual_x, f'Expected {x}, got {actual_x}'
    print('Decryption passed')
    sum_pt = djeg_decrypt(crs, sk, sum_ct)
    # print('sum_pt =', sum_pt)
    assert (x + x2) % (N**w) == sum_pt
    print('Homomorphic operations passed!')
    scalar = 13371337043
    scalar_prod_ct = (powmod(ct[0], scalar, N**(w + 1)), powmod(ct[1], scalar, N**(w + 1)), w)
    scalar_prod_pt = djeg_decrypt(crs, sk, scalar_prod_ct)
    assert (x * scalar) % (N**w) == scalar_prod_pt
    print('Scalar multiplication passed!')

    # Flipped encryption test
    # flipped_ct = djeg_flip_encrypt(crs, pk, x, w)
    # flipped_pt = djeg_decrypt(crs, sk, flipped_ct)
    # # flipped_ct2 = djeg_flip_encrypt(crs, pk, x, w)
    # # flipped_pt2 = djeg_decrypt(crs, sk, flipped_ct2)
    # # assert flipped_pt == flipped_pt2
    # assert (-x * sk) % (N**w) == flipped_pt # A flipped plaintext decrypts to a circular encryption of -x * sk. (Or x * sk if we defined the public key as f = g^{-sk} as in the technical overview of the paper)
    print('All tests passed!')
