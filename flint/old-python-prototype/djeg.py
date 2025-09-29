from random import randint, randrange
from common import generate_generator, keygen
from damgard_jurik import dlog_z_sp1
from gmpy2 import powmod # type: ignore

# DJEG: Damgard-Jurik-style "ElGamal" cryptosystem
# https://eprint.iacr.org/2025/094.pdf#page=58.99
# Figure 18
# DJEG.Setup(1λ):
# 1 : (p, q) ← GenPQ(1λ)
# 2 : N := pq
# 3 : g0 ←$ Z∗
# N 2
# 4 : g := (g0)2N ∈ Z∗
# N 2
# 5 : return crs := (N, g)
# DJEG.KeyGen(crs):
# 1 : parse crs = (N, g)
# 2 : s ←$ [N ]
# 3 : f := gs
# 4 : (pk, sk) := (f, s)
# 5 : return (pk, sk)
# DJEG.Encrypt(crs, pk, x, w):
# 1 : parse crs = (N, g)
# 2 : parse pk = f
# 3 : r ←$ {0, 1, . . . , N }
# 4 : c0 := gr mod N w+1
# 5 : c1 := (N + 1)xf r mod N w+1
# 6 : return ct := (c0, c1)
# DJEG.FlipEncrypt(crs, pk, x, w):
# 1 : parse crs = (N, g)
# 2 : parse pk = f
# 3 : r ←$ {0, 1, . . . , N }
# 4 : c0 := (N + 1)xgr mod N w+1
# 5 : c1 := f r mod N w+1
# 6 : return ct := (c0, c1)
# DJEG.Decrypt(sk, ct):
# 1 : parse ct = (c0, c1)
# 2 : c′ := c1/(c0)sk
# 3 : x := c′ − 1
# N w
# 4 : return x
def djeg_setup(w_max: int):
    """Return crs"""
    (p, q), N = keygen()
    # g_0 = randrange(0, N ** (w_max + 1))
    # g = powmod(g_0, 2 * (N ** w_max), N ** (w_max + 1))
    g = generate_generator(N, w_max)
    return N, g

def djeg_keygen(crs, w_max: int):
    """Returns (pk, sk). Need w_max which is the maximum w supported"""
    N, g = crs
    s = randrange(0, N)
    f = powmod(g, s, N ** (w_max + 1))
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
    c0 = powmod(g, r, modulus)
    # TODO: perhaps can optimize (N+1)^x mod N^w+1 using some binomial expansion (See page 16 and 17 of https://brics.dk/RS/00/45/BRICS-RS-00-45.pdf)
    c1 = (powmod(N + 1, x, modulus) * powmod(f, r, modulus)) % modulus
    return c0, c1, w

def djeg_encrypt_raw(crs, pk, x, w: int, r: int | None = None):
    c0, c1, _ = djeg_encrypt(crs, pk, x, w, r)
    return c0, c1

# Not checked, will be good to verify what is the correctness definition for this
def djeg_flip_encrypt(crs, pk, x, w: int, r: int | None = None):
    """Returns ct"""
    N, g = crs
    modulus = N ** (w + 1)
    f = pk
    if r is None:
        r = randrange(0, N)
    c0 = (powmod(N + 1, x, modulus) * powmod(g, r, modulus)) % modulus
    c1 = powmod(f, r, modulus)
    return c0, c1, w

def djeg_flip_encrypt_raw(crs, pk, x, w: int, r: int | None = None):
    c0, c1, _ = djeg_flip_encrypt(crs, pk, x, w, r)
    return c0, c1

# For now, we assume w = 1
def djeg_decrypt(crs, sk, ct):
    """Returns x"""
    N, _ = crs
    c0, c1, w = ct
    modulus = N ** (w + 1)
    # N_w = N ** w
    # c_prime = (c1 * inverse_mod(powmod(c0, sk, N**2), modulus)) % modulus
    c_prime = (c1 * powmod(c0, -sk, modulus)) % modulus
    # print('c_prime =', c_prime)
    return dlog_z_sp1(N, w, c_prime)

if __name__ == "__main__":
    print('Initialized')
    w = 5
    crs = djeg_setup(w)
    # print('crs =', crs)
    pk, sk = djeg_keygen(crs, w)
    # print('(pk, sk) =', (pk, sk))

    print('Done with keygen')

    N, g = crs
    # x = 1100
    # x2 = 11
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
    flipped_ct = djeg_flip_encrypt(crs, pk, x, w)
    flipped_pt = djeg_decrypt(crs, sk, flipped_ct)
    # flipped_ct2 = djeg_flip_encrypt(crs, pk, x, w)
    # flipped_pt2 = djeg_decrypt(crs, sk, flipped_ct2)
    # assert flipped_pt == flipped_pt2
    assert (-x * sk) % (N**w) == flipped_pt # A flipped plaintext decrypts to a circular encryption of -x * sk. (Or x * sk if we defined the public key as f = g^{-sk} as in the technical overview of the paper)
    print('All tests passed!')
