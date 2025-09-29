from random import randint, randrange

from common import generate_generator, keygen
from ddlog import nidls_ddlog_damgard_jurik
from gmpy2 import powmod # type: ignore


# Definition 16 (Non-Interactive Multiplication; Adapted from [BDSS25]). Let λ be a security parameter, R be a finite ring. A non-interactive multiplication (NIM) scheme consists of three algorithms
# NIM = (Setup, Encode, Decode) with the following syntax:
# – Setup(1λ) → crs. The randomized setup algorithm takes as input the security parameter and
# outputs a common reference string crs.
# – Encode(crs, x) → (peσ , stσ ). The randomized encoding algorithm takes as input the CRS crs and
# a ring element x ∈ R. It outputs a public encoding peσ and secret state stσ .
# – Decode(crs, pe1−σ , stσ ) → ⟨z⟩σ . The deterministic decoding algorithm takes as input the CRS crs,
# another party’s public encoding pe1−σ , and secret state stσ . It outputs a subtractive secret share
# of z.
# The above functionality must satisfy correctness and security, which are defined as follows:
# Correctness. For all security parameters λ ∈ N and every pair of elements x, y ∈ R, a NIM scheme
# is said to be correct if there exists a negligible function negl(·) such that:
# Pr
# 
# 
# 
# 
# 
# 
# 
# 
# ⟨z⟩A − ⟨z⟩B = xy :
# crs ← Setup(1λ)
# (peA, stA) ← Encode(crs, x)
# (peB , stB ) ← Encode(crs, y)
# ⟨z⟩A := Decode(crs, peB , stA)
# ⟨z⟩B := Decode(crs, peA, stB )
# ≥ 1 − negl(λ).

# An insecure mock implementation of the NIM scheme:
# NIM = (Setup, Encode, Decode)
def nim_setup(N: int | None = None, w: int = 1):
    if N is None:
        (_, _), N = keygen()
    # g and h are generated in the same way as in DJEG.Setup, i.e., a 2N-th residue that generates a subgroup of order phi(N)/4
    # TODO: optimize by reusing the same g and just generating a new h
    # g_0 = randrange(0, N ** (w + 1))
    # g = powmod(g_0, 2 * (N ** w), N ** (w + 1))
    g = generate_generator(N, w)
    # h_0 = randrange(0, N ** (w + 1))
    # h = powmod(h_0, 2 * (N ** w), N ** (w + 1))
    h = generate_generator(N, w)
    crs = (N, w, g, h)
    return crs

# NIM from DCR. Let N be a suitable composite modulus and let g and h be random generators of
# Z∗
# N 2 that are part of the CRS. The protocol is instantiated over the ring R = Zℓ, where for correctness
# we need ℓ < 2−λ · √N . The high level idea behind the NIM construction is to have:
# – Alice’s public encoding consist of a Pedersen-like commitment grA hx to her element x and
# – Bob’s public encoding consist of an encryption (grB , (N + 1)y hrB ) of his element y,
# where rA and rB are random elements of ZN .
# Then, given Alice’s encoding peA := grA hx, Bob derives ZB := (grA hx)−rB = g−rArB h−xrB .
# Similarly, given Bob’s encoding peB := (grB , (N + 1)y hrB ), Alice derives ZA := (grB )rA · ((N +
# 1)y hrB )x. It’s not hard to see that ZA and ZB form multiplicative shares of (N + 1)xy mod N since:
# ZA · ZB = ((grB )rA · ((N + 1)y hrB )x) · (grA hx)−rB
# = (grArB · (N + 1)xy hxrB ) · (g−rArB h−xrB )
# = (N + 1)xy .
# Therefore, by applying the DDLog procedure to ZA and ZB , the parties recover subtractive shares
# of xy mod N . Moreover, because x, y < 2−λ · √N , we have that, with all but negligible probability,
# the shares ⟨xy⟩A and ⟨xy⟩B are subtractive shares over the integers, by the correctness of the DDLog
# algorithm
def nim_encode_alice(crs, x):
    N, w, g, h = crs
    r_A = randint(0, N - 1)
    # r_A = 42
    pe_A = (powmod(g, r_A, N ** (w + 1)) * powmod(h, x, N ** (w + 1))) % (N ** (w + 1))
    st_A = (x, r_A)
    return pe_A, st_A

def nim_encode_bob(crs, y):
    N, w, g, h = crs
    r_B = randint(0, N - 1)
    # r_B = 42
    # TODO: can potentially optimize (N+1)^y using some binomial expansion
    pe_B = (powmod(g, r_B, N ** (w + 1)), powmod(N + 1, y, N ** (w+1)) * powmod(h, r_B, N ** (w+1)) % (N ** (w + 1)))
    st_B = r_B
    return (pe_B, st_B)

def nim_encode(crs, sigma: str, x):
    N, w, g, h = crs
    # Insecure scheme
    # r = randint(0, N ** w - 1) if sigma == 'A' else 0
    # return (x, r), (x, r)
    # Secure scheme
    if sigma == 'A':
        return nim_encode_alice(crs, x)
    else:
        assert sigma == 'B'
        return nim_encode_bob(crs, x)

def nim_decode_alice(crs, pe_B, st_A):
    N, w, g, h = crs
    c0, c1 = pe_B
    x, r_A = st_A
    # ZA := (grB )rA · ((N + 1)y hrB )x
    Z_A = powmod(c0, r_A, N ** (w + 1)) * powmod(c1, x, N ** (w + 1)) % (N ** (w + 1))
    return nidls_ddlog_damgard_jurik(N, w, Z_A)

def nim_decode_bob(crs, pe_A, st_B):
    N, w, g, h = crs
    c0 = pe_A
    r_B = st_B
    # ZB := (grA hx)−rB = g−rArB h−xrB
    # To obtain subtractive shares, we compute 1 / Z_B instead of Z_B
    inverse_Z_B = powmod(c0, r_B, N ** (w + 1))
    return nidls_ddlog_damgard_jurik(N, w, inverse_Z_B)

def nim_decode(crs, sigma: str, pe, st):
    # Insecure version:
    # N, w, g, h = crs
    # x, r1 = pe
    # y, r2 = st
    # if sigma == 'A':
    #     return (x * y + r1 + r2) % N ** w
    # else:
    #     return r1 + r2 % N ** w
    # Secure version
    if sigma == 'A':
        return nim_decode_alice(crs, pe, st)
    else:
        assert sigma == 'B'
        return nim_decode_bob(crs, pe, st)

# Test the NIM scheme
if __name__ == "__main__":
    crs = nim_setup(w = 5)
    x = randint(0, 50)
    y = randint(0, 50)
    peA, stA = nim_encode(crs, 'A', x)
    peB, stB = nim_encode(crs, 'B', y)
    zA = nim_decode(crs, 'A', peB, stA)
    zB = nim_decode(crs, 'B', peA, stB)
    print(f'{zA = }')
    print(f'{zB = }')
    assert (zA - zB) == (x * y), f'Expected {x} * {y} = {x * y}, got {zA - zB}'
    print("NIM scheme is correct")
