# Note: this is missing an important experimental optimization. See mkhss-optimized.py for the optimized version
# Notation: {{x}}: input share (encryption), <<x>>: memory share (secret share)
import time
from prf import prf
from typing import Any
import random
from random import randint
from ddlog import nidls_ddlog_damgard_jurik
from djeg import djeg_decrypt, djeg_encrypt_raw, djeg_flip_encrypt_raw, djeg_setup
from nim import nim_decode, nim_encode, nim_setup
from gmpy2 import powmod # type: ignore

OPTIMIZE_EXPONENT_SIZE = True
INTERMEDIATE_INTEGER_BOUND = 10 ** 15 # If optimizing exponent size, this is the bound on all intermediate memory shares of the RMS program. Smaller value lead to better performance

SEC_PARAM = 128 # lambda
SHORT_EXPONENT_ASSUMPTION = True
SECRET_BITS = 2 * SEC_PARAM # If using short exponent assumption, use 256 bits for 128 bits of security
if SHORT_EXPONENT_ASSUMPTION:
    W_MAX = 3
else:
    W_MAX = 5
print(f'{SEC_PARAM = }')
print(f'{SHORT_EXPONENT_ASSUMPTION = }')
print(f'{W_MAX = }')
print('#### Optimizations:')
print(f'{OPTIMIZE_EXPONENT_SIZE = }')
if OPTIMIZE_EXPONENT_SIZE:
    print(f'{INTERMEDIATE_INTEGER_BOUND = }')

verify_exponent_linear = True

def sigma_to_int(sigma: str):
    """Returns 0 if sigma is 'A' and 1 if sigma is 'B'"""
    assert sigma in ('A', 'B')
    return 0 if sigma == 'A' else 1

# ExpLinEncS(skσ , pk1−σ , [[x]]^σ_σ ):
# 1 : parse [[x]]^σ_σ = ((x, r, r′), (_,_))
# 2 : parse skσ := (_, sσ )
# 3 : parse pk1−σ = (_, f1−σ )
# 4 : (c0, c1) := DJEG.FlipEncrypt(f1−σ , x, w; r)
# 5 : (c′
# 0, c′
# 1) := DJEG.Encrypt(f1−σ , x, 2; r′)
# 6 : {{x}} := ((c0, (c1)sσ ), (c′
# 0, (c′
# 1)sσ ))
# 7 : return {{x}}
# ExpLinEncR(skσ , pk1−σ , [[x]]^σ_(1-σ) ):
# 1 : parse [[x]]^σ_(1-σ) = ((c0, c1), (c′
# 0, c′
# 1))
# 2 : parse skσ = ( , sσ )
# 3 : {{x}} := ((c0, (c1)sσ ), (c′
# 0, (c′
# 1)sσ ))
# 4 : return {{x}}
def explinenc_s(crs, sk, pk, x):
    """Returns {{x}}"""
    N, g, _, _, _ = crs
    # ((x, r, r_prime), (_, _)) = x
    ((x, r, r_prime), (ct, ct_prime)) = x
    _, s = sk
    _, f = pk
    # Unoptimized version:
    # c0, c1 = djeg_flip_encrypt_raw((N, g), f, x, W_MAX, r=r)
    # c0_prime, c1_prime = djeg_encrypt_raw((N, g), f, x, 1, r=r_prime)
    # Optimization: Use the existing ciphertexts
    _, c1 = djeg_flip_encrypt_raw((N, g), f, x, W_MAX, r=r)
    _, c1_prime = djeg_encrypt_raw((N, g), f, x, 1, r=r_prime)
    c0 = ct[0]
    c0_prime = ct_prime[0]
    x_1 = (c0, powmod(c1, s, N ** (W_MAX + 1)))
    x_2 = (c0_prime, powmod(c1_prime, s, N ** 2))
    return x_1, x_2

def explinenc_r(crs, sk, pk, x):
    """Returns {{x}}"""
    N = crs[0]
    ((c0, c1), (c0_prime, c1_prime)) = x
    _, s = sk
    x_1 = (c0, powmod(c1, s, N ** (W_MAX + 1)))
    x_2 = (c0_prime, powmod(c1_prime, s, N ** 2))
    return x_1, x_2

# Public Parameters. Let Ssk := {i · N + 1 | 1 ≤ i ≤ N − 1} be the secret key space and let B be a bound
# on the message space. Let NIM = (Setup, Encode, Decode) be a NIM scheme. We will use the algorithms
# ExpLinEncS and ExpLinEncR defined in Figure 20.
# MKHSS.Setup(1λ, w):
# 1 : (N, g) ← DJEG.Setup(1λ)
# 2 : crsnim ← NIM.Setup(1λ)
# 3 : kprf
# 1 , kprf
# 2 ←$ {0, 1}λ
# 4 : crs := (N, g, crsnim, kprf
# 1 , kprf
# 2 )
# 5 : return crs
# MKHSS.KeyGen(crs):
# 1 : parse (N, g, crsnim) from crs
# 2 : s ←$ Ssk, f := g−s
# 3 : (pe, st) ← NIM.Encode(crsnim, s)
# 4 : pk := (pe, f )
# 5 : sk := (st, s)
# 6 : return (pk, sk)
# MKHSS.Share(crs, σ, pkσ , x):
# 1 : parse crs = (N, g, crsnim)
# 2 : parse pkσ = (peσ , fσ )
# 3 : r, r′ ←$ ZN
# 4 : ct ← DJEG.FlipEncrypt(fσ , x, 5; r)
# 5 : ct′ ← DJEG.Encrypt(fσ , x, 1; r′)
# 6 : [[x]]^σ_σ := ((x, r, r′), (ct, ct′))
# 7 : [[x]]^σ_(1-σ) := := (ct, ct′)
# 8 : return ([[x]]^σ_σ, [[x]]^σ_(1-σ))
# MKHSS.Eval(crs, σ, skσ , pk1−σ , JxAKA
# σ , JxB KB
# σ , P ):
# 1 : parse (crsnim, kprf
# 1 , kprf
# 2 ) from crs
# 2 : parse skσ = (stσ , sσ )
# 3 : parse pk1−σ = (pe1−σ , f1−σ )
# 4 : f := (f1−σ )sσ
# 5 : ⟨z⟩σ := NIM.Decode(crsnim, pe1−σ , stσ )
# 6 : kσ := (⟨z⟩σ , 1) if σ = A else kσ := (⟨z⟩σ , 0)
# 7 : for i ∈ [m] :
# 8 : {{x(i)
# σ }} := ExpLinEncS(skσ , pk1−σ , Jx(i)
# σ Kσ
# σ )
# 9 : {{x(i)
# 1−σ }} := ExpLinEncR(skσ , pk1−σ , Jx(i)
# 1−σ K1−σ
# σ )
# 10 : ekσ := (kprf
# 1 , kprf
# 2 , kσ )
# 11 : {{x}} := ({{x(1)
# A }}, . . . , {{x(m)
# A }}, {{x(1)
# B }}, . . . , {{x(m)
# B }})
# 12 : return DEval(σ, ekσ , {{x}}, P )
def mkhss_setup():
    """Return crs"""
    N, g = djeg_setup(W_MAX)
    crsnim = nim_setup(N, W_MAX)
    nim_rand_shift = random.randrange(0, N ** W_MAX)
    kprf = random.randbytes(16)
    return N, g, crsnim, nim_rand_shift, kprf

def mkhss_keygen(crs, sigma: str):
    """Returns (pk, sk)"""
    N, g, crsnim, _, _ = crs
    if SHORT_EXPONENT_ASSUMPTION:
        # 256 bits of secret achieve 128 bits of security due to Pollar's rho
        s = randint(1, 2 ** SECRET_BITS - 1) * N + 1
    else:
        s = randint(1, N - 1) * N + 1
    f = powmod(g, -s, N ** (W_MAX + 1))
    pe, st = nim_encode(crsnim, sigma, s)
    pk = (pe, f)
    sk = (st, s)
    return pk, sk

def mkhss_share(crs, sigma: str, pk_sigma, x):
    """Returns ([[x]]^sigma_sigma, [[x]]^sigma_1-sigma)"""
    N, g, crsnim, _, _ = crs
    _, f_sigma = pk_sigma
    r = randint(0, N)
    r_prime = randint(0, N)
    ct = djeg_flip_encrypt_raw((N, g), f_sigma, x, W_MAX, r=r)
    ct_prime = djeg_encrypt_raw((N, g), f_sigma, x, 1, r=r_prime)
    x_sigma_sigma = (x, r, r_prime), (ct, ct_prime)
    x_sigma_1_sigma = ct, ct_prime
    if sigma == 'A':
        return x_sigma_sigma, x_sigma_1_sigma
    elif sigma == 'B':
        return x_sigma_1_sigma, x_sigma_sigma
    else:
        raise ValueError(f'Invalid party identifier: {sigma}')

result_A_A: Any = None
result_A_B: Any = None
result_B_A: Any = None
result_B_B: Any = None

def debug_set_result(sigma, input_share_sigma, input_share_1_sigma):
    global result_A_A, result_A_B, result_B_A, result_B_B
    if sigma == 'A':
        result_A_A, result_A_B = input_share_sigma, input_share_1_sigma
    elif sigma == 'B':
        result_B_B, result_B_A = input_share_sigma, input_share_1_sigma
    else:
        raise ValueError(f'Invalid party identifier: {sigma}')

# DEval(σ, ekσ , ({{x1}}, . . . , {{xm}}), P ) :
# Parse ekσ = (kprf
# 1 , kprf
# 2 , ⟨⟨1⟩⟩σ ).
# For each id ∈ [|P |], evaluate the id-th instruction as follows:
# • Convert : Mx ← Ix :
# 1: Execute the Mult({{x}}, ⟨⟨1⟩⟩σ ) instruction to compute ⟨⟨x⟩⟩σ .
# • Mult : Mxy ← Ix · My :
# 1: Parse {{x}} = (c1, . . . , cℓ).
# 2: For i ∈ [ℓ] :
# 2.1: f (i)
# σ := ci, ⟨⟨y⟩⟩σ .
# 2.2: ⟨zi⟩σ := DDLog(f (i)
# σ ; F2(kprf
# 2 , id∥i) mod t.
# 3: ⟨⟨xy⟩⟩σ := (⟨z1⟩σ , . . . , ⟨zi⟩σ ).
# • Add : Mx+y ← Mx + My :
# 1: ⟨⟨x + y⟩⟩σ := ⟨⟨x⟩⟩σ + ⟨⟨y⟩⟩σ .
# • Output : z ← Mz :
# 1: Parse ⟨⟨z⟩⟩σ = (⟨z1⟩σ , . . . , ⟨zi⟩σ ).
# 2: Return ⟨zℓ⟩
def distributed_eval(N, g, sigma, ek_sigma, input_shares: list, P):
    (_, kprf, memshare_1) = ek_sigma
    # DEBUG: just run a multiplication assuming input_shares is a list of size 2: [{{x_A}}, {{x_B}}]
    input_share_xA, input_share_xB = input_shares
    # print(f'{memshare_1 = }')
    # print(f'{input_share_xA = }')
    # print(f'{input_share_xB = }')
    def mult(input_share_x, mem_share_y, id: int):
        # id is the instruction identifier
        # Parse {{x}} = (c1, . . . , cℓ)
        c = input_share_x
        # For i ∈ [ℓ]:
        ws = [W_MAX, 1]
        xy_sigma = []
        # print(f'id: {id}, Exponent widths: {mem_share_y[0].bit_length()}/{(N ** (W_MAX)).bit_length()} and {mem_share_y[1].bit_length()}/{(N).bit_length()}')
        for i, ci in enumerate(c):
            w = ws[i]
            # f(i)_sigma := <ci, ⟨⟨y⟩⟩sigma>
            # <.>: Exponent linear inner product
            # start = time.time()
            f_i_sigma = (powmod(ci[0], mem_share_y[0], N ** (w + 1)) * powmod(ci[1], mem_share_y[1], N ** (w + 1))) % (N ** (w + 1))
            # end = time.time()
            # print(f'Time taken for exponent linear inner product {id = }, {i = }: {end - start}')
            # ⟨zi⟩sigma := DDLog(f(i)_sigma) + F2(kprf, id∥i) mod t
            rand_shift = prf(t=N ** w, key=kprf, data = id.to_bytes(8, 'big') + i.to_bytes(8, 'big'))
            zi_sigma = (nidls_ddlog_damgard_jurik(N, w, f_i_sigma) + rand_shift) % (N ** w) # TODO: Optimize with subtraction instead of mod reduction
            if OPTIMIZE_EXPONENT_SIZE:
                # TODO: for even further optimization, should round the modulus to the next power of 2
                if i == 0:
                    if SHORT_EXPONENT_ASSUMPTION:
                        zi_sigma %= (N ** 2 * 2 ** (2 * SECRET_BITS + SEC_PARAM) * INTERMEDIATE_INTEGER_BOUND) # We can bound <z> = <z>_a - <z>_b (using subtractive shares over Z) <= N^2 * 2^(2 * |s|) = B, so reduction modulo B * 2^lambda is safe
                    else:
                        zi_sigma %= (N ** 4 * 2 ** (SEC_PARAM) * INTERMEDIATE_INTEGER_BOUND)
                elif i == 1:
                    zi_sigma %= (2 ** (SEC_PARAM) * INTERMEDIATE_INTEGER_BOUND)
            # zi_sigma %= (N ** (w - 1) * 10 ** 10)
            xy_sigma.append(zi_sigma)
        # print(f'{xy_sigma = }')
        return tuple(xy_sigma)
    start = time.time()
    memshare_xA = mult(input_share_xA, memshare_1, 0)
    # memshare_xB = mult(input_share_xB, memshare_1)
    memshare_xA_xB = mult(input_share_xB, memshare_xA, 1)
    memshare_xA_2_xB = mult(input_share_xA, memshare_xA_xB, 2)
    end = time.time()
    print('Time taken for 3 mults:', end - start)
    memshare_z = memshare_xA_2_xB
    # memshare_z = memshare_xA
    return memshare_z[-1]

def mkhss_eval(crs, sigma: str, sk_sigma, pk_1_sigma, x_A_sigma, x_B_sigma, P):
    # TODO: make x_A and x_B arrays
    """Returns DEval(sigma, ek_sigma, {{x}}, P)"""
    sigma_int = sigma_to_int(sigma)
    x_i_sigma = [x_A_sigma, x_B_sigma]
    # print(x_i_sigma)
    N, g, crsnim, nim_rand_shift, kprf = crs
    st_sigma, s_sigma = sk_sigma
    pe_1_sigma, f_1_sigma = pk_1_sigma
    f = powmod(f_1_sigma, s_sigma, crs[0] ** (W_MAX + 1))
    z_sigma = (nim_decode(crsnim, sigma, pe_1_sigma, st_sigma) + nim_rand_shift) % (N ** W_MAX) # Randomize shares to ensure output is a share over Z (not just N**w) with high probability
    # TODO: Optimize the above with subtraction instead of mod reduction
    # print(f'{z_sigma = }')
    if OPTIMIZE_EXPONENT_SIZE:
        # TODO: for even further optimization, should round the modulus to the next power of 2
        if SHORT_EXPONENT_ASSUMPTION:
            z_sigma %= (N ** 2 * 2 ** (2 * SECRET_BITS + SEC_PARAM)) # We can bound <z> = <z>_a - <z>_b (using subtractive shares over Z) <= N^2 * 2^(2 * |s|) = B, so reduction modulo B * 2^lambda is safe
        else:
            z_sigma %= (N ** 4 * 2 ** (SEC_PARAM))
    k_sigma = (z_sigma, 1) if sigma == 'A' else (z_sigma, 0) # Memory share: <<1>>_{sigma}
    ek_sigma = (nim_rand_shift, kprf, k_sigma)
    # input_share_sigma = {{x_sigma}}
    input_share_sigma = explinenc_s(crs, sk_sigma, pk_1_sigma, x_i_sigma[sigma_int])
    input_share_1_sigma = explinenc_r(crs, sk_sigma, pk_1_sigma, x_i_sigma[1 - sigma_int])
    debug_set_result(sigma, input_share_sigma, input_share_1_sigma) # DEBUG: Make sure both parties receive the same input shares
    if sigma == 'A':
        input_shares = [input_share_sigma, input_share_1_sigma]
    elif sigma == 'B':
        input_shares = [input_share_1_sigma, input_share_sigma]
    # DEBUG, for now: input_shares = [{{x_A}}, {{x_B}}]
    return distributed_eval(N, g, sigma, ek_sigma, input_shares, P) # TODO: allow the application to apply functions to the output shares themselves

# An MKHSS scheme satisfies the following correctness and security properties.
# Correctness. An MKHSS scheme is said to be ε-correct, for some ε ∈ [0, 1), if for all λ ∈ N, all
# 2m-input programs P ∈ P, and all xA, xB ∈ Mm, we have
# Pr
# 
# 
# 
# 
# 
# 
# ⟨z⟩R
# A − ⟨z⟩R
# B̸
# =
# P (xA, xB )
# :
# crs ← Setup(1λ)
# (pkσ , skσ ) ← KeyGen(crs), , ∀σ ∈ {A, B}
# (Jxσ Kσ
# A, Jxσ Kσ
# B ) ← Share(crs, σ, pkσ , xσ ), ∀σ ∈ {A, B}
# ⟨z⟩R
# σ ← Eval(crs, σ, skσ , pk1−σ , JxAKA
# σ , JxB KB
# σ , P ), ∀σ ∈ {A, B}
# 
# 
# 
# 
# 
#  ≤ ε + negl(λ),
# where xσ = (x(1)
# σ , . . . , x(m)
# σ )

# Test the MKHSS scheme
print('==== UNOPTIMIZED VERSION ====')
print('Testing MKHSS scheme...')
print('mkhss_setup()...')
crs = mkhss_setup()
N = crs[0]
print(f'{N = }')
print('A: mkhss_keygen()...')
pk_A, sk_A = mkhss_keygen(crs, 'A')
print('B: mkhss_keygen()...')
pk_B, sk_B = mkhss_keygen(crs, 'B')
x_A = 12345
x_B = 56789
print(f'A: mkhss_share(x_A = {x_A})...')
x_A_A, x_A_B = mkhss_share(crs, 'A', pk_A, x_A)
print(f'B: mkhss_share(x_B = {x_B})...')
x_B_A, x_B_B = mkhss_share(crs, 'B', pk_B, x_B)
P = None # TODO: implement P
def eval_P_plain(x_A, x_B):
    return x_A * x_A * x_B
    # return x_A

print('A: mkhss_eval()...')
z_A = mkhss_eval(crs, 'A', sk_A, pk_B, x_A_A, x_B_A, P)
# Expect: result_A_A = result_B_A = {{x_A}}
# Expect: result_A_B = result_B_B = {{x_B}}
print('B: mkhss_eval()...')
z_B = mkhss_eval(crs, 'B', sk_B, pk_A, x_A_B, x_B_B, P)

if verify_exponent_linear:
    print('Verifying exponent-linear decodability...')
    assert result_A_A is not None
    assert result_A_B is not None
    assert result_B_A is not None
    assert result_B_B is not None
    assert result_A_A == result_B_A
    assert result_A_B == result_B_B
    print('\t[v] {{x}} is equal for both parties')

    # Verify that the shares are exponenti-linear decodable (See page 63 of correctness proof)
    # joint sk := s_A * s_B
    joint_sk = -sk_A[1] * sk_B[1] # Negative sign because djeg_decrypt uses -sk
    # <c_0, k>
    # (k_1, k_2) = (s_A * s_B, 1)
    if x_A * sk_A[1] * sk_B[1] >= N ** W_MAX:
        print('\t[!] Warning: x_A * s_A * s_B >= N ** W_MAX')
    if x_B * sk_A[1] * sk_B[1] >= N ** W_MAX:
        print('\t[!] Warning: x_B * s_A * s_B >= N ** W_MAX')
    assert djeg_decrypt(crs[:2], joint_sk, result_A_A[0] + (W_MAX,)) == (x_A * sk_A[1] * sk_B[1])
    assert djeg_decrypt(crs[:2], joint_sk, result_B_B[0] + (W_MAX,)) == (x_B * sk_A[1] * sk_B[1])
    # <c_1, k>
    if x_A >= N:
        print('\t[!] Warning: x_A >= N')
    if x_B >= N:
        print('\t[!] Warning: x_B >= N')
    assert djeg_decrypt(crs[:2], joint_sk, result_A_A[1] + (1,)) == x_A
    assert djeg_decrypt(crs[:2], joint_sk, result_B_B[1] + (1,)) == x_B
    print('\t[v] Exponent-linear decodability verified')

print(f'{z_A = }')
print(f'{z_B = }')
z = z_A - z_B
print('z = z_A - z_B =', z)
expected_z = eval_P_plain(x_A, x_B)
assert z == expected_z, f'{z} != {expected_z}'
print('[v] Correctness verified')
print('[v] All tests passed!')
print('==== UNOPTIMIZED VERSION ====')
