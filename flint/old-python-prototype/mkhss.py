# Notation: {{x}}: input share (encryption), <<x>>: memory share (secret share)
import dataclasses
import pprint
import time

from prf import prf
from typing import Any
import random
from random import randint
from ddlog import nidls_ddlog_damgard_jurik
from djeg import djeg_decrypt, djeg_flip_encrypt_raw, djeg_setup
from nim import nim_decode, nim_encode, nim_setup
from gmpy2 import powmod, mpz, mp_version # type: ignore

OPTIMIZE_EXPONENT_SIZE = True
INTERMEDIATE_INTEGER_BOUND = 10 ** 15 # If optimizing exponent size, this is the bound on all intermediate memory shares of the RMS program. Smaller value lead to better performance

SEC_PARAM = 128 # lambda
STAT_SEC_PARAM = 40 # tau
SHORT_EXPONENT_ASSUMPTION = True
SECRET_BITS = 2 * SEC_PARAM # If using short exponent assumption, use 256 bits for 128 bits of security
if SHORT_EXPONENT_ASSUMPTION:
    W_MAX = 1
else:
    W_MAX = 3
SECRET_MODULUS = INTERMEDIATE_INTEGER_BOUND * 2 ** STAT_SEC_PARAM
print('GMP version:', mp_version())
print(f'{SEC_PARAM = }')
print(f'{STAT_SEC_PARAM = }')
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
    # _, c1 = djeg_flip_encrypt_raw((N, g), f, x, W_MAX, r=r)
    c0 = ct[0]
    # c1 = powmod(f, r, N ** (W_MAX + 1))
    # c1_s = powmod(c1, s, N ** (W_MAX + 1))
    c1_s = powmod(f, r * s, N ** (W_MAX + 1))
    x_1 = (c0, c1_s)
    return x_1, 'garbage'

def explinenc_r(crs, sk, pk, x):
    """Returns {{x}}"""
    N = crs[0]
    # ((c0, c1), (c0_prime, c1_prime)) = x
    ((c0, c1), _) = x
    _, s = sk
    x_1 = (c0, powmod(c1, s, N ** (W_MAX + 1)))
    # x_2 = (c0_prime, powmod(c1_prime, s, N ** 2))
    return x_1, 'garbage'

def mkhss_setup():
    """Return crs"""
    N, g = djeg_setup(W_MAX)
    crsnim = nim_setup(N, W_MAX)
    nim_rand_shift = random.randrange(0, N ** W_MAX)
    kprf = random.randbytes(16)
    return N, g, crsnim, nim_rand_shift, kprf

# Please regenerate, W_MAX has changed
# TEST_N = mpz('?')
# TEST_G = mpz('?')
# TEST_H = mpz('?')
# CRS_TEST = (
#     TEST_N,
#     TEST_G,
#     (TEST_N, W_MAX, TEST_G, TEST_H),
#     mpz("?"),
#     b'\x8b\x8e\x8f\x8f\x8e\x8e\x8f\x8e\x8e\x8e\x8e\x8e\x8e\x8e\x8e\x8e',
# )

def mkhss_keygen(crs, sigma: str):
    """Returns (pk, sk)"""
    N, g, crsnim, _, _ = crs
    if SHORT_EXPONENT_ASSUMPTION:
        # 256 bits of secret achieve 128 bits of security due to Pollard's lambda algorithm
        s = randint(1, 2 ** SECRET_BITS - 1) * SECRET_MODULUS + 1
    else:
        s = randint(1, N - 1) * SECRET_MODULUS + 1
    # if sigma == 'A': # Malicious attack: glitch all shares to 0
    #     s = 0
    f = powmod(g, -s, N ** (W_MAX + 1))
    pe, st = nim_encode(crsnim, sigma, s)
    pk = (pe, f)
    sk = (st, s)
    return pk, sk

def mkhss_share(crs, sigma: str, pk_sigma, x, r: int | None = None):
    """Returns ([[x]]^sigma_sigma, [[x]]^sigma_1-sigma)"""
    N, g, crsnim, _, _ = crs
    _, f_sigma = pk_sigma
    if r is None:
        r = randint(0, N)
    r_prime = 'garbage'
    ct = djeg_flip_encrypt_raw((N, g), f_sigma, x, W_MAX, r=r)
    ct_prime = 'garbage'
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

@dataclasses.dataclass
class RMSProgramState:
    # sigma: str
    k_sigma: tuple
    instruction_count: int # nonrepeating counter to ensure independent PRF output (this is id in the paper). Caveat: both parties need to execute instructions in the same order for this to agree

# "Initializes" The RMS program by  doing NIM.Decode. Essentially what the first part of `mkhss_eval` does. Prepare for MKHSS.Eval
def mkhss_init_rms(crs, sigma: str, sk_sigma, pk_1_sigma):
    # TODO: make x_A and x_B arrays
    """"Initializes" The RMS program doing NIM.Decode."""
    N, g, crsnim, nim_rand_shift, kprf = crs
    if SHORT_EXPONENT_ASSUMPTION:
        secret_bits = SECRET_BITS
    else:
        secret_bits = N.bit_length()
    # sigma_int = sigma_to_int(sigma)
    # print(x_i_sigma)
    st_sigma, s_sigma = sk_sigma
    pe_1_sigma, f_1_sigma = pk_1_sigma
    # f = powmod(f_1_sigma, s_sigma, crs[0] ** (W_MAX + 1))
    z_sigma = (nim_decode(crsnim, sigma, pe_1_sigma, st_sigma) + nim_rand_shift) % (N ** W_MAX) # Randomize shares to ensure output is a share over Z (not just N**w) with high probability
    # TODO: Optimize the above with subtraction instead of mod reduction
    # print(f'{z_sigma = }')
    if OPTIMIZE_EXPONENT_SIZE:
        # TODO: for even further optimization, should round the modulus to the next power of 2
        exponent_cutoff = SECRET_MODULUS ** 2 * 2 ** (2 * secret_bits + STAT_SEC_PARAM)
        print('NIM exponent cutoff: ', exponent_cutoff.bit_length())
        z_sigma %= exponent_cutoff # We can bound <z> = <z>_a - <z>_b (using subtractive shares over Z) <= N^2 * 2^(2 * |s|) = B, so reduction modulo B * 2^lambda is safe
    k_sigma = (z_sigma, 1) if sigma == 'A' else (z_sigma, 0) # Memory share: <<1>>_{sigma}
    return RMSProgramState(k_sigma, 0)

# mu is a party identifier, it may or may not be the same as sigma (mu != sigma: A -> B, mu == sigma: A -> A)
# x_mu_sigma: [x]^mu_sigma: share from mu intended for sigma
def mkhss_sync_share(crs, sigma: str, sk_sigma, pk_1_sigma, x_mu_sigma, mu: str):
    if mu == sigma:
        synced_share = explinenc_s(crs, sk_sigma, pk_1_sigma, x_mu_sigma)
    else:
        synced_share = explinenc_r(crs, sk_sigma, pk_1_sigma, x_mu_sigma)
    return synced_share

def mkhss_rms_mult(crs, rms_state: RMSProgramState, input_share_x, mem_share_y):
    """I * M -> M"""
    N, _, _, _, kprf = crs
    if SHORT_EXPONENT_ASSUMPTION:
        secret_bits = SECRET_BITS
    else:
        secret_bits = N.bit_length()
    # id is the instruction identifier
    id = rms_state.instruction_count
    rms_state.instruction_count += 1
    print(f'RMS mult {id = }: Exponent widths: {mem_share_y[0].bit_length()}/{(N ** (W_MAX)).bit_length()} and {mem_share_y[1].bit_length()}/{(N).bit_length()}')
    # TODO: for even further optimization, should round the modulus to the next power of 2
    exponent_mod = (SECRET_MODULUS ** 2 * 2 ** (2 * secret_bits + STAT_SEC_PARAM) * INTERMEDIATE_INTEGER_BOUND)
    print(f'Exponent bit cutoff: {exponent_mod.bit_length()}')

    start = time.time()
    # sigma = rms_state.sigma
    c, _ = input_share_x # second element is garbage
    # start_ip = time.time()
    part1 = powmod(c[0], mem_share_y[0], N ** (W_MAX + 1))
    # mid_ip = time.time()
    part2 = powmod(c[1], mem_share_y[1], N ** (W_MAX + 1))
    f_1_sigma = (part1 * part2) % (N ** (W_MAX + 1))
    # end_ip = time.time()
    # print(f'Time taken for exponent linear inner product {id = }: {end_ip - start_ip}')
    # print(f'part1 time: {mid - start}, part2 time: {end - mid}')
    w = W_MAX
    rand_shift = prf(t=N ** w, key=kprf, data = id.to_bytes(8, 'big'))
    z_s_sigma = (nidls_ddlog_damgard_jurik(N, w, f_1_sigma) + rand_shift) % (N ** w) # TODO: Optimize with subtraction instead of mod reduction
    if OPTIMIZE_EXPONENT_SIZE:
        z_s_sigma %= exponent_mod
    # z_s_sigma: <z * s>_sigma
    # <z>_sigma = <z * s>_sigma mod m_s (s === 1 mod m_s)
    result_mem = tuple([z_s_sigma, z_s_sigma % SECRET_MODULUS])
    end = time.time()
    print(f'RMS mult {id = }: {(end - start) * 1000 :.2f} ms')
    print()
    return result_mem

def mkhss_rms_add(mem_share_x, mem_share_y):
    """M + M -> M"""
    return tuple([mem_share_x[0] + mem_share_y[0], mem_share_x[1] + mem_share_y[1]])

# TODO: RMS add for input shares (need to test this out since many optimizations depend on this)

def mkhss_rms_add_input(crs, input_share_x, input_share_y):
    """I + I -> I"""
    N, _, _, _, _ = crs
    c, _ = input_share_x
    c_prime, _ = input_share_y
    return (((c[0] * c_prime[0]) % N ** (W_MAX + 1), (c[1] * c_prime[1]) % N ** (W_MAX + 1)), 'garbage')
    # return ((c[0] * c_prime[0], c[1] * c_prime[1]), 'garbage')

def mkhss_rms_sub(mem_share_x, mem_share_y):
    """M + M -> M"""
    return tuple([mem_share_x[0] - mem_share_y[0], mem_share_x[1] - mem_share_y[1]])

# Convert input to memory share
def mkhss_rms_convert(crs, rms_state: RMSProgramState, input_share):
    """I -> M"""
    return mkhss_rms_mult(crs, rms_state, input_share, rms_state.k_sigma)

def mkhss_rms_output(crs, mem_share):
    return mem_share[1]

# Test the MKHSS scheme
print('Testing MKHSS scheme...')
print('mkhss_setup()...')
crs = mkhss_setup()
# crs = CRS_TEST
N = crs[0]
print(f'{N = }')
# start = time.time()
print('A: mkhss_keygen()...')
pk_A, sk_A = mkhss_keygen(crs, 'A')
print('B: mkhss_keygen()...')
pk_B, sk_B = mkhss_keygen(crs, 'B')
# end = time.time()
# print('Time taken for keygen:', end - start)
# print(f'pk_A: ')
# pprint.pprint(pk_A)
# print(f'sk_A: ')
# pprint.pprint(sk_A)
# print(f'pk_B: ')
# pprint.pprint(pk_B)
# print(f'sk_B: ')
# pprint.pprint(sk_B)
# print(f'{pk_A = }')
# print(f'{sk_A = }')
# print(f'{pk_B = }')
# print(f'{sk_B = }')


x_A = 12
y_A = -100
x_B = -34
print(f'A: mkhss_share(x_A = {x_A})...')
x_A_A, x_A_B = mkhss_share(crs, 'A', pk_A, x_A)
print(f'A: mkhss_share(y_A = {y_A})...')
y_A_A, y_A_B = mkhss_share(crs, 'A', pk_A, y_A)
print(f'B: mkhss_share(x_B = {x_B})...')
x_B_A, x_B_B = mkhss_share(crs, 'B', pk_B, x_B)

def eval_P_plain(x_A, y_A, x_B):
    # return y_A * (x_A * x_B + y_A)
    return (x_A + x_B) * (y_A + x_B)

def rms_program_eval(rms_st, x_A_synced, y_A_synced, x_B_synced):
    # x_A_mem = mkhss_rms_convert(crs, rms_st, x_A_synced)
    # # x_B_mem = mkhss_rms_convert(crs, rms_st_a, x_B_synced)
    # x_A_x_B_mem = mkhss_rms_mult(crs, rms_st, x_B_synced, x_A_mem)
    # y_A_mem = mkhss_rms_convert(crs, rms_st, y_A_synced)
    # x_A_x_B_plus_yA_mem = mkhss_rms_add(x_A_x_B_mem, y_A_mem)
    # final = mkhss_rms_mult(crs, rms_st, y_A_synced, x_A_x_B_plus_yA_mem)
    # z_sigma = mkhss_rms_output(crs, final)
    # return z_sigma
    x = mkhss_rms_add_input(crs, x_A_synced, x_B_synced)
    # x_mem = mkhss_rms_convert(crs, rms_st, x)
    y_A = mkhss_rms_add_input(crs, y_A_synced, x_B_synced)
    y = mkhss_rms_convert(crs, rms_st, y_A)
    x_mem = mkhss_rms_mult(crs, rms_st, x, y)
    final = mkhss_rms_output(crs, x_mem)
    return final

### Alice
print('A: mkhss_init_rms()...')
rms_st_a = mkhss_init_rms(crs, 'A', sk_A, pk_B)
x_A_synced = mkhss_sync_share(crs, 'A', sk_A, pk_B, x_A_A, 'A')
y_A_synced = mkhss_sync_share(crs, 'A', sk_A, pk_B, y_A_A, 'A')
x_B_synced = mkhss_sync_share(crs, 'A', sk_A, pk_B, x_B_A, 'B')
debug_set_result('A', x_A_synced, x_B_synced)
z_A = rms_program_eval(rms_st_a, x_A_synced, y_A_synced, x_B_synced)

### Bob
print('B: mkhss_init_rms()...')
rms_st_b = mkhss_init_rms(crs, 'B', sk_B, pk_A)
x_A_synced = mkhss_sync_share(crs, 'B', sk_B, pk_A, x_A_B, 'A')
y_A_synced = mkhss_sync_share(crs, 'B', sk_B, pk_A, y_A_B, 'A')
x_B_synced = mkhss_sync_share(crs, 'B', sk_B, pk_A, x_B_B, 'B')
debug_set_result('B', x_B_synced, x_A_synced)
z_B = rms_program_eval(rms_st_b, x_A_synced, y_A_synced, x_B_synced)

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
    assert djeg_decrypt(crs[:2], joint_sk, result_A_A[0] + (W_MAX,)) == (x_A * sk_A[1] * sk_B[1]) % (N ** W_MAX)
    assert djeg_decrypt(crs[:2], joint_sk, result_B_B[0] + (W_MAX,)) == (x_B * sk_A[1] * sk_B[1]) % (N ** W_MAX)
    # <c_1, k>
    if x_A >= N:
        print('\t[!] Warning: x_A >= N')
    if x_B >= N:
        print('\t[!] Warning: x_B >= N')
    # assert djeg_decrypt(crs[:2], joint_sk, result_A_A[1] + (1,)) == x_A
    # assert djeg_decrypt(crs[:2], joint_sk, result_B_B[1] + (1,)) == x_B
    print('\t[v] Exponent-linear decodability verified')

print(f'{z_A = }')
print(f'{z_B = }')
z = z_A - z_B
print('z = z_A - z_B =', z)
expected_z = eval_P_plain(x_A, y_A, x_B)
assert z == expected_z, f'{z} != {expected_z}'
print('[v] Correctness verified')
print('[v] All tests passed!')