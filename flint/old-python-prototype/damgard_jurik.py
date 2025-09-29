import math
from random import randrange
from gmpy2 import powmod, invert, lcm # type: ignore
from common import keygen
# https://brics.dk/RS/00/45/BRICS-RS-00-45.pdf
# Damgård-Jurik cryptosystem

# See page 4-5 of the paper linked above
# This equation leads to the following algorithm:
# i := 0;
# for j:= 1 to s do
# begin
# t1 := L(a mod nj+1);
# t2 := i;
# for k:= 2 to j do
# begin
# i := i − 1;
# t2 := t2 ∗ i mod nj ;
# t1 := t1 − t2∗n^{k−1}/k! mod nj ;
# end
# i := t1;
# end
# define the function L() by L(b) = (b − 1)/n
# See also: https://github.com/cryptovoting/damgard-jurik/blob/5471ec2eb098381dd4dc37fac6b041a010290960/damgard_jurik/crypto.py#L277
def dlog_z_sp1(n: int, s: int, a: int) -> int:
    """Computes discrete log of a w.r.t. base (1+n) in Z_{n^{s+1}}. It is sometimes referred to as the L_s(.) function."""
    i = 0
    for j in range(1, s+1):
        t1, rem = divmod((a % (n**(j+1))) - 1, n)
        if rem != 0:
            return -1
        # t1 = ((a % (n**(j+1))) - 1) // n
        t2 = i
        for k in range(2, j+1):
            i = i - 1
            t2 = (t2 * i) % (n**j)
            # assert (t2 * n**(k-1)) % math.factorial(k) == 0 # Not true for s >= 3
            # t1 = (t1 - (t2 * n**(k-1) // math.factorial(k))) % (n**j)
            # Use modular inverse, not division here
            t1 = (t1 - (t2 * n**(k-1) * invert(math.factorial(k), n**j))) % (n**j)
        i = t1
    return i # We are guaranteed that i is already reduced modulo n^s since it is assigned from t1, which is reduced mod n^j and j <= k <= s

def in_coset(N, s, a):
    """Check whether a = (N + 1)^x for some x"""
    result = dlog_z_sp1(N, s, a)
    if result == -1:
        return False
    return powmod(N + 1, result, N**(s + 1)) == a

def test_dlog_z_sp1(s: int, prime_len: int = 16):
    (_, _), n = keygen(prime_len=prime_len)
    assert s >= 1
    print('n =', n)
    print('s =', s)
    for _ in range(10):
        x = randrange(0, n**s) # Order of n + 1 is n^s.
        # print('x =', x)
        a = powmod((n + 1), x, n**(s+1))
        actual_x = dlog_z_sp1(n, s, a)
        assert actual_x == x, f'{actual_x} != {x}'
    print('All tests passed!')

# The original Damgård-Jurik paper cryptosystem
# https://brics.dk/RS/00/45/BRICS-RS-00-45.pdf#page=10.32
def dj_keygen(w: int):
    (p, q), N = keygen()
    # Can choose j and x randomly
    # j = randrange(1, N ** w)
    # g0 = randrange(1, N ** (w + 1))
    # x = powmod(g0, N ** w, N ** (w + 1))
    # Also works (and is still secure):
    j = 1
    x = 1
    g = (powmod((N + 1), j, N ** (w + 1)) * x) % (N ** (w + 1))
    lamb = lcm(p-1, q-1)
    d = lamb
    pk = (N, w, g)
    sk = (d, j)
    return pk, sk

def dj_encrypt(pk, m):
    N, w, g = pk
    r = randrange(1, N ** (w + 1))
    c = (powmod(g, m, N ** (w + 1)) * powmod(r, N ** w, N ** (w + 1))) % (N ** (w + 1))
    return c

def dj_decrypt(pk, sk, c, w_override: int | None = None):
    d, j = sk
    N, w, g = pk
    if w_override is not None:
        w = w_override
    m = (dlog_z_sp1(N, w, powmod(c, d, N ** (w + 1))) * invert(j * d, N ** w)) % (N ** w)
    return m

if __name__ == "__main__":
    # test_dlog_z_sp1(s = 1) # s = 1 is Paillier
    # test_dlog_z_sp1(s = 2)
    # test_dlog_z_sp1(s = 3)
    # test_dlog_z_sp1(s = 6)
    # test_dlog_z_sp1(s = 7)
    # test_dlog_z_sp1(s = 13)

    # test_dlog_z_sp1(s = 1, prime_len=1536)
    # test_dlog_z_sp1(s = 5, prime_len=1536)

    print('Testing Damgård-Jurik encryption and decryption...')
    # Experiment 1
    # pk, sk = dj_keygen(w=1)
    # N, _, g = pk
    # w = 5
    # # m = randrange(0, 100)
    # m = randrange(0, N)
    # mult = N ** (w - 1)
    # print('m =', m)
    # c = dj_encrypt(pk, m)
    # # c = powmod(c, N, N ** (w + 1))
    # m3 = dj_decrypt(pk, sk, c)
    # # m4 = dj_decrypt(pk, sk, powmod(c, mult, N ** (w + 1)))
    # c += randrange(0, N ** (w - 1)) * N ** 2 # Still works even if we add a multiple of N^2
    # m2 = dj_decrypt(pk, sk, c, w_override=w)
    # m5 = dj_decrypt(pk, sk, powmod(c, mult, N ** (w + 1)), w_override=w)
    # assert m == m3, f'{m} != {m3}'
    # print('w = ', w)
    # print('m = ', m)
    # print('m2 =', m2)
    # print('m2 // N', m2 // (N))
    # mred = m2 % N
    # print('m2 (reduced mod N) =', mred)
    # assert mred == m, f'{mred} != {m}'
    # # assert m4 % mult == 0
    # # actual_m = m4 // mult
    # # assert actual_m == m
    # assert m5 % mult == 0
    # actual_m = m5 // mult
    # print('actual_m =', actual_m)
    # print('c =', c)

    # Experiment 2
    # w = 5
    # pk, sk = dj_keygen(w=w-1)
    # N, _, g = pk
    # # m = randrange(0, 100)
    # m = randrange(0, N ** (w - 1))
    # mult = N
    # print('m =', m)
    # c = dj_encrypt(pk, m)
    # # c = powmod(c, N, N ** (w + 1))
    # m3 = dj_decrypt(pk, sk, c)
    # # m4 = dj_decrypt(pk, sk, powmod(c, mult, N ** (w + 1)))
    # c += randrange(0, N) * N ** w # Still works even if we add a multiple of N^w
    # m2 = dj_decrypt(pk, sk, c, w_override=w)
    # m5 = dj_decrypt(pk, sk, powmod(c, mult, N ** (w + 1)), w_override=w)
    # assert m == m3, f'{m} != {m3}'
    # print('w = ', w)
    # print('m = ', m)
    # print('m2 =', m2)
    # print('m2 //  N^(w-1)', m2 // (N ** (w - 1)))
    # mred = m2 % (N ** (w - 1))
    # print('m2 (reduced mod N^(w-1)) =', mred)
    # assert mred == m, f'{mred} != {m}'
    # # assert m4 % mult == 0
    # # actual_m = m4 // mult
    # # assert actual_m == m
    # assert m5 % mult == 0
    # actual_m = m5 // mult
    # print('actual_m =', actual_m)
    # print('c =', c)

    # Experiment 3
    # # Bijection between k1 and k2 in ((1+N)^x + k1 * N^(w+1))^y === (1 + N)^{xy + k2 * N^w} mod N^s
    # # (_, _), n = keygen(prime_len=64)
    # n = 259285936971992382932621424347725005529
    # w = 1
    # s = 4
    # x = 1
    # # y = 80
    # y = 1
    # print('n =', n)
    # print('s =', s)
    # print('w =', w)
    # print('x =', x)
    # print('y =', y)
    # for k1 in [0, 1, 2, 3, 4, n - 1]:
    #     print('k1 =', k1)
    #     lhs = powmod(powmod((n + 1), x, n**(s)) + k1 * (n**(w + 1)), y, (n**(s)))
    #     rhs_exponent = dlog_z_sp1(n, s - 1, lhs)
    #     assert rhs_exponent % (n**w) == x * y, f'{rhs_exponent} != {x * y}'
    #     k2 = rhs_exponent // (n**w)
    #     print('k2 =', k2)

    # Bijection between k1 and k2 in (g^x + k1 * N^(w+1)) === g^{x + k2 * N^w} mod N^s
    # For g = (1 + N)^j where gcd(j, N) = 1
    # (_, _), n = keygen(prime_len=64)
    n = 259285936971992382932621424347725005529
    w = 1
    s = 4
    x = 1
    # y = 80
    y = 1
    print('n =', n)
    print('s =', s)
    print('w =', w)
    print('x =', x)
    print('y =', y)
    for k1 in [0, 1, 2, 3, 4, n - 1]:
        print('k1 =', k1)
        lhs = powmod(powmod((n + 1), x, n**(s)) + k1 * (n**(w + 1)), y, (n**(s)))
        assert in_coset(n, s - 1, lhs), f'{lhs} not in coset'
        rhs_exponent = dlog_z_sp1(n, s - 1, lhs)
        assert rhs_exponent % (n**w) == x * y, f'{rhs_exponent} != {x * y}'
        k2 = rhs_exponent // (n**w)
        print('k2 =', k2)

    # (_, _), n = keygen(prime_len=16)
    # print('n =', n)
    # s = 5
    # print('s =', s)
    # x = n // 2
    # # print(dlog_z_sp1(n, s, powmod((n + 1), x, n**(s+1)))) # Good, =x
    # val = powmod((n + 1), x, n**2)
    # dl = dlog_z_sp1(n, s, val) # Unknown
    # print(dl // (n ** 2))
    # print(dl // n)
    # print(dl % n)
    # assert dl % n == x
    # assert powmod((n + 1), dl, n**(s+1)) == val
