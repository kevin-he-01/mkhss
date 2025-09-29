# Verifying the claim in the NIDLS framework paper:
# Lemma 2 + Section 4.1: that phi(g) = g mod N is a coset labeling function
# https://eprint.iacr.org/2022/363.pdf#page=14.64
from sage.all import *
from common import keygen
from damgard_jurik import dlog_z_sp1
import tqdm

# Experiments
def in_coset(N, s, a):
    """Check whether a = (N + 1)^x for some x"""
    result = dlog_z_sp1(N, s, a)
    if result == -1:
        return False
    return power_mod(N + 1, result, N**(s + 1)) == a

# Testing in_coset
# (_, _), N = keygen()
# w = 5
# print(f'{w = }')
# for _ in tqdm.tqdm(range(10)):
#     x = randint(0, N ** w)
#     assert in_coset(N, w, power_mod(N + 1, x, N ** (w + 1)))
# for _ in tqdm.tqdm(range(10000)):
#     a = randint(0, N ** (w + 1))
#     assert not in_coset(N, w, a) # With high probability

# Conjecture: h = g0 % N = g1 % N is in the same coset as g0 and g1
# (_, _), N = keygen()
# w = 5
# print(f'{w = }')
# g0 = randint(1, N ** (w + 1) - 1)
# x = randint(0, N ** w - 1)
# g1 = (g0 * power_mod(N + 1, x, N ** (w + 1))) % N ** (w + 1)
# h0 = g1 % N
# h1 = g0 % N
# assert h0 == h1, "h doesn't match"
# h = h0
# assert in_coset(N, w, (g0 * inverse_mod(h, N ** (w + 1))) % N ** (w + 1)), 'h and g0 not in same coset'
# assert in_coset(N, w, (g1 * inverse_mod(h, N ** (w + 1))) % N ** (w + 1)), 'h and g1 not in same coset'
# print('Passed!')

# Sub-conjecture: For all g, there exists a unique h such that 0 <= h < N and h = g * (N + 1)^x % N^w+1 for some x (i.e., h is in the same coset as g)
# Experimentally verify this
(_, _), N = keygen()
print(f'{N = }')
w = 5
print(f'{w = }')
# for _ in range(10):
#     g = 0
#     while gcd(g, N) != 1:
#         g = randint(1, N ** (w + 1) - 1)
#     print(f'{g = }')
#     cosets = {(g * power_mod(N + 1, x, N ** (w + 1))) % N ** (w + 1) for x in range(N ** w)}
#     assert len(cosets) == N ** w
#     print('Coset size:', len(cosets))
#     less_than_n = [x for x in cosets if x < N]
#     print('Number of elements < N:', len(less_than_n))
#     assert len(less_than_n) == 1
#     print('Unique element < N:', less_than_n[0])
#     assert less_than_n[0] == g % N
# print('Passed!')

# Sub-sub-conjecture: The cosets represented by [0, 1, ..., N-1] are disjoint and cover the entire coset space
# |Z_N^*| = phi(N) = (p - 1)(q - 1), and |Z_N^w+1| = N^w * phi(N), so disjointness implies covering
# i.e., for all g0 != g1 in Z_N^*, g1 * g0^{-1} is not in the coset of 1 (i.e., g1 * g0^{-1} != (N + 1)^x for some x)
def random_invertible(N):
    while True:
        g = randint(1, N - 1)
        if gcd(g, N) == 1:
            return g
for _ in tqdm.tqdm(range(1000)):
    g0 = random_invertible(N)
    g1 = random_invertible(N)
    # print(g0, g1)
    if g0 == g1:
        continue
    assert not in_coset(N, w, (g1 * inverse_mod(g0, N ** (w + 1))) % N ** (w + 1))
print('Passed!')
