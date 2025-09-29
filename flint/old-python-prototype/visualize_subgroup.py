# Visualize the subgroups of Z_{n^(s+1)} via CRT decomposition
# This allows empirically verifying many properties of the group, especially what happens
# when you reduce a number mod n, mod n^2, mod n^3, etc.
# or when you extend a number mod n^2 to number mod n^3, "zero-padding" high order bits.

# Properties to test:
# Reduction mod n: H subgroup is unchanged. G is undefined
# - This follows from element mod n being the NIDLS "representative element" of a coset
# - Used in the proof of additive homomorphism of the length-flexible Damgard-Jurik encryption scheme
# (suspected) Reduction mod n^2: H subgroup is unchanged. G mod n is unchanged
# Expansion from mod n^w to n^s: H subgroup is undefined. G mod n^(w-1) is unchanged
# - Used in the proof of security of Damgard-Jurik from DCR.

# Motivation: to get a deeper understanding of all the possible efficient operations (beyond trivial ones applicable to a generic group, e.g., multiplication) on Z_{n^(s+1)}, a very structured group
# so I can find optimizations in cryptosystems

from common import keygen
from damgard_jurik import dlog_z_sp1
from sage.all import discrete_log, Mod
from random import randrange, seed
from gmpy2 import gcd, powmod, invert # type: ignore

s = 1

# (p, q), n = keygen(prime_len=6)
# g0 = 0
# while gcd(g0, n) != 1:
#     g0 = randrange(1, n ** (s+1))
# g = powmod(g0, n ** s, n ** (s+1))

# p =  3356845487
# q =  2672170703
# n = 8970064164859167361
# g =  613444066233896609994442350869584093509432251326150551411399830571830001492

p =  47
q =  59
n = 2773
g =  7540652

print('p = ', p)
print('q = ', q)
print('n =', n)
print('g = ', g)
lamb = (p - 1) * (q - 1) // 2 # lambda = lcm(p-1, q-1) because p and q are safe primes (p = 2q + 1)
lamb_inv = invert(lamb, n ** s)
print('lamb(n) =', lamb)
assert powmod(g, lamb, n ** (s + 1)) == 1

# Format: For group element h = (1+N)^j * (g0^(N^s))^i, decompose as
# "N-adic representation of j" x (i mod p-1, i mod q-1). Maybe even further decompose i mod p-1 since p-1=2q for prime q.

def extract_j_gi(h):
    """Returns (j, g^i), when h = (1+N)^j * g^i"""
    # h = (1+N)^j * g^i
    a = powmod(h, lamb, n ** (s + 1))
    l = dlog_z_sp1(n, s, a)
    assert l != -1
    j = l * lamb_inv % (n ** s)
    first_term = powmod(n + 1, j, n ** (s + 1))
    gi = (h * invert(first_term, n ** (s + 1))) % (n ** (s + 1))
    return (j, gi)

def extract_j_i(h):
    j, gi = extract_j_gi(h)
    return (j, discrete_log(Mod(gi, n), Mod(g, n)))

def n_adic_decompose(j):
    """Returns the N-adic decomposition of j, i.e., (j mod n, j // n mod n, j // n^2 mod n, ...)"""
    result = []
    for _ in range(s):
        result.append(j % n)
        j //= n
    return result

def n_adic_decompose_string(j):
    result = n_adic_decompose(j)
    return ' + '.join([f'{result[i]} * N^{i}' for i in range(len(result))][::-1])

def visualize_element(h):
    try:
        j, i = extract_j_i(h)
        return f'(1+N)^{{{n_adic_decompose_string(j)}}} * g^{i} mod N^{s+1}'
    except ValueError: # TODO: handle quadratic nonresidues better (C2 x C2)
        j, gi = extract_j_gi(h)
        return f'(1+N)^{{{n_adic_decompose_string(j)}}} * {gi} mod N^{s+1}'

# expected_j = randrange(0, n**s)
# print('expected j =', expected_j)
# expected_i = randrange(0, lamb // 2)
# h = (powmod(1 + n, expected_j, n**(s+1)) * powmod(g, expected_i, n**(s+1))) % (n**(s+1))
# print('h =', h)
# j, i = extract_j_i(h)
# # j, gi = extract_j_gi(h)
# print('j =', j)
# # print('gi =', gi)
# assert j == expected_j
# print('i =', i)
# assert i == expected_i
# print('All tests passed!')

# for h in range(1, 100):
#     print('h =', h, '=', visualize_element(h))

for w in range(1, s + 2):
    print(f'g mod N^{w} =', visualize_element(g % (n ** w)))
    j, _ = extract_j_gi(g % (n ** w))
    g_prime = (g * powmod(1 + n, j, n ** (s + 1))) % (n ** (s + 1))
    print(g_prime)
    assert g_prime == g % (n ** w)
    expected_g = ((g % (n ** w)) * powmod(1 + n, -j, n ** (s + 1))) % (n ** (s + 1))
    assert g == expected_g, f'Expected {expected_g}, got {g}'
    # print(f'g mod N^{w} =', visualize_element(g_prime))

# for j in range(3):
#     for i in range(3):
#         print('j = ', j, 'i =', i)
#         h = (powmod(1 + n, j, n**(s+1)) * powmod(g, i, n**(s+1))) % (n**(s+1))
#         for w in range(1, s + 2):
#             print(f'h mod N^{w} =', visualize_element(h % (n ** w)))

# h = (powmod(1 + n, 1337, n**(s+1)) * powmod(g, 42, n**(s+1))) % (n**(s+1))
# for w in range(1, s + 2):
#     print(f'h mod N^{w} =', visualize_element(h % (n ** w)))

# h = (powmod(1 + n, n + 1337, n**(s+1)) * powmod(g, 42, n**(s+1))) % (n**(s+1))
# for w in range(1, s + 2):
#     print(f'h mod N^{w} =', visualize_element(h % (n ** w)))

# for w in range(2, s + 2):
#     h = (powmod(1 + n, 1337, n**w) * powmod(g, 42, n**w)) % (n**w)
#     print(f'h mod N^{w} => (mod N^{s+1})', visualize_element(h))
