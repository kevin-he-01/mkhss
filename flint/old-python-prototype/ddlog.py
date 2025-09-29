from random import randint
from common import keygen
from damgard_jurik import dlog_z_sp1
from gmpy2 import powmod, invert # type: ignore

## DDLog over Paillier group
# https://eprint.iacr.org/2021/262.pdf#page=13.40
# Algorithm 3.2
def ddlog_paillier(N, g):
    """N = pq, g \\in \\Z_{N^2}^*"""
    h_prime, h = divmod(int(g), N)
    z = (h_prime * invert(h, N)) % N
    return z

def test_ddlog_paillier():
    (_, _), N = keygen()
    g0 = randint(1, N**2 - 1)
    for _ in range(10):
        x = randint(0, N - 1)
        print('Testing x =', x)
        g1 = (g0 * powmod(N + 1, x, N**2)) % N**2
        x0 = ddlog_paillier(N, g0)
        x1 = ddlog_paillier(N, g1)
        actual_x = (x1 - x0) % N
        assert x == actual_x
    print('DDLog Paillier: All tests passed!')

# Using the NIDLS framework, Section 4.1 instantiation of DDLog:
# https://eprint.iacr.org/2022/363.pdf#page=14.64
def nidls_ddlog_damgard_jurik(N, w, g):
    """N = pq, g \\in \\Z_{N^{w+1}}^*"""
    # Potential optimization: Call ddlog_paillier when w = 1
    # nidls_ddlog_damgard_jurik is empirically verified to be very fast so no need to optimize
    h = g % N # h = phi(g) = g mod N, h is the coset label
    # print(f'{h = }')
    z = dlog_z_sp1(N, w, (g * invert(h, N**(w+1))) % N**(w+1))
    assert z != -1
    # z = (h_prime * inverse_mod(h, N)) % N
    return z

# test_ddlog_paillier()

# Test the NIDLS instantiation
def test_ddlog_damgard_jurik(w: int = 1):
    (_, _), N = keygen()
    print(f'{w = }')
    g0 = randint(1, N**(w + 1) - 1)
    for _ in range(10):
        x = randint(0, N ** w - 1)
        # print('Testing x =', x)
        g1 = (g0 * powmod(N + 1, x, N**(w + 1))) % N**(w + 1)
        x0 = nidls_ddlog_damgard_jurik(N, w, g0)
        x1 = nidls_ddlog_damgard_jurik(N, w, g1)
        actual_x = (x1 - x0) % (N ** w)
        assert x == actual_x
    print('DDLog DJ: All tests passed!')

if __name__ == "__main__":
    # test_ddlog_damgard_jurik(w = 1)
    # test_ddlog_damgard_jurik(w = 2)
    # test_ddlog_damgard_jurik(w = 5)
    (p, q), N = keygen()
    # g = randint(1, N**2 - 1)
    # h = randint(1, N**2 - 1)
    g = randint(1, 10 * p - 1) * powmod(N+1, 5, N**2) % N**2
    h = randint(1, 10 * q - 1) * powmod(N+1, 6, N**2) % N**2
    print((ddlog_paillier(N, g) + ddlog_paillier(N, h)) % N)
    print(ddlog_paillier(N, g * h))
