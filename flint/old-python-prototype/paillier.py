from sage.all import *

def keygen(prime_len=512):
    """Outputs a public key (n, g) and a private key (lamb, mu) for the Paillier cryptosystem."""
    p = random_prime(2**prime_len, lbound=2**(prime_len - 1))
    q = random_prime(2**prime_len, lbound=2**(prime_len - 1))
    # print(p, q)
    n = p * q
    lamb = lcm(p - 1, q - 1)
    mu = None
    while mu is None:
        g = randint(1, n**2 - 1)
        assert (power_mod(g, lamb, n**2) - 1) % n == 0
        mu_inv = (power_mod(g, lamb, n**2) - 1) // n
        if gcd(mu_inv, n) == 1:
            mu = inverse_mod(mu_inv, n)
    return (n, g), (lamb, mu)

def keygen_alt(prime_len=512): # Simplified version in the Boneh-Shoup book
    """Outputs a public key (n, g) and a private key (lamb, mu) for the Paillier cryptosystem."""
    p = random_prime(2**prime_len, lbound=2**(prime_len - 1))
    q = random_prime(2**prime_len, lbound=2**(prime_len - 1))
    # print(p, q)
    n = p * q
    lamb = lcm(p - 1, q - 1)
    # lamb = (p - 1) * (q - 1) // 2
    g = n + 1 # Fix g to be n + 1
    # mu_inv = (power_mod(g, lamb, n**2) - 1) // n
    # assert mu_inv == lamb
    # assert gcd(mu_inv, n) == 1
    # mu = inverse_mod(mu_inv, n)
    mu = inverse_mod(lamb, n)
    return (n, g), (lamb, mu)


def enc(pk, m):
    n, g = pk
    assert 0 <= m < n
    r = randint(1, n - 1)
    return power_mod(g, m, n**2) * power_mod(r, n, n**2) % n**2

def dec(pk, sk, c):
    n, _ = pk
    lamb, mu = sk
    return (((power_mod(c, lamb, n**2) - 1) // n) * mu) % n

(n, g), (lamb, mu) = keygen()
print(f'pk=(n,g)={(n, g)}, sk=(lambda, mu)={(lamb, mu)}')
for pt in [42, 1337, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
    ct = enc((n, g), pt)
    pt_actual = dec((n, g), (lamb, mu), ct)
    assert pt_actual == pt
