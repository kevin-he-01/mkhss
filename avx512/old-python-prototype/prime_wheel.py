# Wheel factorization
from sage.all import *
from gmpy2 import mpz, is_prime # type: ignore

num_tests = 0
num_mods = 0

# Relative to primality test
MOD_COST = 1/3000

def is_prime_(n: int) -> bool:
    global num_tests
    num_tests += 1
    # Check if n is prime using the Miller-Rabin test
    return is_prime(n)

def is_2q_1_safe_prime(q: int) -> bool:
    # Check if p is a safe prime
    return is_prime_(q) and is_prime_(2 * q + 1)

# TODO: exclude 2 

SAFE_PRIMES_ONLY = True

def primes_with_product_under(n: int, start: int = 2):
    # Generate primes until their product exceeds n
    primes = []
    product = 1
    p = start
    while True:
        if product * p > n:
            break
        primes.append(p)
        product *= p
        p = next_prime(p)
    return primes

# TODO: Might be able to speed up (finding just the possible primes) via precomputing CRT bases w.r.t. each prime, skipping 0 and (p-1)/2 mod p, and enumerating all combinations
def build_wheel(max_size: int = 1000, safe_primes_only: bool = SAFE_PRIMES_ONLY):
    primes = primes_with_product_under(max_size)
    actual_size = 1
    for p in primes:
        actual_size *= p
    print('Safe primes only:', safe_primes_only)
    print('Wheel size:', actual_size)
    print('Primes:', primes)
    wheel = [0] * actual_size
    for i in range(actual_size):
        wheel[i] = any(i % p == 0 for p in primes)
        # TODO: consider safe primes
        if safe_primes_only:
            wheel[i] = wheel[i] or any((2*i+1) % p == 0 for p in primes if p > 2)
    print('Statistics:')
    print('  Wheel size:', actual_size)
    print('  Number of composites:', wheel.count(True))
    print('  Number of possible primes:', wheel.count(False))
    print('  Possible prime ratio:', wheel.count(False) / float(actual_size))
    print()
    return wheel, primes

def build_wheel_list(max_size: int = 1000, safe_primes_only: bool = SAFE_PRIMES_ONLY):
    # Generate a list of prime offsets using the wheel method
    prime_offsets = []
    wheel, _ = build_wheel(max_size, safe_primes_only)
    for i in range(len(wheel)):
        if not wheel[i]:
            prime_offsets.append(i)
    return prime_offsets

def next_safe_prime_naive(n: int) -> int:
    # Find the next safe prime greater than n
    if n % 2 == 0:
        n += 1
    while True:
        if is_2q_1_safe_prime(n):
            return n
        n += 2

def next_safe_prime_wheel_slow(n: int, wheel: list, divisors: list[int] = []) -> int:
    global num_tests, num_mods
    # Find the next safe prime greater than n using a wheel
    # Start at wheel[n%len(wheel)]
    index = n % len(wheel)
    while True:
        if wheel[index % len(wheel)]: # Composite for sure
            n += 1
            index += 1
            continue
        to_continue = False
        for divisor in divisors:
            # TODO: use one GMP % not two
            num_mods += 1
            if n % divisor == 0 or (2 * n + 1) % divisor == 0:
                n += 1
                index += 1
                to_continue = True
                break
        if to_continue:
            continue
        if is_2q_1_safe_prime(n):
            return n
        n += 1
        index += 1

# start = 2
# for _ in range(10):
#     prime_list = primes_with_product_under(2 ** 64, start=start)
#     print(prime_list)
#     start = prime_list[-1] + 1

# for i in range(1, 7):
#     print('Building wheel with max size:', 10 ** i)
#     wheel = build_wheel(10 ** i)
#     primes_ = []
#     for p in range(len(wheel)):
#         if not wheel[p]:
#             primes_.append(p)
#     if len(primes_) <= 100:
#         print('Primes:', primes_)
#     else:
#         print('Primes:', primes_[:10], '...', primes_[-10:])

NUM_64_BIT_MS = 10

NAIVE = False

# wheel, wheel_primes = build_wheel(10 ** 6)
# wheel, wheel_primes = build_wheel(10 ** 5)
# wheel, wheel_primes = build_wheel(10 ** 4)
# wheel, wheel_primes = build_wheel(10 ** 3)
# wheel, wheel_primes = build_wheel(10 ** 2)
# wheel, wheel_primes = build_wheel(10)
# trial_divisors = []
# last = wheel_primes[-1]

last = 5 - 1
print('constexpr uint64_t divisor_product_table[] = {')
for _ in range(NUM_64_BIT_MS):
    divisors_ = primes_with_product_under(2 ** 64, start=next_prime(last))
    # print('Divisors:', divisors_)
    product_ = 1
    for d in divisors_:
        product_ *= d
    # print('Product:', product)
    print(str(product_)+'ULL,')
    last = divisors_[-1]
print('};')

last = 5 - 1
print('std::vector<uint64_t> divisor_table[] = {')
for _ in range(NUM_64_BIT_MS):
    divisors_ = primes_with_product_under(2 ** 64, start=next_prime(last))
    # print('Divisors:', divisors_)
    product_ = 1
    for d in divisors_:
        product_ *= d
    # print('Product:', product)
    # Print this as a C++ vector
    print('{', end='')
    for d in divisors_:
        print(d, end=', ')
    print('},')
    last = divisors_[-1]
    # print(str(product_)+'ULL,')
print('};')


# print('NUM_64_BIT_MS:', NUM_64_BIT_MS)
# for _ in range(NUM_64_BIT_MS):
#     divisors_ = primes_with_product_under(2 ** 64, start=next_prime(last))
#     # print('Divisors:', divisors_)
#     product_ = 1
#     for d in divisors_:
#         product_ *= d
#     trial_divisors.extend(divisors_)
#     last = trial_divisors[-1]

# print('Trial divisors:', trial_divisors)
# for _ in range(5):
#     # n = mpz(randint(0, 2 ** 1536))
#     n = mpz(randint(0, 2 ** 256))
#     print('Finding next safe prime after', n)
#     num_tests = 0
#     if NAIVE:
#         naive = next_safe_prime_naive(n)
#         print('Naive:', num_tests, 'tests')
#         num_tests_naive = num_tests

#     num_tests = 0
#     wheel_slow = next_safe_prime_wheel_slow(n, wheel)
#     print('Wheel:', num_tests, 'tests')
#     num_tests_wheel = num_tests

#     num_tests = 0
#     num_mods = 0
#     wheel_divisors = next_safe_prime_wheel_slow(n, wheel, trial_divisors)
#     print('Wheel divisors:', num_tests, 'tests', num_mods, 'mods')
#     num_tests_wheel_divisors = num_tests
#     effective_num_tests_wheel_divisors = num_tests + num_mods * MOD_COST

#     if NAIVE:
#         assert naive == wheel_slow, f'Naive: {naive}, Wheel slow: {wheel_slow}'
#     assert wheel_divisors == wheel_slow, f'Wheel divisors: {wheel_divisors}, Wheel slow: {wheel_slow}'
#     print('Next safe prime:', wheel_slow)
#     if NAIVE:
#         print('Speedup:', num_tests_naive / num_tests_wheel)
#         print('Speedup (with divisors):', num_tests_naive / num_tests_wheel_divisors)
#         print('Effective Speedup (with divisors):', num_tests_naive / effective_num_tests_wheel_divisors)
#     print('Relative speedup:', num_tests_wheel / num_tests_wheel_divisors)
#     print('Effective relative speedup:', num_tests_wheel / effective_num_tests_wheel_divisors)
#     print()
