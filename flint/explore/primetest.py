# Experiments that guide parameters for next prime function

from gmpy2 import next_prime, mpz, is_prime # type: ignore
import random

def list_primes_under(n: int):
    """List all primes under n using gmpy2"""
    primes = []
    nxt_prime = mpz(2)
    while nxt_prime < n:
        primes.append(nxt_prime)
        nxt_prime = next_prime(nxt_prime)
    return primes

def experiment(prime_ub: int, num_bits: int | None = None, iters: int = 1000000):
    primes = list_primes_under(prime_ub)
    print(sum(1 / p for p in primes))
    print('There are', len(primes), 'primes under', prime_ub)
    print('Iters:', iters)
    print('Number of bits of the starting number:', num_bits)
    if num_bits is None: # Infinite
        start_mods = [random.randrange(prime) for prime in primes]
    else:
        start = random.randint(0, 2 ** num_bits - 1)
        start_mods = [start % prime for prime in primes]
    prev_i = [0] * len(primes)
    num_trial_subs = [0] * len(primes) # Number of trial subtractions needed if doing trial subtraction instead of division
    num_hits = [0] * len(primes) # Number of hits for each prime
    # print(start_mods[:100], primes[:100])
    num_primetests = 0
    possible_primes = []
    for i in range(iters):
        for j in range(len(primes)):
            num = (start_mods[j] + i) % primes[j]
            # num_trial_subs[j] += (i - prev_i[j]) // primes[j]
            # prev_i[j] = i
            num_trial_subs[j] += (i) // primes[j]
            # prev_i[j] = i
            num_hits[j] += 1
            if num == 0:
                break
        else:
            num_primetests += 1
            if num_bits is not None:
                possible_primes.append(start + i)
    ratio = num_primetests / iters
    print('We need', ratio, 'primes to test per integer')
    if ratio != 0:
        print('We need a primality test per', 1 / ratio, 'integers')
    for j, trial_subs in enumerate(num_trial_subs):
        if trial_subs == 0:
            continue
        print('Prime', primes[j], ': ', trial_subs, 'trial subtractions')
        print('Prime', primes[j], ': ', num_hits[j], 'hits')
        print('Prime', primes[j], ': ', trial_subs / num_hits[j], 'subtractions per hit')
    intervals = list(range(0, len(primes), len(primes) // 10))
    for i, j in enumerate(intervals):
        if i == 0:
            continue
        prev_j = intervals[i - 1]
        print('Prime', primes[prev_j], '-', primes[j], ': ', sum(num_trial_subs[k] for k in range(prev_j, j)), 'trial subtractions')
    print(num_trial_subs[:10])
    print(num_hits[:10])
    if num_bits is not None:
        if len(possible_primes) == 0:
            print('[!] No possible primes found')
            return
        print('Testing', len(possible_primes), 'possible primes')
        count_primes = 0
        for possible_prime in possible_primes:
            count_primes += is_prime(possible_prime)
        print(count_primes, 'of', len(possible_primes), 'are primes (', count_primes / len(possible_primes), ')')
        theoretically_minimum_ratio = count_primes / iters
        print('Theoretically minimum primetest count:', theoretically_minimum_ratio)
        if theoretically_minimum_ratio != 0:
            print('Inverse: ', 1 / theoretically_minimum_ratio)

experiment(3000000, num_bits=1536)