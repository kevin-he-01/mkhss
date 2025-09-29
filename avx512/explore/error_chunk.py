# Problem: binary string of length n, e bit errors in random shuffled locations, chunked into b bins of size n / b.
# f_{n,t,b}(e) = Probability that there exists at least one bin whose number of errors exceed t / b
# Experiments to estimate f(e)

import random
from typing import Callable
import numpy as np
from scipy.stats import beta

def get_ci(k: int, n: int, alpha: float = 0.05):
    p_u, p_o = beta.ppf([alpha/2, 1 - alpha/2], [k, k + 1], [n - k + 1, n - k])
    if np.isnan(p_o):
        p_o=1
    if np.isnan(p_u):
        p_u=0
    return [p_u, p_o]

def ci2str(p_u, p_o):
    p_u = round(p_u, 3)
    p_o = round(p_o, 3)
    return '({}, {})'.format(p_u, p_o)

def report_prob_success(trial: Callable[[], bool], invert=False, samples: int=10000):
    success = 0
    for _ in range(samples):
        success += (invert != trial())
    # return success / samples
    return success / samples, get_ci(success, samples)

def f(n, t, b, e):
    assert n % b == 0
    assert t % b == 0
    k = n // b
    tt = t // b
    bits = [1] * e + [0] * (n - e)
    random.shuffle(bits)
    for chunk_start in range(0, n, k):
        chunk = bits[chunk_start:chunk_start + k]
        if sum(chunk) >= tt:
            return True
    return False

n = 2048
t = 256 * 3
# t = 256 * 4
for b in [4, 8, 16, 32]:
    start = int(t * 0.5)
    end = int(t)
    for e in range(start, end, int((end - start) / 20)):
        prob, (p_u, p_o) = report_prob_success(lambda: f(n, t, b, e), samples=1000)
        print('n={}, t={}, b={}, e={}, prob={:.3f} {}'.format(n, t, b, e, prob, ci2str(p_u, p_o)))
