# Simulate exponentiations via addition chains.
# This allows rapid benchmarking of exponentiation algorithms.

# lg_B: number of bits in the exponent

from dataclasses import dataclass
import functools
import math
import pprint
import random
from typing import Callable
from matplotlib import pyplot as plt

@dataclass
class PerformanceMetric:
    num_squares: int
    num_mults: int

    def reset(self):
        """Reset the performance metric."""
        self.num_squares = 0
        self.num_mults = 0

@dataclass
class ConcreteCostModel:
    cost_per_square: float = 1
    cost_per_mult: float = 1

    def cost(self, metric: PerformanceMetric) -> float:
        """Calculate the concrete cost of an exponentiation."""
        return metric.num_squares * self.cost_per_square + metric.num_mults * self.cost_per_mult

def print_cost(metric: PerformanceMetric):
    """Print the cost of an exponentiation."""
    print(f"Squares: {metric.num_squares}, Mults: {metric.num_mults}, Total: {metric.num_squares + metric.num_mults}")

def print_concrete_cost(metric: PerformanceMetric, model: ConcreteCostModel):
    """Print the concrete cost of an exponentiation."""
    cost = model.cost(metric)
    print(f'Concrete cost: {cost}')

@dataclass
class SimulatedElement:
    exponent: int
    metric: PerformanceMetric

    def __mul__(self, other: 'SimulatedElement') -> 'SimulatedElement':
        """Simulate multiplication of two elements."""
        if not (self.exponent == 0 or other.exponent == 0): # multiplication by 1 is free
            if self.exponent == other.exponent:
                self.metric.num_squares += 1
            else:
                self.metric.num_mults += 1
        return SimulatedElement(self.exponent + other.exponent, self.metric)

@dataclass
class SimulatedVectorElement:
    exponent: list[int]
    metric: PerformanceMetric

    def __mul__(self, other: 'SimulatedVectorElement') -> 'SimulatedVectorElement':
        """Simulate multiplication of two elements."""
        if self.exponent == other.exponent:
            self.metric.num_squares += 1
        else:
            self.metric.num_mults += 1
        return SimulatedVectorElement([self.exponent[i] + other.exponent[i] for i in range(len(self.exponent))], self.metric)


def square_and_multiply_mean(lg_B):
    """Naive square-and-multiply exponentiation."""
    return PerformanceMetric(
        num_squares=lg_B,
        num_mults=lg_B // 2,
    )

def square_and_multiply(lg_B):
    exponent_bits = [random.randint(0, 1) for _ in range(lg_B)]
    return PerformanceMetric(
        num_squares=lg_B,
        num_mults=exponent_bits.count(1),
    )

# Unoptimized algorithm
def __brauer_naive(exponent: int, lg_B: int, k: int, base, identity):
    """Brauer's algorithm for exponentiation."""
    # print('Exponent:', exponent)

    # Precompute 2^k-size window
    current = base
    precomputed_cache = [base]
    for i in range(2, 2 ** k):
        if i % 2 == 0:
            current = precomputed_cache[i // 2 - 1] * precomputed_cache[i // 2 - 1]
        else:
            current = current * base
        # assert current.exponent == i, "Exponentiation error in Brauer's algorithm"
        precomputed_cache.append(current)
    # Now exponentiate
    bits_lsb_first = [int(x) for x in bin(exponent)[2:].zfill(lg_B)][::-1]
    current = identity
    first = True
    for window in reversed(range(0, lg_B, k)):
        # Keep squaring until we reach the next window
        if not first:
            for _ in range(k):
                current = current * current
        first = False
        bits = bits_lsb_first[window:window + k]
        # Convert bits to integer
        window_value = sum(bit * (2 ** i) for i, bit in enumerate(bits))
        # print(window_value, bits)
        if window_value != 0:
            current = current * precomputed_cache[window_value - 1]
    return current

def __brauer(exponent: int, lg_B: int, k: int, base, identity):
    """Brauer's algorithm for exponentiation."""
    # print('Exponent:', exponent)

    # Precompute 2^k-size window
    current = base
    precomputed_cache = [base]
    base_sq = base * base
    for i in range(3, 2 ** k, 2):
        current = current * base_sq
        # assert current.exponent == i, f'expected {i}, got {current.exponent}'
        precomputed_cache.append(current)
    # Now exponentiate
    bits_lsb_first = [int(x) for x in bin(exponent)[2:].zfill(lg_B)][::-1]
    current = identity
    first = True
    for window in reversed(range(0, lg_B, k)):
        bits = bits_lsb_first[window:window + k]
        # Convert bits to integer
        window_value = sum(bit * (2 ** i) for i, bit in enumerate(bits))
        # print(window_value, bits)
        if window_value != 0:
            value = window_value
            post_squares = 0
            while value % 2 == 0:
                value //= 2
                post_squares += 1
            if not first:
                for _ in range(k - post_squares):
                    current = current * current
            current = current * precomputed_cache[(value - 1) // 2]
            for _ in range(post_squares):
                current = current * current
        else:
            # Keep squaring until we reach the next window
            if not first:
                for _ in range(k):
                    current = current * current
        first = False
        
    return current

def __sliding_window(exponent: int, lg_B: int, k: int, base, identity):
    """Sliding window exponentiation. Algorithm 14.85 from Handbook of Applied Cryptography."""
    # print('Exponent:', exponent)

    # Precompute 2^k-size window
    current = base
    precomputed_cache = [base]
    base_sq = base * base
    for i in range(3, 2 ** k, 2):
        current = current * base_sq
        # assert current.exponent == i, f'expected {i}, got {current.exponent}'
        precomputed_cache.append(current)
    # Now exponentiate
    e = [int(x) for x in bin(exponent)[2:].zfill(lg_B)][::-1] # Bits (LSB first)
    A = identity
    i = lg_B - 1
    while i >= 0:
        if e[i] == 0:
            # Square for the next bit
            A = A * A
            i -= 1
            continue
        # Find the longest e_{i} e_{i-1} ... e_{l} such that i - l + 1 <= k and e_{l} = 1
        l = max(0, i + 1 - k)
        while e[l] == 0:
            l += 1
        for _ in range(i - l + 1):
            A = A * A
        # Now we have A = A^{2^{i-l+1}} and we can multiply by the precomputed value
        bits = e[l:i+1]  # e_{i} e_{i-1} ... e_{l}
        # Convert bits to integer
        window_value = sum(bit * (2 ** i) for i, bit in enumerate(bits))
        assert window_value % 2 == 1, f'expected odd window value, got {window_value}'
        A = A * precomputed_cache[(window_value - 1) // 2]
        i = l - 1
        
    return A

def brauer_analytical_cost_simple(lg_B: int, k: int) -> int:
    """Calculate the analytical cost of the unoptimized Brauer's algorithm which includes squares."""
    # Cost is lg_B + ceil(lg_B / k) + 2^k - 2
    return lg_B + math.ceil(lg_B / k) + ((2 ** k) - 2) + 1

def brauer_analytical_cost(lg_B: int, k: int) -> int:
    """Calculate the analytical cost of Brauer's algorithm."""
    # Cost is lg_B + ceil(lg_B / k) + 2^k - 2
    return lg_B + math.ceil(lg_B / k) + ((2 ** k) - 2) // 2 + 1

def brauer_concrete_cost(lg_B: int, k: int, model: ConcreteCostModel) -> float:
    """Calculate the analytical cost of Brauer's algorithm."""
    # Cost is lg_B + ceil(lg_B / k) + 2^k - 2
    return lg_B * model.cost_per_square + math.ceil(lg_B / k) * model.cost_per_mult + (((2 ** k) - 2) // 2) * model.cost_per_mult + 1 * model.cost_per_square

def brauer_optimize_k(lg_B: int, cost_function: Callable[[int, int], float] = brauer_analytical_cost) -> int:
    k = 1
    prevcost = float('inf')
    while True:
        cost = cost_function(lg_B, k)
        if cost >= prevcost:
            break
        prevcost = cost
        k += 1
    return k - 1  # Return the last k that was better

def optimal_brauer_table(ub: int = 10000, cost_function: Callable[[int, int], float] = brauer_analytical_cost):
    """Generate a table of optimal k values for Brauer's algorithm."""
    # Generate bitcnt_table
    bitcnt_table = []
    prev_k = 1
    for i in range(1, ub + 1):
        new_k = brauer_optimize_k(i, cost_function=cost_function)
        if new_k > prev_k:
            bitcnt_table.append(i)
            prev_k = new_k
    return bitcnt_table

bitcnt_table_sliding = [7,25,81,241,673,1793,4609,11521,28161]

bitcnt_table = [4, 15, 52, 165, 486, 1351, 3592, 9225]

def brauer_k(lg_B: int, table: list[int] = bitcnt_table) -> int:
    """Return the k value for Brauer's algorithm based on lg_B."""
    k = 0
    while lg_B > table[k]:
        k += 1
    return k + 1  # +1 because we want the next power of two

def _brauer(exponent: int, lg_B: int, base, identity):
    k = brauer_k(lg_B)
    # print(f"Using Brauer's algorithm with k={k} for lg_B={lg_B}")
    return __brauer(exponent, lg_B, k, base, identity)

def brauer(lg_B: int):
    """Brauer's algorithm for exponentiation."""
    metric = PerformanceMetric(num_squares=0, num_mults=0)
    exponent = random.randrange(2 ** lg_B)
    base = SimulatedElement(1, metric)
    identity = SimulatedElement(0, metric)
    result = _brauer(exponent, lg_B, base, identity)
    print('Expected cost: ', brauer_analytical_cost(lg_B, brauer_k(lg_B)))
    print('Actual cost: ', metric.num_squares + metric.num_mults)
    assert result.exponent == exponent, f'expected {exponent}, got {result.exponent}'
    return metric

def _sliding_window(exponent: int, lg_B: int, base, identity):
    k = brauer_k(lg_B, table=bitcnt_table_sliding)
    return __sliding_window(exponent, lg_B, k, base, identity)

def sliding_window(lg_B: int):
    """Sliding window exponentiation."""
    metric = PerformanceMetric(num_squares=0, num_mults=0)
    exponent = random.randrange(2 ** lg_B)
    base = SimulatedElement(1, metric)
    identity = SimulatedElement(0, metric)
    result = _sliding_window(exponent, lg_B, base, identity)
    # print('Expected cost: ', brauer_analytical_cost(lg_B, brauer_k(lg_B)))
    # print('Actual cost: ', metric.num_squares + metric.num_mults)
    assert result.exponent == exponent, f'expected {exponent}, got {result.exponent}'
    return metric

# Fixed-base windowing method
class FBEWindow:
    precomputed_powers: list[SimulatedElement]

    def __init__(self, identity: SimulatedElement, base: SimulatedElement, lg_B: int, k: int):
        """Precompute powers of the base for fixed-base exponentiation."""
        # Compute g, g^{2^k}, g^{2^{2*k}}, ...
        self.base = base
        self.identity = identity
        self.precomputed_powers = [base]
        self.lg_B = lg_B
        self.k = k
        current = base
        t_plus_1 = lg_B // k + (lg_B % k != 0)
        for _ in range(1, t_plus_1):
            for _ in range(k):
                current = current * current
            self.precomputed_powers.append(current)
    
    def compute(self, exponent: int) -> SimulatedElement:
        """Compute the exponentiation using the precomputed powers."""
        bits_lsb_first = [int(x) for x in bin(exponent)[2:].zfill(self.lg_B)][::-1]
        chunks = [sum(bit * (2 ** j) for j, bit in enumerate(bits_lsb_first[i:i + self.k])) for i in range(0, len(bits_lsb_first), self.k)]
        A = self.identity
        B = self.identity
        for j in range(1 << self.k, 0, -1):
            for i, chunk in enumerate(chunks):
                if chunk == j:
                    # We have a match, compute the exponentiation
                    B = B * self.precomputed_powers[i]
            A = A * B
            # TODO: micro optimization: don't multiply if A or B is already 1
        return A

# Fixed-base comb method
class FBEComb:
    precomputed_powers: list[SimulatedElement]

    # t is like log(B) here (not log(B)/k)
    def __init__(self, identity: SimulatedElement, base: SimulatedElement, t: int, h: int, v: int, a: int, b: int):
        """Precompute powers of the base for fixed-base exponentiation."""
        self.identity = identity
        gi = [base]
        for i in range(1, h):
            last = gi[-1]
            for _ in range(a):
                last = last * last # Intermediate powers are wasted
            gi.append(last)
        assert len(gi) == h
        # v * 2^h table
        self.G = [[None] * (2**h) for _ in range(v)]
        for i in range(1, 2**h):
            # Bit decompose i
            i_i = [int(x) for x in bin(i)[2:].zfill(h)][::-1]
            self.G[0][i] = functools.reduce(lambda x, y: x * y, [gi[j] if i_i[j] else identity for j in range(h)], identity)
            for j in range(1, v):
                last = self.G[j-1][i]
                for _ in range(b):
                    last = last * last
                self.G[j][i] = last
        # TODO: Some computation looks wasteful. do we know faster ways to compute these powers?
        # Analyze what these powers look like to find out if there's a better way
        self.t = t
        self.h = h
        self.v = v
        self.a = a
        self.b = b

    def extract_bit(self, exponent: int, i: int) -> int:
        # assert 0 <= i <= self.t, f'i={i} out of range: [0, t={self.t}]'
        # Out of range may naturally happen because of imperfect division of t + 1 by h
        assert 0 <= i
        return (exponent >> i) & 1

    def extract_ea_bit(self, exponent: int, col: int, row: int) -> int:
        assert 0 <= col < self.a, f'col={col} out of range: [0, a={self.a})'
        assert 0 <= row < self.h, f'row={row} out of range: [0, h={self.h})'
        # column 0 is the rightmost column, corresponding to the least significant bit
        return self.extract_bit(exponent, row * self.a + col)
    
    def extract_I_jk(self, exponent: int, j: int, k: int) -> int:
        assert 0 <= j < self.v, f'j={j} out of range: [0, v={self.v})'
        assert 0 <= k < self.b, f'k={k} out of range: [0, b={self.b})'
        col = j * self.b + k
        if col >= self.a:
            return 0 # Out of range may naturally happen because of imperfect division of a by v
        value = 0
        for row in range(self.h):
            value |= self.extract_ea_bit(exponent, col, row) << row
        return value
    
    def compute(self, exponent: int) -> SimulatedElement:
        """Compute the exponentiation using the precomputed powers."""
        assert exponent.bit_length() <= self.t + 1, f'Exponent too large: {exponent.bit_length()} > {self.t + 1}'
        A = self.identity
        for k in range(self.b - 1, -1, -1):
            A = A * A
            for j in range(self.v - 1, -1, -1):
                I_jk = self.extract_I_jk(exponent, j, k)
                if I_jk != 0: # As if self.G[j][0] == identity
                    A = A * self.G[j][I_jk]
        return A

############ FIXED-BASE WINDOW OPTIMIZATION ############
def fbe_analytical_comp_cost(lg_B: int, k: int) -> int:
    """Calculate the analytical cost of fixed-base window exponentiation in number of multiplications. There are no squarings in this algorithm."""
    # t = lg_B / k
    # we multiply t by (2^k - 1) / 2^k to account for no-op when e_i = 0
    # and add 2^k - 2 for the outer A loop
    return math.ceil(lg_B / k) * (2 ** k - 1) / (2 ** k) + 2 ** k - 2

def fbe_optimal_k(lg_B: int) -> int:
    k = 1
    prevcost = float('inf')
    while True:
        cost = fbe_analytical_comp_cost(lg_B, k)
        if cost >= prevcost:
            break
        prevcost = cost
        k += 1
    return k - 1  # Return the last k that was better

def fbe_threshold_table(lg_B_max: int = 100000) -> list[int]:
    """Generate an array containing thresholds of lg(B) where the optimal k changes."""
    threshold_table = []
    prev_k = 1
    for lg_B in range(1, lg_B_max + 1):
        new_k = fbe_optimal_k(lg_B)
        if new_k > prev_k:
            assert new_k - prev_k == 1, f'Expected k to increase by 1, got {new_k - prev_k}'
            threshold_table.append(lg_B)
            prev_k = new_k
    return threshold_table
############ FIXED-BASE WINDOW OPTIMIZATION ############

def msm_naive(lg_es: list[int], p: int):
    """Naive algorithm for multi-exponentiation."""
    metric = PerformanceMetric(num_squares=0, num_mults=0)
    exponents = [random.randrange(2 ** lg_es[i]) for i in range(p)]
    identity = SimulatedVectorElement([0] * p, metric)
    basis = [SimulatedVectorElement([0] * i + [1] + [0] * (p-i-1), metric) for i in range(p)]
    # Naive approach
    current = identity
    for i in range(p):
        current = current * _brauer(exponents[i], lg_es[i], basis[i], identity)
    assert current.exponent == exponents, f'expected {exponents}, got {current.exponent}'
    return metric

def __simple_straus(exponents: list[int], lg_es: list[int], k: list[int], basis: list[SimulatedVectorElement], identity: SimulatedVectorElement) -> SimulatedVectorElement:
    print('k:', k)
    lg_B = max(lg_es)
    # Precompute all p 2^k-size windows
    precomputed_caches = []
    for i, base in enumerate(basis):
        current = base
        precomputed_cache = [base]
        for i in range(2, 2 ** k[i]):
            if i % 2 == 0:
                current = precomputed_cache[i // 2 - 1] * precomputed_cache[i // 2 - 1]
            else:
                current = current * base
            # assert current.exponent == i, "Exponentiation error in Brauer's algorithm"
            precomputed_cache.append(current)
        precomputed_caches.append(precomputed_cache)
        del precomputed_cache
    # Now exponentiate
    bits_lsb_first = [[int(x) for x in bin(exponent)[2:].zfill(lg_B)][::-1] for exponent in exponents]
    current = identity
    for window in range(lg_B - 1, -1, -1):
        for i in range(len(basis)):
            if window % k[i] != 0:
                continue
            bits = bits_lsb_first[i][window:window + k[i]]
            # Convert bits to integer
            window_value = sum(bit * (2 ** i) for i, bit in enumerate(bits))
            # print(window_value, bits)
            if window_value != 0:
                current = current * precomputed_caches[i][window_value - 1]
        if window != 0:
            current = current * current  # Square for the next window
    return current

def _simple_straus(exponents: list[int], lg_es: list[int], basis: list[SimulatedVectorElement], identity: SimulatedVectorElement) -> SimulatedVectorElement:
    return __simple_straus(exponents, lg_es, [brauer_k(lg_e) for lg_e in lg_es], basis, identity)

def msm_straus(lg_es: list[int], p: int):
    """Straus's algorithm for multi-exponentiation."""
    metric = PerformanceMetric(num_squares=0, num_mults=0)
    exponents = [random.randrange(2 ** lg_es[i]) for i in range(p)]
    identity = SimulatedVectorElement([0] * p, metric)
    basis = [SimulatedVectorElement([0] * i + [1] + [0] * (p-i-1), metric) for i in range(p)]
    result = _simple_straus(exponents, lg_es, basis, identity)
    assert result.exponent == exponents, f'expected {exponents}, got {result.exponent}'
    return metric

# Obtained via benchmarking
MUL_COST = 1751e-6 # In ms
SQR_COST = 1220e-6
MOD_COST = 3000e-6
concrete_cost = ConcreteCostModel(cost_per_square=SQR_COST + MOD_COST, cost_per_mult=MUL_COST + MOD_COST)

# model = ConcreteCostModel(cost_per_square=SQR_COST, cost_per_mult=MUL_COST)
model_uniform = ConcreteCostModel(cost_per_square=1, cost_per_mult=1)

# BASELINE_COST_MULTIPLIER = .93
# # BASELINE_COST_MULTIPLIER = 1
# for lg_B in [128, 256, 899, 3072, 3072 * 3]:
#     baseline = sliding_window(lg_B)
#     baseline_cost = model.cost(baseline) * BASELINE_COST_MULTIPLIER
#     k = fbe_optimal_k(lg_B)
#     print('lg(B): ', lg_B, ' k: ', k, ' baseline cost: ', baseline_cost)
#     for choice in [(1, 4, 1), (1, 3, 2), (2, 2, 2)]:
#         print('Choice: ', choice)
#         total_cost = 0
#         for num_invocations in choice:
#             if num_invocations == 1:
#                 total_cost += baseline_cost
#             else:
#                 metric = PerformanceMetric(num_squares=0, num_mults=0)
#                 base = SimulatedElement(1, metric)
#                 identity = SimulatedElement(0, metric)
#                 fixed_base_exp = FBEWindow(identity, base, lg_B, k)
#                 precomp_cost = model.cost(metric)
#                 metric.reset()
#                 expected_exponent = random.randrange(2 ** lg_B)
#                 result = fixed_base_exp.compute(expected_exponent)
#                 comp_cost = model.cost(metric)
#                 total_cost += precomp_cost + comp_cost
#         print('Total cost: ', total_cost)

def theory_total_muls(t: int, h: int, a: int, b: int, v: int, m: int) -> float:
    # a = math.ceil((t + 1) / h)
    assert check_params(t, h, a, b, v), 'Parameters are not optimal'
    precomp_muls = (2**h - 1)*(h/2 - 1)
    comp_muls = a - 1
    total_muls = precomp_muls + m * comp_muls
    return total_muls

def theory_total_squares(t: int, h: int, a: int, b: int, v: int, m: int) -> float:
    # a = math.ceil((t + 1) / h)
    # v = math.ceil(a / b)
    assert check_params(t, h, a, b, v), 'Parameters are not optimal'
    precomp_squares = h * a + (2**h - 1)*(v-1)*b
    comp_squares = b - 1
    total_squares = precomp_squares + m * comp_squares
    return total_squares

def theory_total_cost(t: int, h: int, a: int, b: int, v: int, m: int, model: ConcreteCostModel) -> float:
    total_muls = theory_total_muls(t, h, a, b, v, m)
    total_squares = theory_total_squares(t, h, a, b, v, m)
    return total_muls * model.cost_per_mult + total_squares * model.cost_per_square

def theory_total_muls_continuous(t: int, h: int, b: int, m: int) -> float:
    a = ((t + 1) / h)
    precomp_muls = (2**h - 1)*(h/2 - 1)
    comp_muls = a - 1
    total_muls = precomp_muls + m * comp_muls
    return total_muls

def theory_total_squares_continuous(t: int, h: int, b: int, m: int) -> float:
    a = ((t + 1) / h)
    v = (a / b)
    precomp_squares = h * a + (2**h - 1)*(v-1)*b
    comp_squares = b - 1
    total_squares = precomp_squares + m * comp_squares
    return total_squares

def theory_total_cost_continuous(t: int, h: int, b: int, m: int, model: ConcreteCostModel) -> float:
    total_muls = theory_total_muls_continuous(t, h, b, m)
    total_squares = theory_total_squares_continuous(t, h, b, m)
    return total_muls * model.cost_per_mult + total_squares * model.cost_per_square

def special_comb_optimal_h(lg_B: int, m: int, model: ConcreteCostModel, max_h: int = 16) -> int:
    # For the special case where b = 1 and v = a, where we optimize h
    t = lg_B - 1
    h = 1
    argmin_h = 1
    mincost = float('inf')
    while h <= min(max_h, t + 1):
        a = math.ceil((t + 1) / h)
        b = 1
        v = a
        if check_params(t, h, a, b, v, warn=False) == False:
            h += 1
            continue
        cost = theory_total_cost(t, h, a, b, v, m, model)
        # print(f'h={h}, a={a}, b={b}, v={v}, cost={cost}')
        if cost <= mincost:
            mincost = cost
            argmin_h = h
        h += 1
    return argmin_h  # Return the last k that was better


def check_params(t: int, h: int, a: int, b: int, v: int, warn: bool = True) -> bool:
    # Assert valid parameters, return true if parameters are optimal
    assert a * h >= t + 1, f'a * h = {a * h} < t + 1 = {t + 1}'
    assert b * v >= a, f'b * v = {b * v} < a = {a}'
    # Optimal: no parameter can be reduced without violating constraints
    optimal = True
    if (a - 1) * h >= t + 1:
        if warn:
            print(f'[!] Parameter a={a} is not optimal, can be reduced to {math.ceil((t + 1) / h)}')
        optimal = False
    if a * (h - 1) >= t + 1:
        if warn:
            print(f'[!] Parameter h={h} is not optimal, can be reduced to {math.ceil((t + 1) / (a))}')
        optimal = False
    if (b - 1) * v >= a:
        if warn:
            print(f'[!] Parameter b={b} is not optimal, can be reduced to {math.ceil(a / v)}')
        optimal = False
    if b * (v - 1) >= a:
        if warn:
            print(f'[!] Parameter v={v} is not optimal, can be reduced to {math.ceil(a / b)}')
        optimal = False
    return optimal

def list_optimal_divisors(n: int) -> list[tuple[int, int]]:
    """List all optimal (a, b) pairs of positive integers such that a * b >= n, where optimal means that (a - 1) * b < n and a * (b - 1) < n."""
    divisors = []
    for a in range(1, n + 1):
        b = math.ceil(n / a)
        if (a - 1) * b < n and a * (b - 1) < n:
            divisors.append((a, b))
    return divisors

def test_params(t: int, h: int, a: int, b: int, v: int, model: ConcreteCostModel):
    print(f'Testing parameters: t={t}, h={h}, a={a}, b={b}, v={v}')
    check_params(t, h, a, b, v)
    print('Cost model: ', model)
    metric = PerformanceMetric(num_squares=0, num_mults=0)
    fbe_comb = FBEComb(SimulatedElement(0, metric), SimulatedElement(1, metric), t, h, v, a, b)
    precomp_cost = model.cost(metric)
    precomp_cost_detail = (metric.num_mults, metric.num_squares)
    storage_cost = v * (2**h - 1)
    # pprint.pprint(fbe_comb.G)
    # for exp in range(2**(t+1)):
    costs = []
    costs_detailed = []
    for _ in range(100):
        exp = random.randrange(2 ** (t + 1))
        # print('Exponent:', exp)
        metric.reset()
        result = fbe_comb.compute(exp)
        costs.append(model.cost(metric))
        costs_detailed.append((metric.num_mults, metric.num_squares))
        assert result.exponent == exp, f'expected {exp}, got {result.exponent}'
    print('Precomputation cost: ', precomp_cost, f' = {precomp_cost_detail[0]} M + {precomp_cost_detail[1]} S')
    theory_muls = (2**h - 1)*(h/2 - 1)
    theory_squares = h * a + (2**h - 1)*(v-1)*b
    theory_precomp_total = theory_muls + theory_squares # h * a + (2**h - 1)*((h/2 - 1) + (v-1)*b)
    print('Theoretical precomputation cost: ', theory_precomp_total, f' = {theory_muls} M + {theory_squares} S')
    print('Storage (# elements): ', storage_cost)
    print()
    print(f'Computation cost: min: {min(costs)}, mean: {sum(costs) / len(costs)}, max: {max(costs)}')
    print('Number of multiplications and squarings (min, mean, max):')
    num_mults = [c[0] for c in costs_detailed]
    num_squares = [c[1] for c in costs_detailed]
    print(f'Mults: min: {min(num_mults)}, mean: {sum(num_mults) / len(num_mults)}, max: {max(num_mults)}')
    print(f'Squares: min: {min(num_squares)}, mean: {sum(num_squares) / len(num_squares)}, max: {max(num_squares)}')
    print(f'Theoretical computation cost: {a + b - 2} = {a - 1} M + {b - 1} S')
    print('All tests passed!')
    # total_cost_theoretical = theory_precomp_total + (a + b - 2) * m
    # print('Total theoretical cost for', m, 'exponentiations: ', total_cost_theoretical)

def test_optimization():
    # TODO: investigate deviation from ideal due to rounding (floor/ceil), especially in the b dimension.
    # Is the cost stable across small changes in t (like not a function of t's factorization or other non relevant properties)? expect monotonic increase
    # If not, Can we do some "over-provisioning", like set a = b * v + "small number" to boost performance?
    uniform_cost = True
    if uniform_cost:
        model = model_uniform
    else:
        model = concrete_cost
    t = 50
    # Threshold for uniform cost at t = 50
    m = 16 # Number of exponentiations per base
    # m = 17 # Number of exponentiations per base
    print(f't={t}, m={m}, cost model: {model}')
    optimal_h_float = math.log2(m + 1)
    print('Optimal h (float): ', optimal_h_float)
    # Define the ranges for h and b
    # min_h = math.floor(optimal_h_float)
    # max_h = math.ceil(optimal_h_float)
    max_h = t + 1
    max_b = t + 1
    max_h_graph = 7  # Limit the number of h values shown in the graph for clarity
    max_b_graph = t + 1  # Limit the number of b values shown in the graph for clarity

    # Build a matrix of total costs for h in [1, max_h] and b in [1, max_b]
    argmin = (-1, -1)
    min_cost = float('inf')
    costs = [[float('nan') for _ in range(max_b + 1)] for _ in range(max_h + 1)]
    costs_continuous = [[float('nan') for _ in range(max_b + 1)] for _ in range(max_h + 1)]
    legal_h = list(h for h, _ in list_optimal_divisors(t + 1))
    # legal_b = list(b for _, b in list_optimal_divisors(t + 1))[::-1]
    all_b = list(range(1, t + 2))
    # print('h\\b', end='\t')
    # for b in all_b:
    #     print(b, end='\t')
    # print()
    for h, a in list_optimal_divisors(t + 1):
        # print(h, end='\t')
        for b, v in list_optimal_divisors(a):
            total_cost = theory_total_cost(t, h, a, b, v, m, model)
            costs[h][b] = total_cost
            total_cost_continuous = theory_total_cost_continuous(t, h, b, m, model)
            costs_continuous[h][b] = total_cost_continuous
            # print(f'{total_cost:.2f}', end='\t')
            if total_cost < min_cost:
                min_cost = total_cost
                argmin = (h, b)
        # print()
    print('h\\b', end='\t')
    for b in all_b:
        print(b, end='\t')
    print()
    for h, a in list_optimal_divisors(t + 1):
        print(h, end='\t')
        for b in all_b:
            total_cost = costs[h][b]
            if math.isnan(total_cost):
                print('-', end='\t')
            else:
                print(f'{total_cost:.2f}', end='\t')
            if total_cost < min_cost:
                min_cost = total_cost
                argmin = (h, b)
        print()
    
    print(f'Minimum cost: {min_cost:.2f} at h={argmin[0]}, b={argmin[1]}')
    print(f'Corresponding a={math.ceil((t + 1) / argmin[0])}, v={math.ceil(math.ceil((t + 1) / argmin[0]) / argmin[1])}')
    # Plot the graph of b versus cost for each h
    plt.figure(figsize=(10, 6))
    for h, a in list_optimal_divisors(t + 1):
        if h >= max_h_graph:
            continue
        legal_b = list(b for _, b in list_optimal_divisors(a))[::-1]
        legal_b = [b for b in legal_b if b <= max_b_graph]
        # print(f'legal_b for h={h}, a={a}: {legal_b}')
        # print([costs[h][b] for b, v in list_optimal_divisors(a)])
        assert all([not math.isnan(costs[h][b]) for b in legal_b]), f'Costs not fully populated for h={h}'
        plt.plot(legal_b, [costs[h][b] for b in legal_b], label=f'h={h}, discrete')
        plt.plot(legal_b, [costs_continuous[h][b] for b in legal_b], linestyle='dashed', label=f'h={h}, continuous')
    plt.xlabel('b (number of blocks)')
    plt.ylabel('Total Cost')
    plt.title(f'Total Cost vs b for different h (t={t}, m={m})')
    plt.legend()
    plt.grid(True)
    # plt.yscale('log')  # Use logarithmic scale for better visibility
    plt.show()

if __name__ == "__main__":
    pass
    # print([fbe_optimal_k(lg_B) for lg_B in range(1, 100)])
    # print(fbe_threshold_table())
    # print([special_comb_optimal_h(128, m, model_uniform) for m in range(1, 1000)])
    # print([special_comb_optimal_h(900, m, model_uniform) for m in range(1, 1000)])
    # print([special_comb_optimal_h(900, m, concrete_cost) for m in range(1, 1000)])
    # print([special_comb_optimal_h(9000, m, model_uniform) for m in range(1, 1000)])
    # print([special_comb_optimal_h(900, m, model_uniform) for m in range(1, 100)] == [special_comb_optimal_h(9000, m, model_uniform) for m in range(1, 100)])
    # print([special_comb_optimal_h(lg_B, 1000, model_uniform) for lg_B in range(1, 100)])
    # test_optimization()
    # # A set of optimal parameters
    # # t = 902
    # # h = 7
    # # a = 129
    # # v = 13
    # # b = 10

    # t = 50
    # h = 9
    # a = 6
    # b = 2
    # v = 3

    # # Suboptimal parameter examples
    # # t = 50
    # # h = 10
    # # a = math.ceil((t + 1) / h)
    # # v = 13
    # # b = math.ceil(a / v)

    # # t = 50
    # # a = 10
    # # h = math.ceil((t + 1) / a)
    # # b = 13
    # # v = math.ceil(a / b)

    # # test_params(t, h, a, b, v, model_uniform)

    # for n in range(7,20):
    #     print(f'\n=== Test {n} ===\n')
    #     print(list_optimal_divisors(n))