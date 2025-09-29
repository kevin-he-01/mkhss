# Simulate exponentiations via addition chains.
# This allows rapid benchmarking of exponentiation algorithms.

# lg_B: number of bits in the exponent

from dataclasses import dataclass
import math
import random
from typing import Callable


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

class FixedBaseExp:
    precomputed_powers: list[SimulatedElement]

    def __init__(self, identity: SimulatedElement, base: SimulatedElement, lg_B: int, k: int):
        """Precompute powers of the base for fixed-base exponentiation."""
        # Compute g, g^{2^k}, g^{2^{2*k}}, ...
        self.base = base
        self.identity = identity
        self.precomputed_powers = [base]
        self.k = k
        self.lg_B = lg_B
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

def fbe_analytical_comp_cost(lg_B: int, k: int) -> int:
    """Calculate the analytical cost of fixed-base exponentiation."""
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

# msm_naive(256, 10)

# for lg_es in [[128, 899], [3072, 3072 * 3]]:
#     p = len(lg_es)
#     print(f"\n\nBenchmarking multi-exponentiation for exponents of {lg_es} bits:")
#     baseline_concrete_cost = concrete_cost.cost(msm_naive(lg_es, p))
#     minimal_cost = concrete_cost.cost(brauer(max(lg_es)))
#     for algorithm in [msm_naive, msm_straus]:
#         metric = algorithm(lg_es, p)
#         print(f"{algorithm.__name__}:")
#         print_cost(metric)
#         print_concrete_cost(metric, concrete_cost)
#         ccc = concrete_cost.cost(metric)
#         speedup = baseline_concrete_cost / ccc
#         slowdown = ccc / minimal_cost
#         print(f"Speedup over baseline: {speedup:.2f}x")
#         print(f"Slowdown over minimal cost: {slowdown:.2f}x")
#         print()

# for lg_B in [256, 899, 3072, 3072 * 3]:
#     print(f"\n\nBenchmarking exponentiation for lg_B={lg_B} bits:")
#     baseline_concrete_cost = concrete_cost.cost(square_and_multiply_mean(lg_B))
#     for algorithm in [square_and_multiply_mean, square_and_multiply, brauer, sliding_window]:
#         metric = algorithm(lg_B)
#         print(f"{algorithm.__name__} for lg_B={lg_B}:")
#         print_cost(metric)
#         print_concrete_cost(metric, concrete_cost)
#         ccc = concrete_cost.cost(metric)
#         speedup = baseline_concrete_cost / ccc
#         print(f"Speedup over baseline: {speedup:.2f}x")
#         print()

# for lg_B in range(18, 30):
#     print(f'lg_B={lg_B}')
#     print([brauer_analytical_cost(lg_B, k) for k in range(1, 10)])
# print(optimal_brauer_table(cost_function=brauer_analytical_cost_simple))
# print(optimal_brauer_table())
# print(optimal_brauer_table(cost_function=lambda lg_B, k: brauer_concrete_cost(lg_B, k, concrete_cost)))

model = ConcreteCostModel(cost_per_square=SQR_COST, cost_per_mult=MUL_COST)
# model = ConcreteCostModel(cost_per_square=1, cost_per_mult=1)

print('Cost model: ', model)
for lg_B in [50, 128, 256, 899, 3072, 3072 * 3]:
    print(f"\n\nBenchmarking fixed-base exponentiation for lg_B={lg_B} bits:")
    baseline = sliding_window(lg_B)
    baseline_cost = model.cost(baseline)
    print(f'Baseline cost: {baseline_cost}')
    costs = []
    fbe_optimal_k_value = fbe_optimal_k(lg_B)
    for k in range(1, 10):
        metric = PerformanceMetric(num_squares=0, num_mults=0)
        base = SimulatedElement(1, metric)
        identity = SimulatedElement(0, metric)
        fixed_base_exp = FixedBaseExp(identity, base, lg_B, k)
        precomp_cost = model.cost(metric)
        metric.reset()
        expected_exponent = random.randrange(2 ** lg_B)
        result = fixed_base_exp.compute(expected_exponent)
        comp_cost = model.cost(metric)
        analytical_comp_cost = fbe_analytical_comp_cost(lg_B, k)
        print(f'k={k}, precomp cost: {precomp_cost}, comp cost: {comp_cost}, analytical comp cost: {analytical_comp_cost}')
        assert result.exponent == expected_exponent, f'expected {expected_exponent}, got {result.exponent}'
        # print(f"k={k}, TC: {metric.num_mults + metric.num_squares}, ", end='')
        costs.append((precomp_cost, comp_cost))
    minimal_k = min(enumerate(costs), key=lambda x: x[1][1])[0] + 1
    print(f'Minimal k: {minimal_k}, Optimal k: {fbe_optimal_k_value}')
    # print(costs)
    precomp_minimal = costs[minimal_k - 1][0]
    comp_minimal = costs[minimal_k - 1][1]
    for num_invocations in [1, 2, 3, 4, 8, 16, 32, 64, 128, 256]:
        amortized_precomp = precomp_minimal / num_invocations
        amortized_comp = comp_minimal
        amortized_cost = amortized_precomp + amortized_comp
        improvement = (baseline_cost - amortized_cost) / baseline_cost * 100
        print(f'For {num_invocations} invocations: '
              f'Amortized precomputation cost: {amortized_precomp:.2f}, '
            #   f'Amortized computation cost: {amortized_comp:.2f}, '
              f'Total amortized cost: {amortized_cost:.2f}, '
              f'Improvement: {improvement:.2f}%')

print([fbe_optimal_k(lg_B) for lg_B in range(1, 100)])
print(fbe_threshold_table())

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
#                 fixed_base_exp = FixedBaseExp(identity, base, lg_B, k)
#                 precomp_cost = model.cost(metric)
#                 metric.reset()
#                 expected_exponent = random.randrange(2 ** lg_B)
#                 result = fixed_base_exp.compute(expected_exponent)
#                 comp_cost = model.cost(metric)
#                 total_cost += precomp_cost + comp_cost
#         print('Total cost: ', total_cost)
