import math

TOTAL_WORDS_FULL = 7776

L_FULL = 4

# passphrase/final_filtered_eff_large_wordlist.txt
TOTAL_WORDS_FILTERED = 6403
NUM_NEIGHBORS_MAX = 10

L = 8
Q = 2
W = 9
T = 2

print('Per-word entropy (full):', math.log2(TOTAL_WORDS_FULL))
print('Total entropy (full):', math.log2(TOTAL_WORDS_FULL) * L_FULL)

print('Filtered per-word entropy:', math.log2(TOTAL_WORDS_FILTERED))
print('Filtered total entropy:', math.log2(TOTAL_WORDS_FILTERED) * L)
print('Conservative estimate:', (L-T) * math.log2(TOTAL_WORDS_FILTERED / NUM_NEIGHBORS_MAX) - math.log2(math.comb(L,T)))

def estimate_alt():
    return sum([math.comb(L, k) * (NUM_NEIGHBORS_MAX / TOTAL_WORDS_FILTERED) ** k * (1 - 1 / TOTAL_WORDS_FILTERED) ** (L - k) for k in range(L-T, L+1)])

print('Alternative estimate:', -math.log2(estimate_alt()))
