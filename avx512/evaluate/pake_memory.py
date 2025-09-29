import sys
from typing import Literal


NIKE_PUBLIC_KEY_SIZE_BYTES = 33 # NIST P-256 public key size, using the compressed form
NIKE_PRIVATE_KEY_SIZE_BYTES = 32 # Order of the curve is a 256-bit prime
LOG_B = 1 # For our evaluation
SEC_PARAM = 128 # Security parameter
STAT_SEC_PARAM = 128 # Statistical security parameter
N_WIDTH = 3072 # Width of N, in bits

def bit2byte(bit: int) -> int:
    """
    Convert bits to bytes
    """
    return (bit + 7) // 8

class MKHSS:
    stat_sec_param: int
    sec_param: int
    n_width_bits: int
    log2_b: int
    w: int
    sub_y_size: int # Size of the subtractive share <y> in bits, assuming |y| = B - 1, this is `bprime_bits` in the C++ code
    sub_ys_size: int # Size of <y * s> in bits, assuming |y| = B - 1, this is `exponent_bits` in the C++ code

    def __init__(self, stat_sec_param: int, sec_param: int, n_width_bits: int, log2_b: int, w: int):
        self.stat_sec_param = stat_sec_param
        self.sec_param = sec_param
        self.n_width_bits = n_width_bits
        self.log2_b = log2_b
        self.w = w
        self.sub_y_size = log2_b + stat_sec_param
        # log_2(B^3 * 2^{4 * sec_par + 3 * stat_sec_par})
        self.sub_ys_size = 3 * log2_b + 4 * sec_param + 3 * stat_sec_param

    def crs_size(self) -> int:
        """
        Calculate the size of the MKHSS CRS structure, in bytes
        """
        # CRS = (N, g mod N^{w+1}, h mod N^{w+1})
        return bit2byte(self.n_width_bits) + 2 * bit2byte(self.n_width_bits * (self.w + 1))

    def public_key_size(self) -> int:
        """
        Calculate the size of the MKHSS public key, in bytes
        """
        # pk = (f mod N^{w+1}, pe_{nim,A} mod N^{w+1}, pe_{nim,B,0} mod N^{w+1}, pe_{nim,B,1} mod N^{w+1})
        return 4 * bit2byte(self.n_width_bits * (self.w + 1))

    def private_key_size(self) -> int:
        """
        Calculate the size of the MKHSS private key, in bytes
        """
        # sk = (s, t, r_A, r_B)
        # s = t * B' + 1, so it can be derived from t, so need not be stored
        # Using SEI, t, r_A, and r_B are 2 * sec_param bits each
        return 3 * bit2byte(2 * self.sec_param)

    def private_share_size(self, x_size_bits: int = 1) -> int:
        """
        Calculate the size of the MKHSS private share, in bytes
        """
        # [[x]]^sigma_sigma = (r, \ct_0 mod N^{w+1})
        # r is 2 * sec_param bits
        return bit2byte(2 * self.sec_param) + bit2byte(self.n_width_bits * (self.w + 1))

    def public_share_size(self) -> int:
        """
        Calculate the size of the MKHSS public share, in bytes
        """
        # [[x]]^sigma_{1-sigma} = (\ct_0 mod N^{w+1}, \ct_1 mod N^{w+1})
        return 2 * bit2byte(self.n_width_bits * (self.w + 1))
    
    def synced_input_share_size(self) -> int:
        """
        Calculate the size of the MKHSS synchronized input share, in bytes
        """
        # {{x}} = (\ct_0 mod N^{w+1}, \ct_1 mod N^{w+1})
        return 2 * bit2byte(self.n_width_bits * (self.w + 1))

    def memory_share_size(self) -> int:
        """
        Calculate the size of the MKHSS memory share, in bytes
        """
        # This depends on log_2(B).
        # Assume this is a fresh share (not a result of arithmetic on multiple memory shares)
        # And the value held in this share is as large as B
        return bit2byte(self.sub_y_size) + bit2byte(self.sub_ys_size)

class BaselineMKHSS:
    stat_sec_param: int
    sec_param: int
    n_width_bits: int
    log2_b: int
    w: int

    def __init__(self, stat_sec_param: int, sec_param: int, n_width_bits: int, log2_b: int, w: int):
        self.stat_sec_param = stat_sec_param
        self.sec_param = sec_param
        self.n_width_bits = n_width_bits
        self.log2_b = log2_b
        self.w = w

    def crs_size(self) -> int:
        """
        Calculate the size of the MKHSS CRS structure, in bytes
        """
        # CRS = (N, g mod N^{w+1}, h mod N^{w+1})
        return bit2byte(self.n_width_bits) + 2 * bit2byte(self.n_width_bits * (self.w + 1))

    def public_key_size(self) -> int:
        """
        Calculate the size of the MKHSS public key, in bytes
        """
        # pk = (f mod N^{w+1}, pe_{nim,A} mod N^{w+1}, pe_{nim,B,0} mod N^{w+1}, pe_{nim,B,1} mod N^{w+1}, pk_{nike})
        return 4 * bit2byte(self.n_width_bits * (self.w + 1)) + NIKE_PUBLIC_KEY_SIZE_BYTES

    def private_key_size(self) -> int:
        """
        Calculate the size of the MKHSS private key, in bytes
        """
        # sk = (s, r_A, r_B)
        # s = t * N + 1, so it can be derived from t, so need not be stored
        # Using SEI, t, r_A, and r_B are 2 * sec_param bits each
        return 3 * bit2byte(2 * self.sec_param) + NIKE_PRIVATE_KEY_SIZE_BYTES

    def private_share_size(self, x_size_bits: int = 1) -> int:
        """
        Calculate the size of the MKHSS private share, in bytes
        """
        # [[x]]^sigma_sigma = (x, r, r', \ct_0 mod N^{w+1}, \ct_0' mod N^{2})
        # r and r' are 2 * sec_param bits
        # By default, assume best case scenario (x is a bit) to give the baseline an advantage
        # Plus it is more realistic in the ANIKE setting
        return bit2byte(x_size_bits) + 2 * bit2byte(2 * self.sec_param) + bit2byte(self.n_width_bits * (self.w + 1)) + bit2byte(self.n_width_bits * 2)

    def public_share_size(self) -> int:
        """
        Calculate the size of the MKHSS public share, in bytes
        """
        # [[x]]^sigma_{1-sigma} = (\ct_0 mod N^{w+1}, \ct_1 mod N^{w+1}, \ct_0' mod N^{2}, \ct_1' mod N^{2})
        return 2 * bit2byte(self.n_width_bits * (self.w + 1)) + 2 * bit2byte(self.n_width_bits * 2)
    
    def synced_input_share_size(self) -> int:
        """
        Calculate the size of the MKHSS synchronized input share, in bytes
        """
        # {{x}} = (\ct_0 mod N^{w+1}, \ct_1 mod N^{w+1}, \ct_0' mod N^{2}, \ct_1' mod N^{2})
        return 2 * bit2byte(self.n_width_bits * (self.w + 1)) + 2 * bit2byte(self.n_width_bits * 2)

    def memory_share_size(self) -> int:
        """
        Calculate the size of the MKHSS memory share, in bytes
        """
        # <<y*s>> is mod N^w
        # <<y>> is mod N
        return bit2byte(self.n_width_bits) + bit2byte(self.n_width_bits * self.w)

HSS = MKHSS | BaselineMKHSS

# Fuzzy PAKE
class FPAKE:
    stat_sec_param: int
    sec_param: int
    n_width_bits: int
    L: int
    W: int
    b: int
    Q: int
    T: int
    total_length: int
    hss: HSS

    def __init__(self, stat_sec_param: int, sec_param: int, n_width_bits: int, L: int, W: int, b: int, Q: int, T: int, hss: Literal['mkhss', 'baseline_mkhss']):
        self.L = L
        self.W = W
        self.b = b
        self.Q = Q
        self.T = T
        self.total_length = L * W * b
        self.stat_sec_param = stat_sec_param
        self.sec_param = sec_param
        self.n_width_bits = n_width_bits
        log_b = 1
        if hss == 'mkhss':
            self.hss = MKHSS(stat_sec_param, sec_param, n_width_bits, log_b, w = 1)
        elif hss == 'baseline_mkhss':
            self.hss = BaselineMKHSS(stat_sec_param, sec_param, n_width_bits, log_b, w = 3)
        else:
            raise ValueError("Invalid HSS type. Must be 'mkhss' or 'baseline_mkhss'.")

    def anike_crs_size(self) -> int:
        """
        Calculate the size of the ANIKE CRS structure, in bytes
        """
        return self.hss.crs_size()

    def anike_pe_alice_size(self) -> int:
        """
        Calculate the size of the ANIKE public encoding, in bytes
        """
        return self.hss.public_key_size() + self.total_length * self.hss.public_share_size() + NIKE_PUBLIC_KEY_SIZE_BYTES
    
    def anike_pe_bob_size(self) -> int:
        """
        Calculate the size of the ANIKE public encoding, in bytes
        """
        return self.hss.public_key_size() + self.total_length * self.hss.public_share_size() + NIKE_PUBLIC_KEY_SIZE_BYTES

    def anike_st_alice_size(self) -> int:
        """
        Calculate the size of the ANIKE private state, in bytes
        """
        return self.hss.private_key_size() \
                + self.total_length * self.hss.private_share_size(x_size_bits=1) + NIKE_PRIVATE_KEY_SIZE_BYTES
    
    def anike_st_bob_size(self) -> int:
        """
        Calculate the size of the ANIKE private state, in bytes
        """
        return self.hss.private_key_size() \
                + self.total_length * self.hss.private_share_size(x_size_bits=1) + NIKE_PRIVATE_KEY_SIZE_BYTES

table = [[''] * 9 for _ in range(4)]

for i, (L, W, b, Q, T) in enumerate([(8, 9, 5, 2, 2), (10, 8, 16, 1, 1), (12, 10, 8, 3, 3)]):
    x_offset = i * 3

    ours_fpake = FPAKE(
        stat_sec_param=STAT_SEC_PARAM,
        sec_param=SEC_PARAM,
        n_width_bits=N_WIDTH,
        L=L,
        W=W,
        b=b,
        Q=Q,
        T=T,
        hss='mkhss'
    )

    print(f'% FPAKE with L={L}, W={W}, b={b}, Q={Q}, T={T}, total length={ours_fpake.total_length}')
    baseline_fpake = FPAKE(
        stat_sec_param=STAT_SEC_PARAM,
        sec_param=SEC_PARAM,
        n_width_bits=N_WIDTH,
        L=L,
        W=W,
        b=b,
        Q=Q,
        T=T,
        hss='baseline_mkhss'
    )

    table[0][x_offset] = f'{ours_fpake.anike_crs_size() / 1000 :.1f} kB'
    table[0][x_offset + 1] = f'{baseline_fpake.anike_crs_size() / 1000 :.1f} kB'
    table[0][x_offset + 2] = f'${baseline_fpake.anike_crs_size() / ours_fpake.anike_crs_size():.1f} \\times$'
    table[1][x_offset] = f'{ours_fpake.anike_pe_alice_size() / 1000 :.1f} kB'
    table[1][x_offset + 1] = f'{baseline_fpake.anike_pe_alice_size() / 1000 :.1f} kB'
    table[1][x_offset + 2] = f'${baseline_fpake.anike_pe_alice_size() / ours_fpake.anike_pe_alice_size():.1f} \\times$'
    table[2][x_offset] = f'{ours_fpake.anike_pe_bob_size() / 1000 :.1f} kB'
    table[2][x_offset + 1] = f'{baseline_fpake.anike_pe_bob_size() / 1000 :.1f} kB'
    table[2][x_offset + 2] = f'${baseline_fpake.anike_pe_bob_size() / ours_fpake.anike_pe_bob_size():.1f} \\times$'
    # sum of pe + st
    table[3][x_offset] = f'{(ours_fpake.anike_pe_alice_size() + ours_fpake.anike_pe_bob_size()) / 1000 :.1f} kB'
    table[3][x_offset + 1] = f'{(baseline_fpake.anike_pe_alice_size() + baseline_fpake.anike_pe_bob_size()) / 1000 :.1f} kB'
    table[3][x_offset + 2] = f'${(baseline_fpake.anike_pe_alice_size() + baseline_fpake.anike_pe_bob_size()) / (ours_fpake.anike_pe_alice_size() + ours_fpake.anike_pe_bob_size()):.1f} \\times$'
# Print the table in LaTeX format
for row in table:
    print(' & '.join(row) + r' \\')

# geoke_table_template = r'''
# % geolocation-based ANIKE (GeoKE) using MKHSS
# \begin{tabular}{llll}
# \toprule
#     Data           & Ours & Baseline & Our saving \\
#     \hline
# [LINES]
# [TOTAL]
# \bottomrule
# \end{tabular}
# '''[:-1]

# geoke_total_row = r'''
#     \TotalComm* & {ours_total} kB & {baseline_total} kB & ${saving_total} \times$ \\
# '''[:-1]

# lines = []
# for (d, l) in [(4, 64)]:
#     lines.append(f'% Geolocation-based ANIKE (GeoKE) with d={d}, l={l}')
#     our_anike = GeoKE(
#         stat_sec_param=STAT_SEC_PARAM,
#         sec_param=SEC_PARAM,
#         n_width_bits=N_WIDTH,
#         l=l,
#         d=d,
#         hss='mkhss'
#     )

#     baseline_anike = GeoKE(
#         stat_sec_param=STAT_SEC_PARAM,
#         sec_param=SEC_PARAM,
#         n_width_bits=N_WIDTH,
#         l=l,
#         d=d,
#         hss='baseline_mkhss'
#     )
#     # anike_schemes: list[ANIKE] = [our_anike, baseline_anike]
#     # for scheme in anike_schemes:
#     #     print('=' * 40)
#     #     print(f"MKHSS Scheme: {scheme.hss.__class__.__name__}")
#     #     # print(f"Statistical security parameter: {scheme.stat_sec_param} bits")
#     #     # print(f"Security parameter: {scheme.sec_param} bits")
#     #     # print(f"Width of N: {scheme.n_width_bits} bits")
#     #     print(f"L: {scheme.l} bits")
#     #     print('-' * 40)
#     #     print(f"CRS size: {scheme.anike_crs_size()} bytes")
#     #     print(f"Public encoding size: {scheme.anike_pe_size()} bytes")
#     #     print(f"Private state size: {scheme.anike_st_size()} bytes")
#     for data_name in ['crs', 'pe_alice', 'pe_bob', 'st_alice', 'st_bob']:
#         our_size = getattr(our_anike, f"anike_{data_name}_size")()
#         baseline_size = getattr(baseline_anike, f"anike_{data_name}_size")()
#         saving = baseline_size / our_size
#         line = r'''
#     {name} & {ours} kB & {baseline} kB & ${saving} \times$ \\
# '''[:-1].format(
#             name='\\DataName' + ''.join(s.capitalize() for s in data_name.split('_')),
#             ours=f'{our_size / 1000:.1f}',
#             baseline=f'{baseline_size / 1000:.1f}',
#             saving=f'{saving:.1f}'
#         )
#         lines.append(line)
#     # lines.append('\n' + r'    \hline')

# lines = ''.join(lines)
# # Replace [LINES] in template with lines
# geoke_table_template = geoke_table_template.replace('[LINES]', lines).replace('[TOTAL]', geoke_total_row.format(
#     ours_total=f'{(our_anike.anike_pe_alice_size() + our_anike.anike_pe_bob_size()) / 1000:.1f}',
#     baseline_total=f'{(baseline_anike.anike_pe_alice_size() + baseline_anike.anike_pe_bob_size()) / 1000:.1f}',
#     saving_total=f'{(baseline_anike.anike_pe_alice_size() + baseline_anike.anike_pe_bob_size()) / (our_anike.anike_pe_alice_size() + our_anike.anike_pe_bob_size()):.1f}'
# ))
# print(geoke_table_template)
