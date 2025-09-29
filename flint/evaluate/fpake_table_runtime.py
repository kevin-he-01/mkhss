import os
import glob
import json
from typing import Literal

os.chdir(os.path.dirname(__file__))


def get_records_fuzzy_pake() -> dict[tuple, dict]:
    records: dict[tuple, dict] = {}
    for file in glob.glob(f'anike/*_fuzzy_pake_*_stdout.txt'):
        with open(file, 'r') as f:
            lines = f.readlines()
            filtered_lines = [line for line in lines if not line.startswith('#')]
            record = json.loads(''.join(filtered_lines))  # Ensure the first line is valid JSON
        # print(record)
        assert record['anike_app'] == 'fuzzy_pake', f'Error: Record is not for fuzzy_pake'
        scheme = record['mkhss_scheme']
        assert scheme in ['optimized', 'baseline'], f'Error: Invalid MKHSS scheme {scheme}'
        L, W, b, T ,Q = map(int, file.split('_')[3:-1])
        assert (scheme, L, W, b, T ,Q) not in records, f'Error: Duplicate record for scheme {scheme}, L {L}, W {W}, b {b}, T {T}, Q {Q}'
        # record['L'] = L
        # record['W'] = W
        # record['b'] = b
        # record['T'] = T
        # record['Q'] = Q
        records[(scheme, L, W, b, T ,Q)] = record
    return records

def format_sigfig(x, sig=3):
    from math import log10, floor
    
    if x == 0:
        return f"{0:.{sig-1}f}"
    
    # order of magnitude
    digits = sig - int(floor(log10(abs(x)))) - 1
    return f"{x:.{max(digits, 0)}f}"

def display_time_table(time_in_seconds: float):
    """
    Display the time in seconds in a human-readable format
    """
    if time_in_seconds < 1e-6:
        return f'{format_sigfig(time_in_seconds * 1e9, sig=3)} ns'
    elif time_in_seconds < 1e-3:
        return f'{format_sigfig(time_in_seconds * 1e6, sig=3)} $\\mu$s'
    elif time_in_seconds < 1:
        return f'{format_sigfig(time_in_seconds * 1e3, sig=3)} ms'
    elif time_in_seconds < 60:
        return f'{format_sigfig(time_in_seconds, sig=3)} s'
    else:
        return f'{format_sigfig(time_in_seconds, sig=3)} s'

parameters = [(8, 9, 5, 2, 2), (10, 8, 16, 1, 1), (12, 10, 8, 3, 3)]

table_cells = [[''] * 9 for _ in range(6)]
totals = {}
records = get_records_fuzzy_pake()
for i, (L, W, b, Q, T) in enumerate(parameters):
    x_offset = i * 3
    # print(records)
    for (column, scheme) in [(x_offset, 'optimized'), (x_offset + 1, 'baseline')]:
        table_cells[0][column] = display_time_table(records[(scheme, L, W, b, Q, T)]['profile']['setup_time'])
        table_cells[1][column] = display_time_table(records[(scheme, L, W, b, Q, T)]['profile']['keygenA_time'])
        table_cells[2][column] = display_time_table(records[(scheme, L, W, b, Q, T)]['profile']['keygenB_time'])
        table_cells[3][column] = display_time_table(records[(scheme, L, W, b, Q, T)]['profile']['keyderA_time'])
        table_cells[4][column] = display_time_table(records[(scheme, L, W, b, Q, T)]['profile']['keyderB_time'])
        # Total column
        totals[scheme] = (
            max(
                records[(scheme, L, W, b, Q, T)]['profile']['keygenA_time'],
                records[(scheme, L, W, b, Q, T)]['profile']['keygenB_time']
            ) +
            max(
                records[(scheme, L, W, b, Q, T)]['profile']['keyderA_time'],
                records[(scheme, L, W, b, Q, T)]['profile']['keyderB_time']
            )
        )
        table_cells[5][column] = display_time_table(
            totals[scheme]
        )
    # Speedup column
    table_cells[0][x_offset + 2] = rf"${format_sigfig(records[('baseline', L, W, b, Q, T)]['profile']['setup_time'] / records[('optimized', L, W, b, Q, T)]['profile']['setup_time'], sig=2)} \times$"
    table_cells[1][x_offset + 2] = rf"${format_sigfig(records[('baseline', L, W, b, Q, T)]['profile']['keygenA_time'] / records[('optimized', L, W, b, Q, T)]['profile']['keygenA_time'], sig=2)} \times$"
    table_cells[2][x_offset + 2] = rf"${format_sigfig(records[('baseline', L, W, b, Q, T)]['profile']['keygenB_time'] / records[('optimized', L, W, b, Q, T)]['profile']['keygenB_time'], sig=2)} \times$"
    table_cells[3][x_offset + 2] = rf"${format_sigfig(records[('baseline', L, W, b, Q, T)]['profile']['keyderA_time'] / records[('optimized', L, W, b, Q, T)]['profile']['keyderA_time'], sig=2)} \times$"
    table_cells[4][x_offset + 2] = rf"${format_sigfig(records[('baseline', L, W, b, Q, T)]['profile']['keyderB_time'] / records[('optimized', L, W, b, Q, T)]['profile']['keyderB_time'], sig=2)} \times$"
    table_cells[5][x_offset + 2] = rf"${format_sigfig(totals['baseline'] / totals['optimized'], sig=2)} \times$"

end = True
for row in table_cells:
    print(' & '.join(f'{cell}' for cell in row) + (' \\\\' if end else '&'))
