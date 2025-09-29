import os
import glob
import json
from typing import Literal

os.chdir(os.path.dirname(__file__))


def get_records_location() -> dict[tuple[str, int, int], dict]:
    records: dict[tuple[str, int, int], dict] = {}
    for file in glob.glob(f'../evaluate_rns/anike/*_location_*_stdout.txt'):
        with open(file, 'r') as f:
            lines = f.readlines()
            filtered_lines = [line for line in lines if not line.startswith('#')]
            record = json.loads(''.join(filtered_lines))  # Ensure the first line is valid JSON
        # print(record)
        assert record['anike_app'] == 'geolocation', f'Error: Record is not for geolocation'
        scheme = file[file.index('../evaluate_rns/anike/') + len('../evaluate_rns/anike/'):file.index('_location_')]
        assert scheme in ['flint', 'rns'], f'Error: Invalid acceleration type {scheme}'
        assert record['mkhss_scheme'] == "optimized", f'Error: Invalid MKHSS scheme {record["mkhss_scheme"]}'
        D = record['D']
        L = record['L']
        assert (scheme, D, L) not in records, f'Error: Duplicate record for acceleration type {scheme}, D {D}, L {L}'
        records[(scheme, D, L)] = record
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

# D = 4
# L = 64
# table_cells = [[''] * 3 for _ in range(6)]
table_cells = [[''] * 9 for _ in range(6)]
totals = {}
records = get_records_location()
for i, (D, L) in enumerate([(2, 32), (3, 48), (4, 64)]):
    x_offset = i * 3
    for (column, scheme) in [(x_offset, 'rns'), (x_offset + 1, 'flint')]:
        table_cells[0][column] = display_time_table(records[(scheme, D, L)]['profile']['setup_time'])
        table_cells[1][column] = display_time_table(records[(scheme, D, L)]['profile']['keygenA_time'])
        table_cells[2][column] = display_time_table(records[(scheme, D, L)]['profile']['keygenB_time'])
        table_cells[3][column] = display_time_table(records[(scheme, D, L)]['profile']['keyderA_time'])
        table_cells[4][column] = display_time_table(records[(scheme, D, L)]['profile']['keyderB_time'])
        # Total column
        totals[scheme] = (
            max(
                records[(scheme, D, L)]['profile']['keygenA_time'],
                records[(scheme, D, L)]['profile']['keygenB_time']
            ) +
            max(
                records[(scheme, D, L)]['profile']['keyderA_time'],
                records[(scheme, D, L)]['profile']['keyderB_time']
            )
        )
        table_cells[5][column] = display_time_table(
            totals[scheme]
        )
    # Speedup column
    table_cells[0][x_offset + 2] = rf"${format_sigfig(records[('flint', D, L)]['profile']['setup_time'] / records[('rns', D, L)]['profile']['setup_time'], sig=2)} \times$"
    table_cells[1][x_offset + 2] = rf"${format_sigfig(records[('flint', D, L)]['profile']['keygenA_time'] / records[('rns', D, L)]['profile']['keygenA_time'], sig=2)} \times$"
    table_cells[2][x_offset + 2] = rf"${format_sigfig(records[('flint', D, L)]['profile']['keygenB_time'] / records[('rns', D, L)]['profile']['keygenB_time'], sig=2)} \times$"
    table_cells[3][x_offset + 2] = rf"${format_sigfig(records[('flint', D, L)]['profile']['keyderA_time'] / records[('rns', D, L)]['profile']['keyderA_time'], sig=2)} \times$"
    table_cells[4][x_offset + 2] = rf"${format_sigfig(records[('flint', D, L)]['profile']['keyderB_time'] / records[('rns', D, L)]['profile']['keyderB_time'], sig=2)} \times$"
    table_cells[5][x_offset + 2] = rf"${totals['flint'] / totals['rns'] :.2g} \times$"

end = True
for row in table_cells:
    print(' & '.join(f'{cell}' for cell in row) + (' \\\\' if end else '&'))
