# First, parse result.json, the output of Google Benchmark Tool

# Then, generate LaTeX-formatted tables for the results

import json
import os
import pprint

os.chdir(os.path.dirname(__file__))

# BENCHMARK_UNIT = 'ns'
# MS_PER_BENCHMARK_UNIT = 1e-6
UNITS = {
    'ns': 1e-9,
    'us': 1e-6,
    'ms': 1e-3,
    's': 1,
}

def parse_summary_stats(result_filename: str):
    """
    Parse the result.json file and for each benchmark,
    extract summary statistics, like mean, median, and standard deviation
    """
    with open(result_filename, 'r') as f:
        data = json.load(f)

    # Extract the mean and standard deviation values
    summary_stats = {}
    for benchmark in data['benchmarks']:
        name = benchmark['run_name']
        if name not in summary_stats:
            summary_stats[name] = {}
        if benchmark['run_type'] == 'aggregate':
            time = benchmark['real_time']
            if benchmark['aggregate_unit'] == 'time':
                # assert benchmark['time_unit'] == BENCHMARK_UNIT
                # time = time * MS_PER_BENCHMARK_UNIT
                time = UNITS[benchmark['time_unit']] * time
            summary_stats[name][benchmark['aggregate_name']] = time
        elif benchmark['run_type'] == 'iteration':
            time = benchmark['real_time']
            time = UNITS[benchmark['time_unit']] * time
            summary_stats[name]['mean'] = time
            summary_stats[name]['stddev'] = 0

    return summary_stats

table_columns = { # Map from namespace to column index
    'mkhss': 1,
    'baseline_mkhss': 2,
}

def parse_benchmark_name(benchmark_name: str):
    """
    Parse the benchmark name to extract the procedure name and log_2(B)
    """
    # Example benchmark name: "mkhss::ProcedureName/1"
    # Split by "::" and "/"
    namespace = benchmark_name[:benchmark_name.find('::')]
    procedure = benchmark_name[benchmark_name.rfind('::') + 2:benchmark_name.find('/')]
    log_b = int(benchmark_name[benchmark_name.find('/') + 1:])
    return namespace, procedure, log_b

def display_time(time_in_seconds: float):
    """
    Display the time in seconds in a human-readable format
    """
    if time_in_seconds < 1e-6:
        return f'{time_in_seconds * 1e9:.3f} ns'
    elif time_in_seconds < 1e-3:
        return f'{time_in_seconds * 1e6:.3f} us'
    elif time_in_seconds < 1:
        return f'{time_in_seconds * 1e3:.3f} ms'
    else:
        return f'{time_in_seconds:.3f} s'


def display_time_table(time_in_seconds: float):
    """
    Display the time in seconds in a human-readable format
    """
    if time_in_seconds < 1e-6:
        return f'{time_in_seconds * 1e9:.1f} ns'
    elif time_in_seconds < 1e-3:
        return f'{time_in_seconds * 1e6:.1f} $\\mu$s'
    elif time_in_seconds < 1:
        return f'{time_in_seconds * 1e3:.1f} ms'
    else:
        return f'{time_in_seconds:.1f} s'

def generate_latency_table(summary_stats, want_log_b = 1):
    # First print the data to the console to do a sanity check
    print('Benchmark Name'.ljust(40) + 'Confidence Interval')
    for benchmark, stats in summary_stats.items():
        mean = stats['mean']
        stddev = stats['stddev']
        print(f"{benchmark :40}[{display_time(mean - 2 * stddev)}, {display_time(mean + 2 * stddev)}]")
    ####

    procedure_names: list[str] = []
    for benchmark_name, stats in summary_stats.items():
        benchmark_name: str
        namespace, procedure, log_b = parse_benchmark_name(benchmark_name)
        if procedure not in procedure_names:
            procedure_names.append(procedure)
        # print(namespace, procedure, log_b)
    
    table_cells: list[list[str]] = [
        # Contents: [Procedure name, Ours, Baseline, Our Speedup]
        [f'$\\Method{procedure}{"~*" if procedure == "Setup" else ""}$', '?', '?', '?'] for procedure in procedure_names
    ]
    table_cells_int: list[list[int]] = [
        # Contents: [Procedure name, Ours, Baseline, Our Speedup]
        [0, 0, 0, 0] for _ in procedure_names
    ]

    for benchmark_name, stats in summary_stats.items():
        benchmark_name: str
        namespace, procedure, log_b = parse_benchmark_name(benchmark_name)
        if log_b == want_log_b:
        # table_cells[procedure_names.index(procedure)][table_columns[namespace]] = f'{stats["mean"] :.1f} ms'
            table_cells[procedure_names.index(procedure)][table_columns[namespace]] = display_time_table(stats["mean"])
            table_cells_int[procedure_names.index(procedure)][table_columns[namespace]] = stats["mean"]
        # else:
        #     print(f'Skipping benchmark with log_b = {log_b} != 1:', benchmark_name)

    for i, row in enumerate(table_cells):
        # Convert the last column to speedup
        if row[1] != '?' and row[2] != '?':
            ours = table_cells_int[i][1]
            baseline = table_cells_int[i][2]
            speedup = baseline / ours
            table_cells[i][3] = f'${speedup:.1f} \\times$' if speedup < 10 else f'${speedup:.0f} \\times$'
        else:
            table_cells[i][3] = '?'
    
    pprint.pprint(table_cells)

    # Print the table in LaTeX format (just the rows for now)
    print('% For log_2(B) =', want_log_b)
    for row in table_cells:
        # Print the row in LaTeX format
        print(' & '.join(row) + r' \\')

summary_stats = parse_summary_stats('result.json')
pprint.pprint(summary_stats)
generate_latency_table(summary_stats, want_log_b=1)
