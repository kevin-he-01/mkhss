# First, parse result.json, the output of Google Benchmark Tool

# Then, generate LaTeX-formatted tables for the results

import json
import os
import sys
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

    return summary_stats

def parse_benchmark_name(benchmark_name: str):
    vals = benchmark_name.split('/')
    name = vals[0]
    type = name[:name.find('::')]
    name = name[name.rfind('::') + 2:]  # Remove namespace prefix
    assert type in ['rns', 'fmpz_bench'], f'Unexpected benchmark type: {type}'
    if type == 'fmpz_bench':
        if len(vals) == 3:
            return name, int(vals[-2]), int(vals[-1]), None
        elif len(vals) == 4:
            return name, int(vals[-3]), int(vals[-1]), int(vals[-2])
    else:
        if len(vals) == 2:
            return name, int(vals[-1]), 1, None
        elif len(vals) == 3:
            return name, int(vals[-2]), 1, int(vals[-1])
    assert False, f'Unexpected benchmark name format: {benchmark_name} for type {type}, {len(vals)} parts'

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

def visualize(summary_stats):
    # First print the data to the console to do a sanity check
    print('Benchmark Name'.ljust(60) + 'Confidence Interval')
    for benchmark, stats in summary_stats.items():
        mean = stats['mean']
        stddev = stats['stddev']
        print(f"{benchmark :60}[{display_time(mean - 2 * stddev)}, {display_time(mean + 2 * stddev)}]")
    ####

    data: dict[tuple[int, int], dict] = {}
    for benchmark_name, stats in summary_stats.items():
        # benchmark_name: str
        name, width, w, window_size = parse_benchmark_name(benchmark_name)
        print(width, w, window_size, name)
        if (width, w) not in data:
            data[(width, w)] = {}
        if name == 'ModExpFixedBasePrecompute':
            if 'precompute' not in data[(width, w)]:
                data[(width, w)]['precompute'] = dict()
            data[(width, w)]['precompute'][window_size] = stats['mean']
        elif name == 'ModExpFixedBase':
            if 'fixed_base' not in data[(width, w)]:
                data[(width, w)]['fixed_base'] = dict()
            data[(width, w)]['fixed_base'][window_size] = stats['mean']
        elif name == 'ModExpFmpz':
            data[(width, w)]['baseline'] = stats['mean']
        # print(namespace, procedure, log_b)
    
    pprint.pprint(data)

    # for (width, w), stats in data.items():
    #     print(f'Width: {width}, w: {w}')
    #     for window_size in range(1, 9):
    #         print(f'  Window size: {window_size}')
    #         for num_invocations in [4, 8, 16, 32, 64]:
    #             amortized_time = stats['precompute'][window_size] / num_invocations + stats['fixed_base'][window_size]
    #             improvement = (stats['baseline'] - amortized_time) / stats['baseline'] * 100
    #             print(f'    Num invocations: {num_invocations}, ' +
    #                   f'Amortized time: {display_time(amortized_time)}, Improvement: {improvement:.2f}%')

    # for (width, w), stats in data.items():
    #     print(f'Width: {width}, w: {w}')
    #     for num_invocations in [4, 8, 16, 32, 64]:
    #         print(f'  Num invocations: {num_invocations}, ')
    #         for window_size in range(1, 9):
    #             amortized_time = stats['precompute'][window_size] / num_invocations + stats['fixed_base'][window_size]
    #             improvement = (stats['baseline'] - amortized_time) / stats['baseline'] * 100
    #             print(f'    Window size: {window_size}, ' +
    #                   f'Amortized time: {display_time(amortized_time)}, Improvement: {improvement:.2f}%')

    # for (width, w), stats in data.items():
    #     print(f'Width: {width}, w: {w}')
    #     for num_invocations in [4, 8, 16, 32, 64, 128, 256]:
    #         print(f'  Num invocations: {num_invocations}, ')
    #         min_amortized_time = float('inf')
    #         argmin = None
    #         for window_size in range(1, 9):
    #             amortized_time = stats['precompute'][window_size] / num_invocations + stats['fixed_base'][window_size]
    #             if amortized_time < min_amortized_time:
    #                 min_amortized_time = amortized_time
    #                 argmin = window_size
    #         improvement = (stats['baseline'] - min_amortized_time) / stats['baseline'] * 100
    #         print(f'    Optimal window size: {argmin}, ' +
    #                 f'Amortized time: {display_time(min_amortized_time)}, Improvement: {improvement:.2f}%')
    
    precomputed_tune_table = ''
    for (width, w), stats in data.items():
        optimality_table = []
        print(f'Width: {width}, w: {w}')
        range_max = 4096
        optimal_window_sizes = [0] * (range_max + 1)
        improvements = [0] * (range_max + 1)
        for num_invocations in range(2, range_max + 1):
            # print(f'  Num invocations: {num_invocations}, ')
            min_amortized_time = float('inf')
            argmin = None
            for window_size in list(range(1, 9)) + [-1]:
                amortized_time = stats['precompute'][window_size] / num_invocations + stats['fixed_base'][window_size]
                if amortized_time < min_amortized_time:
                    min_amortized_time = amortized_time
                    argmin = window_size
            assert argmin is not None
            optimal_window_sizes[num_invocations] = argmin
            improvements[num_invocations] = (stats['baseline'] - min_amortized_time) / stats['baseline'] * 100
            # improvement = (stats['baseline'] - min_amortized_time) / stats['baseline'] * 100
            # print(f'    Optimal window size: {argmin}, ' +
            #         f'Amortized time: {display_time(min_amortized_time)}, Improvement: {improvement:.2f}%')
        # print(f'  Optimal window sizes: {optimal_window_sizes[3:]}')
        # Print out the first indices where the optimal window size changes
        last_window_size = optimal_window_sizes[0]
        for num_invocations in range(1, range_max + 1):
            if optimal_window_sizes[num_invocations] != last_window_size:
                print(f'    For {num_invocations} invocations, optimal window size is {optimal_window_sizes[num_invocations]}, ' +
                      f'with improvement of {improvements[num_invocations]:.2f}%')
                optimality_table.append((num_invocations, optimal_window_sizes[num_invocations]))
                last_window_size = optimal_window_sizes[num_invocations]
        optimality_table = optimality_table[::-1] # To help C code
        precomputed_tune_table += ('{' + str(width) + ', {' + ', '.join([f'{{{x[0]}, {x[1]}}}' for x in optimality_table]) + '}},\n')
    print('Precomputed tune table:')
    print(precomputed_tune_table)

if len(sys.argv) != 2:
    print(f'Usage: {sys.argv[0]} <[rns/flint]_bench_out.json>')
    sys.exit(1)
filename = str(sys.argv[1])
summary_stats = parse_summary_stats(filename)
# pprint.pprint(summary_stats)
visualize(summary_stats)
