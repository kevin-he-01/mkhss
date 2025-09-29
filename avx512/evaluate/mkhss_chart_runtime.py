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

template = r'''
\begin{tikzpicture}
    \begin{axis}[
        width=\plotwidth,
        height=\plotheight,
        ybar,
        bar width=\barwidth,
        enlargelimits=\enlargelimits,
        legend style={at={(0.5,1.05)}, anchor=south,legend columns=-1},
        ylabel={\ylabel},
        symbolic x coords={<GROUPS>},
        xtick=data,
        xticklabels=<PLSFIX>,
        nodes near coords align={vertical},
        x tick label style={rotate=45, anchor=east},
    ]

    <ADDPLOTS>

    \legend{<LEGENDS>}
    \end{axis}
\end{tikzpicture}
'''

# template_scheme = r'''
# \addplot+[
#     draw=<SCHEME>,
#     fill=<SCHEME>,
#     nodes near coords,
#     every node near coord/.append style={font=\small, color=<SCHEME>}
# ] coordinates {<GROUP_COORDS>};
# '''

template_scheme = r'''
\addplot+[
    draw=<SCHEME>,
    fill=<SCHEME>
] coordinates {<GROUP_COORDS>};
'''

def format_legends(schemes: list[str]):
    return ', '.join(schemes)

def addplot_scheme(scheme: str, group_coords: list[tuple[str, float]]):
    r"""
    Generate the LaTeX code for an \addplot command with the given scheme and group coordinates.
    """
    coords_str = ' '.join([f'({group}, {value})' for group, value in group_coords])
    return template_scheme.replace('<SCHEME>', scheme).replace('<GROUP_COORDS>', coords_str)

def generate_barplot(summary_stats, groups: list[str], eval_targets: list[tuple[str, int]]):
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
    
    if not groups:
        groups = procedure_names
    
    # table_cells: list[list[str]] = [
    #     # Contents: [Procedure name, Ours, Baseline, Our Speedup]
    #     [f'$\\Method{procedure}{"~*" if procedure == "Setup" else ""}$', '?', '?', '?'] for procedure in procedure_names
    # ]
    # table_cells_int: list[list[int]] = [
    #     # Contents: [Procedure name, Ours, Baseline, Our Speedup]
    #     [0, 0, 0, 0] for _ in procedure_names
    # ]

    # Parallel arrays
    # eval_targets = [
    #     ('mkhss', 1),
    #     # ('mkhss', 512),
    #     ('mkhss', 1024),
    #     ('mkhss', 2048),
    #     # ('baseline_mkhss', 1),
    #     # ('baseline_mkhss', 4096),
    # ]
    colornames = {'mkhss': 'ours', 'baseline_mkhss': 'baseline'}
    legend_names = {'mkhss': 'Ours', 'baseline_mkhss': 'Baseline'}
    # legends = [r'Ours ($\log_2(B) = 1$)',
    #         #    r'Ours ($l_B = 512$)',
    #             r'Ours ($\log_2(B) = 1024$)',
    #             r'Ours ($\log_2(B) = 2048$)',
    #         #    r'Baseline ($1 \leq \log_2(B) \leq 2048$)',
    #         #    r'Baseline ($l_B = 4096$)',
    #            ]
    legends = [rf'{legend_names[scheme_name]} ($\logtwo(B) = {log_b}$)' for (scheme_name, log_b) in eval_targets]
    assert len(eval_targets) == len(legends)

    coords: dict[tuple[str, int], list[tuple[str, float]]] = {}
    for benchmark_name, stats in summary_stats.items():
        benchmark_name: str
        namespace, procedure, log_b = parse_benchmark_name(benchmark_name)
        # print(namespace, procedure, log_b)
        if procedure in groups:
            # table_cells[procedure_names.index(procedure)][table_columns[namespace]] = display_time_table(stats["mean"])
            # table_cells_int[procedure_names.index(procedure)][table_columns[namespace]] = stats["mean"]
            key = (namespace, log_b)
            if key not in eval_targets:
                continue
            if key not in coords:
                coords[key] = []
            coords[key].append((procedure, stats['mean'] / UNITS['ms']))

    output = template.replace('<GROUPS>', ', '.join(groups)).replace('<ADDPLOTS>', ''.join(
        addplot_scheme(colornames[scheme_name] + '_' + str(log_b), coords[(scheme_name, log_b)]) for (scheme_name, log_b) in eval_targets
    )).replace('<LEGENDS>', format_legends(legends))
    print(output)

    # for i, row in enumerate(table_cells):
    #     # Convert the last column to speedup
    #     if row[1] != '?' and row[2] != '?':
    #         ours = table_cells_int[i][1]
    #         baseline = table_cells_int[i][2]
    #         speedup = baseline / ours
    #         table_cells[i][3] = f'${speedup:.1f} \\times$' if speedup < 10 else f'${speedup:.0f} \\times$'
    #     else:
    #         table_cells[i][3] = '?'
    
    # pprint.pprint(table_cells)

    # # Print the table in LaTeX format (just the rows for now)
    # for row in table_cells:
    #     # Print the row in LaTeX format
    #     print(' & '.join(row) + r' \\')

summary_stats = parse_summary_stats('result.json')
pprint.pprint(summary_stats)
eval_targets = [
    ('mkhss', 1),
    ('mkhss', 512),
    ('mkhss', 1024),
    ('mkhss', 2048),
    ('baseline_mkhss', 1),
    # ('baseline_mkhss', 4096),
]
# generate_barplot(summary_stats, groups=[])
# generate_barplot(summary_stats, groups=['Setup', 'Keygen', 'Share', 'RMSBootstrapAlice', 'RMSBootstrapBob', 'SyncShareSelf', 'SyncShareOther', 'RMSIAdd', 'RMSISub', 'RMSMAdd', 'RMSMult', 'RMSConvert'])
generate_barplot(summary_stats, groups=['Setup', 'Keygen', 'Share', 'RMSBootstrapAlice', 'RMSBootstrapBob', 'SyncShareSelf', 'SyncShareOther', 'RMSMult', 'RMSConvert'], eval_targets=eval_targets)
# eval_targets = [
#     ('mkhss', 1),
#     # ('mkhss', 512),
#     ('mkhss', 1024),
#     ('mkhss', 2048),
#     # ('baseline_mkhss', 1),
#     # ('baseline_mkhss', 4096),
# ]
# generate_barplot(summary_stats, groups=['Keygen', 'Share', 'RMSBootstrapAlice', 'RMSBootstrapBob', 'SyncShareSelf', 'SyncShareOther'], eval_targets=eval_targets)
# eval_targets = [
#     ('mkhss', 1),
#     # ('mkhss', 512),
#     ('mkhss', 1024),
#     ('mkhss', 2048),
#     ('baseline_mkhss', 1),
#     # ('baseline_mkhss', 4096),
# ]
# generate_barplot(summary_stats, groups=['RMSIAdd', 'RMSISub', 'RMSMAdd'], eval_targets=eval_targets)
# generate_barplot(summary_stats, groups=['Keygen', 'Share', 'RMSBootstrapAlice', 'RMSBootstrapBob', 'SyncShareSelf', 'SyncShareOther', 'RMSMult', 'RMSConvert'])
