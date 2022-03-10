import os
import numpy as np
import matplotlib.pyplot as plt


def get_results(results_dir, solve_methods, Gammas):
    """
    Reads specified results files from directory and returns dictionary of results data.

    :param results_dir: Absolute filepath of directory containing results files. Directory must contain results files
        for all combinations of solve methods and values of Gamma.
    :type results_dir: str
    :param solve_methods: List of solution methods to read results for.
    :type solve_methods: list
    :param Gammas: List of Gamma values to read results for.
    :type Gammas: list
    :return: Dictionary of results for each combination of solve method and Gamma value.
    :rtype: dict
    """
    # check results files exist for each combination solve_method and Gamma
    missing_files = []
    results_files = [f for f in os.listdir(results_dir) if os.path.splitext(f)[1] == '.txt']
    for method in solve_methods:
        for g in Gammas:
            required_file = f"{method}_results_{g}.txt"
            if required_file not in results_files:
                missing_files.append(required_file)
    if len(missing_files) > 0:
        print("The following results files have not been found:")
        for v in missing_files:
            print(f"\t{v}")
        raise Exception("Please add the missing files to the specified results directory and try again.")

    # get data from results files
    results_data = {}
    for method in solve_methods:
        results_data[method] = {}
        for g in Gammas:
            file_name = os.path.join(results_dir, f"{method}_results_{g}.txt")
            file_data = {}
            with open(file_name, 'r') as f:
                raw_lines = f.read().splitlines()
                for line in raw_lines[1:]:
                    splitline = line.split()
                    instance = splitline[0]
                    if 'benders' in method:
                        file_data[instance] = {'objval': float(splitline[2]), 'objbound': float(splitline[3]),
                                               'gap': float(splitline[4]), 'n_iterations': int(splitline[5]),
                                               'avg_iteration': float(splitline[6]), 'runtime': float(splitline[7])}
                    elif 'compact_reformulation' in method:
                        file_data[instance] = {'objval': float(splitline[2]), 'objbound': float(splitline[3]),
                                               'gap': float(splitline[4]), 'runtime': float(splitline[5])}
            results_data[method][g] = file_data
    return results_data


def plot_gaps(results_data, solve_methods, Gammas, save_path, legend=True):
    """
    Plots optimality gaps (x-axis) against the proportion of instances solved to within that gap (y-axis) for specified
    solve methods.

    :param results_data: Dictionary of results data to plot. Must contain results for each combination of
        solve method and Gamma specified in the input.
    :type results_data: dict
    :param solve_methods: List of solution methods to plot results for.
    :type solve_methods: list
    :param Gammas: List of Gamma values to use in plot data.
    :type Gammas: list
    :param save_path: Absolute filepath of directory in which to save plot.
    :type save_path: str
    :param legend: Indicates whether or not to add legend to plot. Defaults to True.
    :type legend: bool
    """
    # Check results_data has data for all combinations of solve method and Gamma
    missing_results = []
    for method in solve_methods:
        for g in Gammas:
            try:
                v = results_data[method][g]
            except KeyError:
                missing_results.append((method, g))
    if len(missing_results) > 0:
        print("Results for the following combinations of solve method and Gamma have not been found:")
        for v in missing_results:
            print(f"\tSolve method = '{v[0]}', Gamma = {v[1]}")
        raise Exception("Please correct this and try again.")

    # Get optimality gaps
    instances = [i for i in results_data[solve_methods[0]][Gammas[0]]]  # Number of instances for each Gamma
    gaps = {}
    for method in solve_methods:
        gaps[method] = [results_data[method][g][i]['gap'] for g in Gammas for i in instances]

    # Get proportion of instances with optimality gap smaller than value
    n_total_instances = len(gaps[solve_methods[0]])  # Total number of instances across all Gamma values
    max_gap = max(max(gaps[method]) for method in solve_methods)
    n_steps = 100
    gap_values = np.linspace(0, max_gap, n_steps)
    less_than_gap = {}
    for method in solve_methods:
        less_than_gap[method] = [len([i for i in range(n_total_instances) if gaps[method][i] <= v]) / n_total_instances
                                 for v in gap_values]

    # Plot figure
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})  # use latex fonts
    plt.rc('text', usetex=True)
    plt.figure()
    for method in solve_methods:
        plt.plot([100 * v for v in gap_values], [100 * v for v in less_than_gap[method]],
                 color=get_solve_method_line_colour(method), label=get_solve_method_legend_name(method))
    plt.xlim([0, 100 * max_gap])
    plt.ylim([0, 105])
    plt.xlabel(r'Gap (\%)')
    plt.ylabel(r'\% of instances')
    if legend is True:
        plt.legend(loc=4, frameon=False)
    plt.savefig(os.path.join(save_path, 'gaps.pdf'))


def plot_performance_profile(results_data, solve_methods, Gammas, save_path, legend=True):
    """
    Plots performance profile of solution times for solve methods.

    :param results_data: Dictionary of results data to plot. Must contain results for each combination of
        solve method and Gamma specified in the input.
    :type results_data: dict
    :param solve_methods: List of solution methods to plot results for.
    :type solve_methods: list
    :param Gammas: List of Gamma values to use to in plot data.
    :type Gammas: list
    :param save_path: Absolute filepath of directory in which to save plot.
    :type save_path: str
    :param legend: Indicates whether or not to add legend to plot. Defaults to True.
    :type legend: bool
    """
    # Check results_data has data for all combinations of solve method and Gamma
    missing_results = []
    for method in solve_methods:
        for g in Gammas:
            try:
                v = results_data[method][g]
            except KeyError:
                missing_results.append((method, g))
    if len(missing_results) > 0:
        print("Results for the following combinations of solve method and Gamma have not been found:")
        for v in missing_results:
            print(f"\tSolve method = '{v[0]}', Gamma = {v[1]}")
        raise Exception("Please correct this and try again.")

    # Get solution times
    instances = [i for i in results_data[solve_methods[0]][Gammas[0]]]  # Number of instances for each Gamma
    runtimes = {}
    gaps = {}
    for method in solve_methods:
        runtimes[method] = [results_data[method][g][i]['runtime'] for g in Gammas for i in instances]
        gaps[method] = [results_data[method][g][i]['gap'] for g in Gammas for i in instances]

    # Get performance ratios, i.e. ratio between solution time and best solution time for that instance
    n_total_instances = len(runtimes[solve_methods[0]])  # Total number of instances across all Gamma values
    perf_ratios = {}
    for method in solve_methods:
        perf_ratios[method] = []
        for i in range(n_total_instances):
            if gaps[method][i] == 0:
                perf_ratios[method].append(runtimes[method][i] / min(runtimes[m][i] for m in solve_methods))
            elif (gaps[method][i] > 0) or (np.isnan(gaps[method][i])):
                # add penalty if no optimal solution was found
                perf_ratios[method].append(1e6)

    # Get performance profile, i.e. proportion of instances with performance ratio smaller than value
    perf_profile = {}
    taus = np.linspace(0, 10, 100)
    for method in solve_methods:
        perf_profile[method] = [
            len([i for i in range(n_total_instances) if np.log(perf_ratios[method][i]) <= t]) / n_total_instances for t
            in taus]

    # Plot figure
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})  # use latex fonts
    plt.rc('text', usetex=True)
    plt.figure()
    for method in solve_methods:
        plt.plot(taus, [100 * p for p in perf_profile[method]], color=get_solve_method_line_colour(method),
                 label=get_solve_method_legend_name(method))
    plt.xlim([0, max(taus)])
    plt.ylim([0, 105])
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'\% of instances with $\log(r_{ia})\leq \tau$')
    if legend is True:
        plt.legend(loc=4, frameon=False)
    plt.savefig(os.path.join(save_path, 'performance_profile.pdf'))


def get_solve_method_line_colour(solve_method):
    """
    Returns hexadeximal colour code to use in plots for given solve method.
    """
    solve_method_colours = {'compact_reformulation_trans': '#ff7f0e', 'benders': '#a6bddb', 'new_benders': '#2b8cbe'}
    try:
        colour = solve_method_colours[solve_method]
    except KeyError:
        print(f"{solve_method} is not a solve method (must be either 'compact_reformulation', 'benders' or "
              f"'new_benders'). Please correct this and try again!")
        raise
    return colour


def get_solve_method_legend_name(solve_method):
    """
    Converts solve method name (as used in the code) into a user-readable solve method name to be printed in legend.
    """
    solve_methods = {'compact_reformulation_trans': 'Compact formulation',
                     'benders': "Balouka \& Cohen (2021) Benders'", 'new_benders': "Alg. 1 Benders'"}
    try:
        legend_name = solve_methods[solve_method]
    except KeyError:
        print(f"{solve_method} is not a recognised results file solve method (must be either "
              f"'compact_reformulation', 'benders' or 'new_benders'). Please correct this and try again!")
        raise
    return legend_name


if __name__ == "__main__":
    results_folder = '/home/boldm1/OneDrive/robust-mrcpsp/code/results/storm_j20_7200/'
    methods = ['benders', 'new_benders', 'compact_reformulation_trans']
    Gammas = [5, 10, 15]
    results = get_results(results_folder, methods, Gammas)
    plot_gaps(results, methods, Gammas, '/home/boldm1/OneDrive/robust-mrcpsp/code/plots/')
    plot_performance_profile(results, methods, Gammas, '/home/boldm1/OneDrive/robust-mrcpsp/code/plots/', legend=False)
