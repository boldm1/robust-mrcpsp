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
                    if method == 'benders':
                        file_data[instance] = {'objval': float(splitline[2]), 'objbound': float(splitline[3]),
                                               'gap': float(splitline[4]), 'n_iterations': int(splitline[5]),
                                               'avg_iteration': float(splitline[6]), 'runtime': float(splitline[7])}
                    elif method == 'compact_reformulation':
                        file_data[instance] = {'objval': float(splitline[2]), 'objbound': float(splitline[3]),
                                               'gap': float(splitline[4]), 'runtime': float(splitline[5])}
            results_data[method][g] = file_data
    return results_data


def plot_gaps(results_data, solve_methods, Gammas, save_path):
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
        method_name = method.replace('_', ' ').capitalize()  # format method name for legend
        plt.plot([100 * v for v in gap_values], [100 * v for v in less_than_gap[method]], label=method_name)
    plt.xlim([0, 100 * max_gap])
    plt.xlabel(r'Gap (\%)')
    plt.ylabel(r'\% of instances')
    plt.legend(loc=4, frameon=False)
    plt.savefig(os.path.join(save_path, 'gaps.pdf'))


def plot_performance_profile(results_data, solve_methods, Gammas, save_path):
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
        method_name = method.replace('_', ' ').capitalize()  # format method name for legend
        plt.plot(taus, perf_profile[method], label=method_name)
    plt.xlim([0, max(taus)])
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$\mathrm{P}(\log(p_{im})\leq \tau)$')
    plt.legend(loc=4, frameon=False)
    plt.savefig(os.path.join(save_path, 'performance_profile.pdf'))


if __name__ == "__main__":
    results = get_results('/home/boldm1/OneDrive/robust-mrcpsp/code/results/storm_j10_7200/',
                           ['benders', 'compact_reformulation'], [3, 5, 7])
    plot_gaps(results, ['benders', 'compact_reformulation'], [3, 5, 7],
              '/home/boldm1/OneDrive/robust-mrcpsp/code/plots/')
    plot_performance_profile(results, ['benders', 'compact_reformulation'], [3, 5, 7],
                             '/home/boldm1/OneDrive/robust-mrcpsp/code/plots/')
