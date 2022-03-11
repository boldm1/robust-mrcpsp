import os
import re
import argparse
from mrcpsp import load_nominal_mrcpsp
from benders import Benders
from compact_reformulation import CompactRefomulation


def create_instances(instance_dir, uncertainty_level):
    """
    Returns list of uncertain MRCPSP instances generated from nominal MRCPSP instances contained in instance_dir. Max
    durational deviations are set according to uncertainty_level.

    :param instance_dir: Absolute filepath of directory containing nominal MRCPSP instances downloaded from PSPLIB.
    :type instance_dir: str
    :param uncertainty_level: Non-negative value representing desired ratio between d_bar[m] / d[m].
    :type uncertainty_level: float
    :return: List of uncertain MRCPSP instances.
    :rtype: list
    """
    instance_files = [instance for instance in os.listdir(instance_dir)]
    # sort instance files into their natural order
    instance_files = sorted(instance_files, key=lambda f: 10 * int(''.join(re.findall(r'\d+', f[3:5]))) + int(
        ''.join(re.findall(r'\d+', f[5:8]))))

    instances = []
    for file in instance_files:
        # load nominal instance from data file
        instance_filepath = os.path.join(instance_dir, file)
        instance = load_nominal_mrcpsp(instance_filepath)
        # set durational uncertainty
        instance.set_dbar_uncertainty_level(uncertainty_level)
        # append to list of uncertain MRCPSP instances
        instances.append(instance)
    return instances


def run_experiment(instances, Gamma, solve_method, time_limit):
    """
    Solves a set of uncertain MRCPSP instances using specified solve method and writes results to an output file.

    :param instances: List of uncertain MRCPSP instances to solve.
    :type instances: list
    :param Gamma: Robustness parameter to use.
    :type Gamma: int
    :param solve_method: Method to run experiment for. Either 'benders' or 'compact_reformulation'.
    :type solve_method: str
    :param time_limit: Time limit in seconds to give to Benders' algorithm for each instance.
    :type time_limit: int
    """
    current_dir = os.path.dirname(os.path.realpath(__file__))

    # Check if results already exist. Get solved instances if it does, or create new results file if not
    results_file = os.path.join(current_dir, f'{solve_method}_results_{Gamma}.txt')
    if os.path.exists(results_file):
        with open(results_file, 'r') as f:
            raw_lines = f.read().splitlines(keepends=True)
        solved_instances = []
        lines_to_keep = []
        first_line = raw_lines[0]
        for line in raw_lines[1:]:
            splitline = line.split()
            instance = splitline[0]
            gap = float(splitline[4])
            if solve_method == 'benders':
                runtime = float(splitline[7])
            elif solve_method == 'compact_reformulation':
                runtime = float(splitline[5])
            if (gap == 0) or (runtime > time_limit - 0.001):
                solved_instances.append(instance)
                lines_to_keep.append(line)
        instances = [i for i in instances if i.name not in solved_instances]  # instances still to be solved
        # re-write results file, removing partial solved instances
        os.remove(results_file)
        with open(results_file, 'w+') as f:
            f.write(first_line)
            for line in lines_to_keep:
                f.write(line)
    else:
        # create results file
        with open(results_file, 'w+') as f:
            if solve_method == 'benders':
                f.write("instance\tGamma\tobjval\tobjbound\tgap\tn_iterations\tavg_iteration\truntime\n")
            elif solve_method == 'compact_reformulation':
                f.write("instance\tGamma\tobjval\tobjbound\tgap\truntime\n")

    # solve each instance using chosen solve method and write results
    for instance in instances:
        if solve_method == 'benders':
            benders_sol = Benders(instance, Gamma, time_limit).solve()
            gap = abs(benders_sol['objbound'] - benders_sol['objval']) / abs(benders_sol['objval'])  # optimality gap
            with open(results_file, 'a') as f:
                f.write(
                    "{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t{:.3f}\n".format(instance.name, Gamma,
                                                                                  benders_sol['objval'],
                                                                                  benders_sol['objbound'], gap,
                                                                                  benders_sol['n_iterations'],
                                                                                  benders_sol['avg_iteration'],
                                                                                  benders_sol['runtime']))
        if solve_method == 'compact_reformulation':
            compact_sol = CompactRefomulation(instance, Gamma, time_limit).solve()
            gap = abs(compact_sol['objbound'] - compact_sol['objval']) / abs(compact_sol['objval'])  # optimality gap
            with open(results_file, 'a') as f:
                f.write("{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n".format(instance.name, Gamma, compact_sol['objval'],
                                                                          compact_sol['objbound'], gap,
                                                                          compact_sol['runtime']))


def is_valid_directory(arg):
    """
    Checks directory path given as input argument exists. Throws an error if it does not.
    """
    if os.path.isdir(arg):
        return arg
    else:
        raise argparse.ArgumentTypeError("{} is not a valid directory.".format(arg))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-instance_dir",
                        type=is_valid_directory,
                        required=True,
                        help="Path to directory containing nominal MRCPSP instances.")
    parser.add_argument("-solve_method",
                        choices={'benders', 'compact_reformulation'},
                        required=True,
                        help="Method to solve instances with. Either 'benders' or 'compact_reformulation'.")
    parser.add_argument("-Gamma",
                        type=int,
                        required=True,
                        help="Robustness parameter.")
    parser.add_argument("-time_limit",
                        type=int,
                        required=True,
                        help="Time limit in seconds to use for each instance.")
    args = parser.parse_args()

    instances = create_instances(args.instance_dir, uncertainty_level=0.7)
    run_experiment(instances, args.Gamma, args.solve_method, args.time_limit)
