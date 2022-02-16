import os
import re
import math
import argparse
import multiprocessing
from filelock import FileLock
from functools import partial
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
    instance_files = [instance for instance in os.listdir(instance_dir)][:10]
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


def benders_solve_and_write(instance, Gamma, time_limit):
    """
    Solves uncertain MRCPSP instance using Bender' decomposition and writes result to an (existing) output file.

    :param instance: Uncertain MRCPSP instance to solve.
    :type instance: Instance
    :param Gamma: Robustness parameter to use.
    :type Gamma: int
    :param time_limit: Time limit in seconds to give to Benders' algorithm for each instance.
    :type time_limit: int
    """
    current_dir = os.path.dirname(os.path.realpath(__file__))
    results_file = os.path.join(current_dir, 'benders_results_{}.txt'.format(Gamma))

    # solve instance using Benders' and write result
    benders_sol = Benders(instance, Gamma, time_limit).solve()
    gap = abs(benders_sol['objbound'] - benders_sol['objval']) / abs(benders_sol['objval'])  # optimality gap
    with FileLock(results_file + '.lock'):  # use FileLock to prevent two processes opening file at same time
        with open(results_file, 'a') as f:
            f.write("{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t{:.3f}\n".format(instance.name, Gamma,
                                                                                  benders_sol['objval'],
                                                                                  benders_sol['objbound'], gap,
                                                                                  benders_sol['n_iterations'],
                                                                                  benders_sol['avg_iteration'],
                                                                                  benders_sol['runtime']))


def parallel_benders_experiment(instances, Gamma, time_limit, num_threads):
    """
    Solves set of uncertain MRCPSP instances using Benders' decomposition. Available CPUs is split into parallel
    processes with each process being allocated 4 cores. Results are written to output file.

    :param instances: List of uncertain MRCPSP instances to solve.
    :type instances: list
    :param Gamma: Robustness parameter, i.e. max number of jobs that can achieve worst-case durations simultaneously.
        Takes an integer value from 0 to n (= number of non-dummy jobs in instance).
    :type Gamma: int
    :param time_limit: Time limit in seconds to allow for solving each instance.
    :type time_limit: int
    :param num_threads: Total number of available cores.
    :type num_threads: int
    """

    # create results file
    current_dir = os.path.dirname(os.path.realpath(__file__))
    results_file = os.path.join(current_dir, 'benders_results_{}.txt'.format(Gamma))
    with open(results_file, 'w+') as f:
        f.write("instance\tGamma\tobjval\tobjbound\tgap\tn_iterations\tavg_iteration\truntime\n")

    # set off multiple processes
    p = multiprocessing.Pool(math.floor(num_threads / 4))
    p.map(partial(benders_solve_and_write, Gamma=Gamma, time_limit=time_limit), instances)

    # delete FileLock .lock file
    os.remove(results_file + '.lock')


def compact_reformulation_solve_and_write(instance, Gamma, time_limit):
    """
    Solves uncertain MRCPSP instance using compact reformulation and writes result to an (existing) output file.

    :param instance: Uncertain MRCPSP instance to solve.
    :type instance: Instance
    :param Gamma: Robustness parameter to use.
    :type Gamma: int
    :param time_limit: Time limit in seconds to give to Benders' algorithm for each instance.
    :type time_limit: int
    """
    print('starting', instance.name)
    # create results file
    current_dir = os.path.dirname(os.path.realpath(__file__))
    results_file = os.path.join(current_dir, 'compact_reformulation_results_{}.txt'.format(Gamma))

    # solve instance using compact reformulation and write results
    compact_sol = CompactRefomulation(instance, Gamma, time_limit).solve()
    gap = abs(compact_sol['objbound'] - compact_sol['objval']) / abs(compact_sol['objval'])  # optimality gap
    with FileLock(results_file + '.lock'):  # use FileLock to prevent two processes opening file at same time
        with open(results_file, 'a') as f:
            f.write("{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n".format(instance.name, Gamma, compact_sol['objval'],
                                                                      compact_sol['objbound'], gap,
                                                                      compact_sol['runtime']))
    print('finishing', instance.name)


def parallel_compact_reformulation_experiment(instances, Gamma, time_limit, num_threads):
    """
    Solves set of uncertain MRCPSP instances using compact reformulation. Available CPUs is split into parallel
    processes with each process being allocated 4 cores. Results are written to output file.

    :param instances: List of uncertain MRCPSP instances to solve.
    :type instances: list
    :param Gamma: Robustness parameter, i.e. max number of jobs that can achieve worst-case durations simultaneously.
        Takes an integer value from 0 to n (= number of non-dummy jobs in instance).
    :type Gamma: int
    :param time_limit: Time limit in seconds to allow for solving each instance.
    :type time_limit: int
    :param num_threads: Total number of available cores.
    :type num_threads: int
    """
    # create results file
    current_dir = os.path.dirname(os.path.realpath(__file__))
    results_file = os.path.join(current_dir, 'compact_reformulation_results_{}.txt'.format(Gamma))
    with open(results_file, 'w+') as f:
        f.write("instance\tGamma\tobjval\tobjbound\tgap\truntime\n")

    # set off multiple processes
    p = multiprocessing.Pool(math.floor(num_threads / 4))
    p.map(partial(compact_reformulation_solve_and_write, Gamma=Gamma, time_limit=time_limit), instances)

    # delete FileLock .lock file
    os.remove(results_file + '.lock')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("solve_method",
                        choices={'benders', 'compact_formulation'},
                        help="Method to solve instances with. Either 'benders' or 'compact_reformulation'.")
    parser.add_argument("Gamma",
                        type=int,
                        help="Robustness parameter.")
    parser.add_argument("time_limit",
                        type=int,
                        help="Time limit in seconds to use for each instance.")
    parser.add_argument("num_threads",
                        type=int,
                        help='Total number of threads to use across all processes. Must be at least 4 since Gurobi uses'
                             '4 threads in solve methods.')
    args = parser.parse_args()

    # check num_threads is at least 4
    if args.num_threads < 4:
        raise argparse.ArgumentTypeError("num_threads must be at least 4.")

    instances = create_instances('/home/boldm1/OneDrive/robust-mrcpsp/instances/j10.mm', uncertainty_level=0.7)

    if args.solve_method == 'benders':
        parallel_benders_experiment(instances, args.Gamma, args.time_limit, args.num_threads)
    elif args.solve_method == 'compact_reformulation':
        parallel_compact_reformulation_experiment(instances, args.Gamma, args.time_limit, args.num_threads)
