import os
import re
from mrcpsp import load_nominal_mrcpsp
from benders import Benders
from compact_reformulation import compact_reformulation


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


def benders_experiment(instances, Gamma, time_limit):
    """
    Solves set uncertain MRCPSP instances using Bender' decomposition and writes results to an output file.

    :param instances: List of uncertain MRCPSP instances to solve.
    :type instances: list
    :param Gamma: Robustness parameter to use.
    :type Gamma: int
    :param time_limit: Time limit in seconds to give to Benders' algorithm for each instance.
    :type time_limit: int
    """
    # create results file
    current_dir = os.path.dirname(os.path.realpath(__file__))
    results_file = os.path.join(current_dir, 'benders_results.txt')
    f1 = open(results_file, 'w+')
    f1.write("instance\tGamma\tobjval\tobjbound\tgap\tn_iterations\tavg_iteration\truntime\n")
    f1.close()

    # solve each instance using Benders' and write results
    for instance in instances:
        benders_sol = Benders(instance, Gamma).solve(time_limit)
        gap = abs(benders_sol['objbound'] - benders_sol['objval']) / abs(benders_sol['objval'])  # optimality gap
        f = open(results_file, 'a')
        f.write(
            "{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t{:.3f}\n".format(instance.name, Gamma, benders_sol['objval'],
                                                                          benders_sol['objbound'], gap,
                                                                          benders_sol['n_iterations'],
                                                                          benders_sol['avg_iteration'],
                                                                          benders_sol['runtime']))
        f.close()


def compact_reformulation_experiment(instances, Gamma, time_limit):
    """
    Solves set uncertain MRCPSP instances using compact reformulation and writes results to an output file.

    :param instances: List of uncertain MRCPSP instances to solve.
    :type instances: list
    :param Gamma: Robustness parameter to use.
    :type Gamma: int
    :param time_limit: Time limit in seconds to give to Benders' algorithm for each instance.
    :type time_limit: int
    """
    # create results file
    current_dir = os.path.dirname(os.path.realpath(__file__))
    results_file = os.path.join(current_dir, 'compact_reformulation_results.txt')
    f = open(results_file, 'w+')
    f.write("instance\tGamma\tobjval\tobjbound\tgap\truntime\n")
    f.close()

    # solve each instance using compact reformulation and write results
    for instance in instances:
        compact_sol = compact_reformulation(instance, Gamma, time_limit)
        gap = abs(compact_sol['objbound'] - compact_sol['objval']) / abs(compact_sol['objval'])  # optimality gap
        f = open(results_file, 'a')
        f.write("{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n".format(instance.name, Gamma, compact_sol['objval'],
                                                                  compact_sol['objbound'], gap, compact_sol['runtime']))
        f.close()


instances = create_instances('/home/boldm1/OneDrive/robust-mrcpsp/instances/j10.mm', uncertainty_level=0.7)
benders_experiment(instances, 1, 5*60)
#compact_reformulation_experiment(instances, 1, 5 * 60)
