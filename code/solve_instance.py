import mrcpsp
from benders import Benders
from compact_reformulation import CompactRefomulation


def solve(nominal_data_file, solve_method, Gamma, time_limit, max_durational_deviations=None, print_log=False):
    """
    Solve instance given in nominal_data_file with Gamma using specified solve method ('compact_reformulation' or
    'benders'). Max durational deviations can be explicitly provided for each activity in the input, but if not provided
    they are set to floor(0.7 * nominal values). This function is used to test code for specific instances.

    :param nominal_data_file: Absolute file path to deterministic MRCPSP instance from which to create robust instance.
    :type nominal_data_file: str
    :param solve_method: Solution method to use, i.e. 'compact_reformulation' or 'benders'.
    :type solve_method: str
    :param Gamma: Robustness parameter, i.e. max number of jobs that can achieve worst-case durations simultaneously.
        Takes integer value from 0 to n (= number of non-dummy jobs in instance).
    :type Gamma: int
    :param time_limit: Time limit (in seconds) to give to chosen solution method.
    :type time_limit: int
    :param max_durational_deviations: Optional parameter. Dictionary with max durational deviation for each mode of each
        activity.
    :type max_durational_deviations: dict
    :param print_log: Indicates whether or not to print solve log to terminal. Defaults to False.
    :type print_log: bool
    """
    # load instance
    instance = mrcpsp.load_nominal_mrcpsp(nominal_data_file)
    if max_durational_deviations:
        instance.set_dbar_explicitly(max_durational_deviations)
    else:
        instance.set_dbar_uncertainty_level(0.7)

    # check solve_method is recognised
    if solve_method not in ['compact_reformulation', 'benders']:
        raise Exception("{} is not a recognised solve method. Available methods are 'compact_reformulation' and "
                        "'benders'. Please check this input and try again.".format(solve_method))
    # solve instance using compact reformulation
    elif solve_method == 'compact_reformulation':
        print("Solving {} using compact reformulation:".format(instance.name))
        print('"""""""""""""""""""""""""""""""""""""""""""""\n')
        compact_sol = CompactRefomulation(instance, Gamma, time_limit).solve(print_log=print_log)
        print("objval:", compact_sol['objval'])
        print("runtime:", compact_sol['runtime'])
        print("modes:", compact_sol['modes'])
        print("network:", compact_sol['network'])
        print("resource flows:", compact_sol['flows'])

    # solve instance using Benders' decomposition
    elif solve_method == 'benders':
        print("Solving {} using Benders' decomposition:".format(instance.name))
        print('"""""""""""""""""""""""""""""""""""""""""""""\n')
        benders_sol = Benders(instance, Gamma, time_limit).solve(print_log=print_log)
        print("objval:", benders_sol['objval'])
        print("runtime:", benders_sol['runtime'])
        print("n_iterations:", benders_sol['n_iterations'])
        print("modes:", benders_sol['modes'])
        print("network:", benders_sol['network'])
        print("resource flows:", benders_sol['flows'])


if __name__ == "__main__":

    # toy example from Balouka & Cohen (2021)
    balouka_toy_example = '/home/boldm1/OneDrive/robust-mrcpsp/code/instances/mrcpsp_toy_example_data.txt'
    balouka_max_durational_deviations = {0: [0], 1: [0, 3], 2: [8, 2], 3: [1], 4: [3], 5: [0], 6: [0]}

    # toy example from Bruni et al. (2017)
    bruni_toy_example = '/home/boldm1/OneDrive/robust-mrcpsp/code/instances/rcpsp_toy_example_data.txt'
    bruni_max_durational_deviations = {0: [0], 1: [2], 2: [1], 3: [3], 4: [2], 5: [1], 6: [1], 7: [1], 8: [0]}

    # toy example for my paper
    paper_toy_example = '/home/boldm1/OneDrive/robust-mrcpsp/code/instances/paper_toy_example_data.txt'
    zero_durational_deviations = {0: [0], 1: [0], 2: [0, 0], 3: [0], 4: [0], 5: [0, 0], 6: [0]}

    test_instance = '/home/boldm1/OneDrive/robust-mrcpsp/code/instances/j10.mm/j1014_2.mm'

    solve(paper_toy_example, 'benders', 3, 60, max_durational_deviations=zero_durational_deviations, print_log=True)
