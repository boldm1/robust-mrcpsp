from benders import Benders
from compact_reformulation import CompactRefomulation


def solve(instance, solve_method, Gamma, time_limit):
    """
    Solves given uncertain MRCPSP instance using solve_method ('compact_reformulation' or 'benders') and prints solve
    log to terminal. Used to test code for specific instances.

    :param instance: Uncertain MRCPSP instance.
    :type instance: Instance
    :param solve_method: Solution method to use, i.e. 'compact_reformulation' or 'benders'.
    :type solve_method: str
    :param Gamma: Robustness parameter, i.e. max number of jobs that can achieve worst-case durations simultaneously.
        Takes integer value from 0 to n (= number of non-dummy jobs in instance).
    :type Gamma: int
    :param time_limit: Time limit (in seconds) to give to chosen solution method.
    :type time_limit: int
    """
    # check solve_method is recognised
    if solve_method not in ['compact_reformulation', 'benders']:
        raise Exception("{} is not a recognised solve method. Available methods are 'compact_reformulation' and "
                        "'benders'. Please check this input and try again.".format(solve_method))

    # solve instance using compact reformulation
    elif solve_method == 'compact_reformulation':
        print("Solving {} using compact reformulation:".format(instance.name))
        print('"""""""""""""""""""""""""""""""""""""""""""""\n')
        sol = CompactRefomulation(instance, Gamma, time_limit).solve(print_log=True)

    # solve instance using Benders' decomposition
    elif solve_method == 'benders':
        print("Solving {} using Benders' decomposition:".format(instance.name))
        print('"""""""""""""""""""""""""""""""""""""""""""""\n')
        sol = Benders(instance, Gamma, time_limit).solve(print_log=True)

    print(sol)
    print("objval:", sol['objval'])
    print("objbound:", sol['objbound'])
    print("runtime:", sol['runtime'])
