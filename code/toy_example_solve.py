import mrcpsp
from compact_reformulation import CompactRefomulation
from benders import Benders


def solve_toy_example(nominal_data_file, Gamma, max_durational_deviations=None):
    """
    Solves toy example given in nominal_data_file for given Gamma. Max durational deviations can be explicitly provided
    for each activity in the input, but if not provided they are set to 0.5 * nominal values.
    """
    # load instance
    instance = mrcpsp.load_nominal_mrcpsp(nominal_data_file)
    if max_durational_deviations:
        instance.set_dbar_explicitly(max_durational_deviations)
    else:
        instance.set_dbar_uncertainty_level(0.5)

    print(f"\nSolving toy example given in {nominal_data_file}...")
    print("\nJob data [r, d, d_bar]:")
    print("-------------------------")
    for i in instance.V:
        print("job {}:".format(i),
              [[instance.jobs[i].r[m][0], instance.jobs[i].d[m], instance.jobs[i].d_bar[m]] for m in
               instance.jobs[i].M])

    print("\nCompact reformulation")
    print("---------------------")
    compact_sol = CompactRefomulation(instance, Gamma, time_limit=20).solve(print_log=False)
    print("objval:", compact_sol['objval'])
    print("runtime:", compact_sol['runtime'])
    print("modes:", compact_sol['modes'])
    print("network:", compact_sol['network'])
    print("resource flows:", compact_sol['flows'])

    print("\nBenders'")
    print("--------")
    benders_sol = Benders(instance, Gamma, time_limit=20).solve(print_log=True)
    print("objval:", benders_sol['objval'])
    print("runtime:", benders_sol['runtime'])
    print("n_iterations:", benders_sol['n_iterations'])
    print("modes:", benders_sol['modes'])
    print("network:", compact_sol['network'])
    print("resource flows:", benders_sol['flows'])


if __name__ == "__main__":

    # toy example from Balouka & Cohen (2021)
    balouka_toy_example = '/home/boldm1/OneDrive/robust-mrcpsp/code/instances/mrcpsp_toy_example_data.txt'
    balouka_max_durational_deviations = {0: [0], 1: [0, 3], 2: [8, 2], 3: [1], 4: [3], 5: [0], 6: [0]}

    # toy example from Bruni et al. (2017)
    bruni_toy_example = '/home/boldm1/OneDrive/robust-mrcpsp/code/instances/rcpsp_toy_example_data.txt'

    # toy example for my paper
    paper_toy_example = '/home/boldm1/OneDrive/robust-mrcpsp/code/instances/paper_toy_example_data.txt'
    zero_durational_deviations = {0: [0], 1: [0], 2: [0, 0], 3: [0], 4: [0], 5: [0, 0], 6: [0]}

    solve_toy_example(paper_toy_example, 2)
