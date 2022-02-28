

def read_benders_results(results_file):
    """
    Reads text file of Benders' results from results_file and returns dictionary of results data.

    :param results_file: Absolute filepath to results file produced by running a Benders' experiment.
    :type results_file: str
    :return: Dictionary of results data.
    :rtype: dict
    """
    print(results_file)
    results = {}
    with open(results_file, 'r') as f:
        raw_lines = f.read().splitlines()

    for line in raw_lines[1:]:
        splitline = line.split()
        print(raw_lines)

read_benders_results('/home/boldm1/OneDrive/robust-mrcpsp/code/results/storm_j20_7200/benders_results_0.txt')
