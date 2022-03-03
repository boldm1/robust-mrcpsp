import re
import os
import math


def load_nominal_mrcpsp(filepath):
    """
    Loads MRCPSP instance from data file and returns Instance object.

    :param filepath: Filepath of MRCPSP instance to load
    :type filepath: str
    :return: MRCPSP instance
    :rtype: Instance
    """
    with open(filepath, 'r') as f:
        raw_lines = f.read().splitlines()

    block = 0
    succ = {}  # dict of successors of each job
    M = {}  # dict of modes for each job
    d = {}  # dict of duration for each job
    r = {}  # dict of resource requirements for each job
    for line in raw_lines:
        ints = list(map(int, re.findall(r'\d+', line)))  # list of integers in line
        if len(ints) == 0:
            if '*****' in line:
                block += 1
            continue
        elif block == 2:
            if 'jobs' in line:
                n = int(ints[0]) - 2  # number of jobs excl. dummy source and sink
            elif 'horizon' in line:
                T = int(ints[0])  # project horizon
            elif '- renewable' in line:
                n_renew = int(ints[0])  # number of renewable resources
            elif '- nonrenewable' in line:
                n_nonrenew = int(ints[0])  # number of non-renewable resources
        elif block == 4:
            jobnr = int(ints[0]) - 1  # jobs are numbered from 1 to n+2 in data files
            M[jobnr] = [m for m in range(int(ints[1]))]
            succ[jobnr] = [int(succ) - 1 for succ in ints[3:]]  # ints[3:] = successors
        elif block == 5:
            n_res = n_renew + n_nonrenew
            if len(ints) == 3 + n_res:  # mode 1 data
                jobnr = int(ints[0]) - 1
                d[jobnr] = [int(ints[2])]
                r[jobnr] = [[int(x) for x in ints[3:]]]
            elif len(ints) == 2 + n_res:  # modes 2,3,4... data
                d[jobnr].append(int(ints[1]))
                r[jobnr].append([int(x) for x in ints[2:]])
        if block == 6:
            R = [int(x) for x in ints]
    jobs = {}
    for j in range(n + 2):
        pred = [i for i in range(j) if j in succ[i]]
        jobs[j] = NominalJob(j, pred, succ[j], M[j], d[j], r[j])
    instance_name = os.path.basename(os.path.normpath(os.path.splitext(filepath)[0]))
    K_renew = [k for k in range(n_renew)]  # indices of renewables resources in R
    K_nonrenew = [k for k in range(n_renew, n_res)]  # indices of non-renewable resources in R
    instance = Instance(instance_name, jobs, R, K_renew, K_nonrenew, T)
    return instance


class Job:
    """
    Class representing multi-mode job from Robust Multi-mode Resource Constrained Project Scheduling Problem instance.
    """

    def __init__(self, idx, pred, succ, M, d, d_bar, r):
        """
        Initialises Job with data.

        :param idx: Job index
        :type idx: int
        :param pred: List of predecessors
        :type pred: list
        :param succ: List of successors
        :type succ: list
        :param M: Activity processing modes
        :type M: list
        :param d: Nominal duration of each mode in M
        :type d: list[int]
        :param d_bar: Maximal deviation from nominal duration, for each mode in M
        :type d_bar: list[int]
        :param r: Requirement for each resource of each mode in M
        :type r: list[list]
        """
        self.idx = idx
        self.pred = pred
        self.succ = succ
        self.M = M
        self.d = d
        self.d_bar = d_bar
        self.r = r

        # Earliest and latest start and finish times. Computed upon project creation.
        self.ES = None
        self.LS = None
        self.EF = None
        self.LF = None


class NominalJob(Job):
    """
    Class representing job with certain duration from MRCPSP. Inherits parameters from Job class.
    """

    def __init__(self, idx, pred, succ, M, d, r):
        """
        Initialises NominalJob with data.

        :param idx: Job index
        :type idx: int
        :param pred: List of predecessors
        :type pred: list
        :param succ: List of successors
        :type succ: list
        :param M: Activity processing modes
        :type M: list
        :param d: Nominal duration of each mode in M
        :type d: list[int]
        :param r: Requirement for each resource of each mode in M
        :type r: list[list]
        """
        super().__init__(idx, pred, succ, M, d, [0 for _ in M], r)


class Instance:
    """
    Class representing Robust Multi-mode Resource Constrained Project Scheduling Problem instance.
    """

    def __init__(self, name, jobs, R, K_renew, K_nonrenew, T):
        """
        Initialises MRCPSP instance with jobs and resource data.

        :param name: Instance name
        :type name: str
        :param jobs: Dictionary of jobs in instance. Keys are job indices, values are Job objects.
        :type jobs: dict
        :param R: Availability of renewable and non-renewable resources
        :type R: list
        :param K_renew: List of indices of renewable resource in R
        :type K_renew: list
        :param K_nonrenew: List of indices non-renewable resources in R
        :type K_nonrenew: list
        :param T: Project horizon
        :type T: int
        """
        self.name = name
        self.jobs = jobs
        self.R = R
        self.K_renew = K_renew
        self.K_nonrenew = K_nonrenew
        self.T = T

        self.V = [i for i in jobs]  # list of jobs, incl. dummy source and sink
        self.n = len(self.V) - 2  # number of jobs, excl. dummy source and sink
        self.N = [i for i in jobs if i != 0 if i != self.n + 1]  # list of jobs, excl. dummy source and sink
        self.E = [(i, j) for i in self.V for j in self.jobs[i].succ]  # list of finish-to-start precedences

        # Compute earliest and latest start and finish for each job
        self.forward_pass()
        self.backward_pass()

    def forward_pass(self):
        """
        Runs through jobs from 0 to n+1 computing ES and EF for each job. Assumes that jobs are numbered topologically
        (i.e. if i -> j, then i < j). This function is called upon instance creation.
        """
        self.jobs[0].ES = 0
        self.jobs[0].EF = 0
        for j in range(1, self.n + 2):
            self.jobs[j].ES = max(self.jobs[i].EF for i in self.jobs[j].pred)  # ES = max EF of predecessors
            self.jobs[j].EF = self.jobs[j].ES + self.jobs[j].d[0]  # shortest duration

    def backward_pass(self):
        """
        Runs through jobs from n+1 to 0 computing LF and LS for each job. Assumes that jobs are numbered topologically
        (i.e. if i -> j, then i < j). This function is called upon instance creation.
        """
        n = self.n
        self.jobs[n + 1].LF = self.T
        self.jobs[n + 1].LS = self.T
        for j in range(n, -1, -1):
            self.jobs[j].LF = min(self.jobs[i].LS for i in self.jobs[j].succ)  # LF = min LS of all successors
            self.jobs[j].LS = self.jobs[j].LF - self.jobs[j].d[0]  # shortest duration

    def set_dbar_explicitly(self, d_bar):
        """
        Function to explicitly set d_bar for each mode of jobs in Instance.

        :param d_bar: Dict with d_bar values for jobs in Instance. Values must be provided for each mode of each job in
            instance.
        :type d_bar: dict
        """
        # check each job has an entry in d_bar
        no_value = []
        for j in self.V:
            try:
                v = d_bar[j]
            except KeyError:
                no_value.append(j)
        if len(no_value) > 0:
            raise Exception(f"d_bar does not have an entry for job {str(no_value)[1:-1]}. Please correct this and try "
                            "again.")

        # check each job in d_bar has values for each of its possible modes
        to_fix = []
        for j in self.V:
            if len(d_bar[j]) != len(self.jobs[j].M):
                to_fix.append(j)
        # throw error if to_fix is non-empty
        if len(to_fix) > 0:
            raise Exception(f"The number of values specified in d_bar for jobs {str(to_fix)[1:-1]} does not match the "
                            "number of modes for those jobs! Please correct this and try again.")

        # set d_bar
        for j in d_bar:
            self.jobs[j].d_bar = d_bar[j]

    def set_dbar_uncertainty_level(self, uncertainty_level):
        """
        Sets d_bar according to target uncertainty_level for each mode of each job.

        :param uncertainty_level: Non-negative value representing desired ratio between d_bar[m] / d[m].
        :type uncertainty_level: float
        """
        for j in self.V:
            self.jobs[j].d_bar = [math.floor(uncertainty_level * self.jobs[j].d[m]) for m in self.jobs[j].M]
