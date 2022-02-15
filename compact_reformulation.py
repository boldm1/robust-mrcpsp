import gurobipy as gp
from gurobipy import GRB


class CompactRefomulation():
    """
    Class to solve uncertain MRCPSP instance using compact reformulation.
    """

    def __init__(self, instance, Gamma, time_limit):
        """
        Initialises CompactReformulation object with instance to solve and solve parameters.

        :param instance: Robust MRCPSP instance to solve
        :type instance: Instance
        :param Gamma: Robustness parameter, i.e. max number of jobs that can achieve worst-case durations
                simultaneously. Takes an integer value from 0 to n.
        :type Gamma: int
        :param time_limit: Time limit for Gurobi, in seconds
        :type time_limit: int
        :return: Dictionary containing solution information.
        :rtype: dict
        """
        self.instance = instance
        self.Gamma = Gamma
        self.time_limit = time_limit

    def solve(self, print_log=False):
        """
        Builds and solves compact reformulation for uncertain MRCPSP. Terminates upon finding an optimal solution or
        reaching the specified time limit.

        :param print_log: Indicates whether or not to print Gurobi solve log to terminal. Defaults to False.
        :type print_log: bool
        :return: Dictionary containing solution information.
        :rtype: dict
        """
        # get instance data
        V = self.instance.V
        n = self.instance.n
        E = self.instance.E
        d = [self.instance.jobs[j].d for j in V]  # d[j][m] = duration of j in mode m
        d_bar = [self.instance.jobs[j].d_bar for j in V]  # d_bar[j][m] = max durational deviation of j in mode m
        M = [self.instance.jobs[j].M for j in V]  # M[j] = list of modes of j
        # r[j][m][k] = amount of resource k required by j in mode m. Dummy source and sink all available resource.
        r = [[self.instance.R]] + [self.instance.jobs[j].r for j in self.instance.N] + [[self.instance.R]]

        # Create model
        model_name = '{}_compact_reformulation'.format(self.instance.name)
        model = gp.Model(model_name)

        # Set Gurobi parameters
        model.setParam('OutputFlag', print_log)
        model.setParam('TimeLimit', self.time_limit)
        model.setParam('Threads', 4)  # always limit number of threads to 4

        G = range(self.Gamma + 1)  # levels in auxiliary graph

        # Create variables
        S = model.addVars([(i, g) for i in V for g in G], name='S', vtype=GRB.CONTINUOUS, lb=0)
        y = model.addVars([(i, j) for i in V for j in V], name="y",
                          vtype=GRB.BINARY)  # y_ij = 1 if i -> j in transitive project precedence network
        x = model.addVars([(i, m) for i in V for m in M[i]], name="x",
                          vtype=GRB.BINARY)  # x_im = 1 if i is executed in mode m
        f = model.addVars([(i, j, k) for i in V for j in V for k in self.instance.K_renew], name="f",
                          vtype=GRB.CONTINUOUS, lb=0)  # f_ijk = flow from i to j of renewable resource k

        # set model objective
        model.setObjective(S[n + 1, self.Gamma], GRB.MINIMIZE)

        # Precedence constraints
        model.addConstr(S[0, 0] == 0)
        BigM = sum(max(d[i][m] + d_bar[i][m] for m in M[i]) for i in V)  # largest possible makespan
        model.addConstrs(
            S[j, g] - S[i, g] >= d[i][m] * x[i, m] - BigM * (1 - y[i, j]) for i in V for j in V for m in M[i] for g in
            G)
        model.addConstrs(
            S[j, g + 1] - S[i, g] >= (d[i][m] + d_bar[i][m]) * x[i, m] - BigM * (1 - y[i, j]) for i in V for j in V for
            m in M[i] for g in G[:-1])
        model.addConstrs(y[e[0], e[1]] == 1 for e in E)
        model.addConstr(y[n + 1, n + 1] == 1)
        # Renewable resource constraints
        BigR = [
            [[min(max(r[i][m][k] for m in M[i]), max(r[j][m][k] for m in M[j])) for k in self.instance.K_renew] for j in
             V] for i in V]
        model.addConstrs(
            f[i, j, k] <= BigR[i][j][k] * y[i, j] for i in V for j in V if i != n + 1 if j != 0 for k in
            self.instance.K_renew)
        model.addConstrs(
            gp.quicksum(f[i, j, k] for j in V if j != 0) == gp.quicksum(r[i][m][k] * x[i, m] for m in M[i]) for i in V
            if i != n + 1 for k in self.instance.K_renew)
        model.addConstrs(
            gp.quicksum(f[i, j, k] for i in V if i != n + 1) == gp.quicksum(
                r[j][m][k] * x[j, m] for m in M[j]) for j in V if j != 0 for k in self.instance.K_renew)
        model.addConstrs(gp.quicksum(x[i, m] for m in M[i]) == 1 for i in V)
        # Non-renewable resource constraints
        model.addConstrs(gp.quicksum(r[i][m][k] * x[i, m] for m in M[i]) <= self.instance.R[k] for i in V for k in
                         self.instance.K_nonrenew)

        # solve model
        model.optimize()

        # dict to store solution info
        modes = {i: int(sum(m * x[i, m].X for m in M[i])) for i in V}
        T_E = [(i, j) for i in V for j in V if i != j if y[i, j].X > 0.99]
        flows = {(i, j): [f[i, j, k].X for k in self.instance.K_renew] for i in V for j in V if i != n + 1 if j != 0 if
                 i != j}
        return {'objval': model.ObjVal, 'objbound': model.ObjBound, 'runtime': model.Runtime, 'modes': modes,
                'network': T_E, 'flows': flows}
