import gurobipy as gp
from gurobipy import GRB


def compact_reformulation(instance, Gamma, time_limit):
    """
    Solves given robust MRCPSP instance using compact mixed-integer programming reformulation.

    :param instance: Robust MRCPSP instance to solve
    :type instance: Instance
    :param Gamma: Robustness parameter, i.e. max number of jobs that can achieve worst-case durations
            simultaneously. Takes an integer value from 0 to n.
    :type Gamma: int
    :param time_limit: Time limit for Gurobi, in seconds
    :type time_limit: int
    :return: Solution object
    :rtype: Solution
    """
    # get instance data
    d = [instance.jobs[j].d for j in instance.V]  # d[j][m] = duration of j in mode m
    d_bar = [instance.jobs[j].d_bar for j in instance.V]  # d_bar[j][m] = max durational deviation of j in mode m
    M = [instance.jobs[j].M for j in instance.V]  # M[j] = list of modes of j
    # r[j][m][k] = amount of resource k required by j in mode m. Dummy source and sink all available resource.
    r = [[instance.R]] + [instance.jobs[j].r for j in instance.N] + [[instance.R]]

    # Create model
    model_name = '{}_compact_reformulation'.format(instance.name)
    model = gp.Model(model_name)

    # Set Gurobi parameters
    model.setParam('TimeLimit', time_limit)

    G = range(Gamma + 1)  # levels in auxiliary graph

    # Create variables
    S = model.addVars([(i, g) for i in instance.V for g in G], name='S', vtype=GRB.INTEGER, lb=0)
    y = model.addVars([(i, j) for i in instance.V for j in instance.V], name="y",
                      vtype=GRB.BINARY)  # y_ij = 1 if i -> j in transitive project precedence network
    x = model.addVars([(i, m) for i in instance.V for m in M[i]], name="x",
                      vtype=GRB.BINARY)  # x_im = 1 if i is executed in mode m
    f = model.addVars([(i, j, k) for i in instance.V for j in instance.V for k in instance.K_renew], name="f",
                      vtype=GRB.CONTINUOUS, lb=0)  # f_ijk = flow from i to j of renewable resource k

    # set model objective
    model.setObjective(S[instance.n + 1, Gamma], GRB.MINIMIZE)

    # Precedence constraints
    model.addConstr(S[0, 0] == 0)
    BigM = sum(max(d[i][m] + d_bar[i][m] for m in M[i]) for i in instance.V)  # largest possible makespan
    model.addConstrs(
        S[j, g] - S[i, g] >= d[i][m] * x[i, m] - BigM * (1 - y[i, j]) for i in instance.V for j in instance.V for m in
        M[i] for g in G)
    model.addConstrs(
        S[j, g + 1] - S[i, g] >= (d[i][m] + d_bar[i][m]) * x[i, m] - BigM * (1 - y[i, j]) for i in instance.V for j in
        instance.V for m in
        M[i] for g in G[:-1])
    model.addConstrs(y[e[0], e[1]] == 1 for e in instance.E)
    model.addConstr(y[instance.n + 1, instance.n + 1] == 1)
    # Renewable resource constraints
    BigR = [[[min(max(r[i][m][k] for m in M[i]), max(r[j][m][k] for m in M[j])) for k in instance.K_renew] for j in
             instance.V] for i in instance.V]
    model.addConstrs(
        f[i, j, k] <= BigR[i][j][k] * y[i, j] for i in instance.V for j in instance.V if i != instance.n + 1 if j != 0
        for k in instance.K_renew)
    model.addConstrs(
        gp.quicksum(f[i, j, k] for j in instance.V) == gp.quicksum(r[i][m][k] * x[i, m] for m in M[i]) for i in
        instance.V if i != instance.n + 1 for k in instance.K_renew)
    model.addConstrs(
        gp.quicksum(f[i, j, k] for i in instance.V) == gp.quicksum(r[j][m][k] * x[j, m] for m in M[j]) for j in
        instance.V if j != 0 for k in instance.K_renew)
    model.addConstrs(gp.quicksum(x[i, m] for m in M[i]) == 1 for i in instance.V)

    # Non-renewable resource constraints
    model.addConstrs(gp.quicksum(r[i][m][k] * x[i, m] for m in M[i]) for i in instance.V for k in instance.K_nonrenew)

    # solve model
    model.optimize()
    model.write('{}.sol'.format(model_name))


class Solution:
    """
    Class representing a solution to the robust MRCPSP.
    """

    def __init__(self, instance, S, M):
        """
        Initialises Solution object with job start times and processing modes.

        :param instance: MRCPSP instance that has been solved
        :type instance: Instance
        :param S: Dictionary of start times for each job
        :type S: dict
        :param M: Dictionary of chosen processing mode for each job
        :type M: dict
        """
        self.instance = instance
        self.S = S
        self.M = M


class Benders:
    """
    Class representing a Bender' style decomposition approach for solving robust MRCPSP.
    """

    def __init__(self, instance, Gamma):
        """

        :param instance: Robust MRCPSP instance to solve
        :type instance: Instance
        :param Gamma: Robustness parameter, i.e. max number of jobs that can achieve worst-case durations
            simultaneously. Takes an integer value from 0 to n.
        :type Gamma: int
        """
        self.instance = instance
        self.Gamma = Gamma

        self.LB = float('-inf')  # LB on objective function of robust MRCPSP instance
        self.UB = float('inf')  # UB on objective function of robust MRCPSP instance
        t = 0  # iteration number
        cuts = []  # cuts generated by subproblem to be applied to master problem

        # initialised in create_master_problem
        self.master_model = None
        self.master_eta = None
        self.master_y = None
        self.master_x = None
        self.master_f = None
        self.master_s = None
        self.master_phi = None

    def benders_solve(self, instance, Gamma, time_limit):
        """
        Iterates between master problem and subproblem, adding cuts after each iteration, until LB >= UB.

        :return:
        :rtype:
        """
        LB = float('-inf')  # LB on objective function of robust MRCPSP instance
        UB = float('inf')  # UB on objective function of robust MRCPSP instance

        t = 0  # iteration number
        cuts = []  # cuts generated by subproblem to be applied to master problem

        while LB < UB:
            master_solution = self.master_problem(t, cuts)
            subproblem_objval = self.subproblem(t, master_solution)
            if subproblem_objval < UB:
                UB = subproblem_objval

    def create_master_problem(self):
        """
        Selects job processing modes and resource allocations to minimise project duration, under nominal job processing
        times. Provides a LB to optimal robust solution.

        :return: Gurobi model populated with variables, constraints and objective function, representing basic master
            problem without optimality cuts
        :rtype: Model
        """
        # get instance data
        V = self.instance.V
        n = self.instance.n
        d = [self.instance.jobs[j].d for j in V]
        d_bar = [self.instance.jobs[j].d_bar for j in V]
        M = [self.instance.jobs[j].M for j in V]
        # Dummy source and sink all available resource.
        r = [[self.instance.R]] + [self.instance.jobs[j].r for j in self.instance.N] + [[self.instance.R]]

        # Create model
        model_name = '{}_benders_master'.format(self.instance.name)
        model = gp.Model(model_name)

        # Create variables
        eta = model.addVar(name="eta")  # variable to represent objective value
        y = model.addVars([(i, j) for i in V for j in V], name="y",
                          vtype=GRB.BINARY)  # y_ij = 1 if i -> j in transitive project precedence network
        x = model.addVars([(i, m) for i in V for m in M[i]], name="x",
                          vtype=GRB.BINARY)  # x_im = 1 if i is executed in mode m
        f = model.addVars([(i, j, k) for i in V for j in V for k in self.instance.K_renew], name="f",
                          vtype=GRB.CONTINUOUS, lb=0)  # f_ijk = flow from i to j of renewable resource k
        s = model.addVars([i for i in V], name="s", vtype=GRB.INTEGER)  # start time of each job
        phi = model.addVars([(i, j) for i in V for j in V], name="phi", vtype=GRB.CONTINUOUS, lb=0)  # longest path vars

        # set model objective
        model.setObjective(eta, GRB.MINIMIZE)

        # Precedence constraints
        model.addConstrs(y[e[0], e[1]] == 1 for e in self.instance.E)
        model.addConstrs(y[i, j] + y[j, i] <= 1 for i in V for j in V if i < j)
        model.addConstrs(
            y[i, p] >= y[i, j] + y[j, p] - 1 for i in V for j in V for p in V if i != j if i != p if j != p)
        # Renewable resource constraints
        model.addConstrs(
            f[i, j, k] <= min(r[i][m_i][k], r[j][m_j][k]) * y[i, j] +
            (1 - x[i, m_i]) * (
                    max(max(r[i][m][k] for m in M[i]), max(r[j][m][k] for m in M[j])) -
                    min(r[i][m_i][k], r[j][m_j][k])) +
            (1 - x[j, m_j]) * (
                    max(max(r[i][m][k] for m in M[i]), max(r[j][m][k] for m in M[j])) -
                    min(r[i][m_i][k], r[j][m_j][k])) for i in
            V if i != n + 1
            for j in V if j != 0 for m_i in M[i] for m_j in M[j] for k in self.instance.K_renew)
        model.addConstrs(gp.quicksum(x[i, m] for m in M[i]) == 1 for i in V)
        model.addConstrs(
            gp.quicksum(f[i, j, k] for j in V if j != i if j != 0) == gp.quicksum(r[i][m][k] * x[i, m] for m in M[i])
            for i in V if i != n + 1 for k in self.instance.K_renew)
        model.addConstrs(
            gp.quicksum(f[i, j, k] for i in V if i != j if i != n + 1) == gp.quicksum(
                r[j][m][k] * x[j, m] for m in M[j]) for j in V if j != 0 for k in self.instance.K_renew)

        # UB on min makespan, i.e. largest possible gap between s[j] and s[i].
        N = sum(max(d[i][m] + d_bar[i][m] for m in M[i]) for i in V)
        model.addConstrs(
            s[j] - s[i] >= d[i][m] - N * (1 - y[i, j]) - N * (1 - x[i, m]) for i in V for j in V if i != j for m in
            M[i])
        model.addConstr(eta >= s[n])

        # longest path with minimal job durations
        model.addConstr(eta >= gp.quicksum(min(d[i][m] for m in M[i]) * phi[i][j] for i in V for j in V))
        model.addConstrs(
            gp.quicksum(phi[j][i] for j in V) == gp.quicksum(phi[i][j] for j in V) for i in V if i != 0 if i != n + 1)
        model.addConstr(gp.quicksum(phi[0][j] for j in V) == 1)
        model.addConstr(gp.quicksum(phi[j][n + 1] for j in V) == 1)
        model.addConstrs(phi[i, j] <= y[i, j] for i in V for j in V if i != j)

        # Non-renewable resource constraints
        model.addConstrs(
            gp.quicksum(r[i][m][k] * x[i, m] for m in M[i]) for i in V for k in self.instance.K_nonrenew)

        # add model and variables to Benders object
        self.master_model = model
        self.master_eta = eta
        self.master_y = y
        self.master_x = x
        self.master_f = f
        self.master_s = s
        self.master_phi = phi


    def solve_master(self):
        """
        Solves master problem and returns optimal objective (which gives LB to robust MRCPSP).

        :return: Optimal objective value of master problem
        :rtype: float
        """

        # solve master_model
        self.master_model.optimize()

        return self.master_model.ObjVal

    def generate_cut(self, master_variables, master_sol, subproblem_objval):
        eta = master_variables['eta']
        y = master_variables['y']
        x = master_variables['x']
        f = master_variables['f']

        self.master_model.addConstr(eta >= (subproblem_objval - self.LB) * gp.quicksum(1/3*(y[i,j] +  ))

    def subproblem(self, t, master_solution):
        """
        Takes job processing modes and resource allocations selected by the master problem and finds longest path
        through the project network, considering uncertain job durations. Produces a feasible solution to the original
        robust MRCPSP problem, and therefore gives an UB.

        :param t: Iteration number
        :type t: int
        :param master_solution: Solution to the master problem
        :type master_solution:
        :return: Objective value of subproblem
        :rtype: float
        """
        # get instance data
        V = self.instance.V
        n = self.instance.n
        d = [self.instance.jobs[j].d for j in V]
        d_bar = [self.instance.jobs[j].d_bar for j in V]
        M = [self.instance.jobs[j].M for j in V]
        # Dummy source and sink all available resource.
        r = [[self.instance.R]] + [self.instance.jobs[j].r for j in self.instance.N] + [[self.instance.R]]
        # compute transitive network from master problem solution
        T_E = [(i, j) for i in V for j in V if master_solution.x[i][j] > 0.99]

        # Create model
        model_name = '{}_benders_sub_{}'.format(self.instance.name, t)
        model = gp.Model(model_name)

        # Add variables to model
        alpha = model.addVars([(i, j) for i in V for j in V], name="alpha", vtype=GRB.BINARY)
        w = model.addVars([(i, j) for i in V for j in V], name="w", vtype=GRB.CONTINUOUS, lb=0)
        xi = model.addVars([i for i in V], name="xi", vtype=GRB.CONTINUOUS, lb=0, ub=1)

        # Add constraints to model
        model.addConstrs(gp.quicksum(alpha[i, n + 1] for i in V if (i, n + 1) in T_E) == 1)
        model.addConstrs(gp.quicksum(alpha[0, i] for i in V if (0, i) in T_E) == 1)
        model.addConstrs(gp.quicksum(alpha[e[0], e[1]] for e in T_E if e[0] == i) - gp.quicksum(
            alpha[e[0], e[1]] for e in T_E if e[1] == i) == 0 for i in V if i != 0 if i != n + 1)
        model.addConstrs(w[e[0], e[1]] <= xi[e[0]] for e in T_E)
        model.addConstrs(w[e[0], e[1]] <= alpha[e[0], e[1]] for e in T_E)
        model.addConstr(gp.quicksum(xi[i] for i in V) <= self.Gamma)

        # solve model
        model.optimize()

        return model.ObjVal
