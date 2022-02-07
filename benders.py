import time
import gurobipy as gp
from gurobipy import GRB


class Benders:
    """
    Class to solve robust MRCPSP instance using Bender' style decomposition.
    """

    def __init__(self, instance, Gamma):
        """
        Initialises Benders object with instance to solve and robustness parameter. Algorithm parameters and master and
        sub-problem models and variables are also initialised.

        :param instance: Robust MRCPSP instance to solve
        :type instance: Instance
        :param Gamma: Robustness parameter, i.e. max number of jobs that can achieve worst-case durations
            simultaneously. Takes an integer value from 0 to n (= number of non-dummy jobs in instance).
        :type Gamma: int

        """
        self.instance = instance
        self.Gamma = Gamma

        self.LB = float('-inf')  # LB on objective function of robust MRCPSP instance
        self.UB = float('inf')  # UB on objective function of robust MRCPSP instance
        self.cuts = []  # cuts applied to master problem

        # Master problem model and variables. Populated in create_master_problem and updated when solved.
        self.master_model = None
        self.master_eta = None
        self.master_y = None
        self.master_x = None
        self.master_f = None
        self.master_s = None
        self.master_phi = None
        # Sub-problem model and objective value. Updated each time sub-problem is generated and solved.
        self.sub_model = None
        self.sub_objval = None

    def solve(self, time_limit, write_master=False, print_log=False):
        """
        Runs Benders' decomposition algorithm to solve robust MRCPSP instance. Terminates upon finding an optimal
        solution or reaching the specified time limit.

        :param time_limit: Max running time of algorithm in seconds. Current UB is returned if time_limit is reached.
        :type time_limit: int
        :param write_master: Indicates whether or not to write solution of master model to output file each iteration.
            Defaults to False.
        :type write_master: bool
        :param print_log: Indicates whether or not to print Gurobi solve log to terminal. Defaults to False.
        :type print_log: bool
        :return: Dictionary containing objective value of solution, overall solve time, number of iterations.
        :rtype: dict
        """
        # Benders' decomposition algorithm
        t = 0  # iteration number
        start_time = time.time()
        self.create_master_problem(print_log)  # create basic master problem with no optimality cuts
        while self.LB < self.UB:
            if print_log is True:
                print("Iteration {}".format(t))
                print("------------")
                print("Solving master problem...")
            # solve master problem and update LB
            self.master_model.optimize()
            if write_master is True:
                self.master_model.write('{}_iter{}.sol'.format(self.master_model.ModelName, t))
            master_objval = self.master_model.ObjVal
            if master_objval > self.LB:
                self.LB = self.master_model.ObjVal
            # solve sub-problem and update UB
            self.get_subproblem(t, print_log)
            if print_log is True:
                print("Solving sub-problem...")
            self.sub_model.optimize()
            self.sub_objval = self.sub_model.ObjVal
            if self.sub_objval < self.UB:
                self.UB = self.sub_objval
            if print_log is True:
                print("LB = {}, UB = {}".format(self.LB, self.UB))
                print("-----------------------------------------------------------------------------\n")
            # check termination conditions
            if self.UB - self.LB < 1e-6:  # allow for numerical imprecision
                solve_time = time.time() - start_time
                return {'objval': self.UB, 'solve_time': solve_time, 'n_iterations': t}  # optimal solution
            elif time.time() - start_time > time_limit:
                return {'objval': self.UB, 'solve_time': time_limit, 'n_iterations': t}  # time-limit reached
            # add cut and go to next iteration
            else:
                self.add_cut(t)
                t += 1

    def create_master_problem(self, print_log=False):
        """
        Creates Gurobi model to represent basic master problem without optimality cuts. The master problem selects job
        processing modes and resource allocations to minimise project duration under nominal job processing times. The
        solution to the master problem provides a LB to optimal robust solution.

        :param print_log: Indicates whether or not to print Gurobi solve log to terminal when the model is solved.
            Defaults to False.
        :type print_log: bool
        """
        # get instance data
        V = self.instance.V
        n = self.instance.n
        d = [self.instance.jobs[j].d for j in V]  # nominal durations
        d_bar = [self.instance.jobs[j].d_bar for j in V]  # max durational deviation
        M = [self.instance.jobs[j].M for j in V]  # available processing modes
        # Dummy source and sink require all available resource
        r = [[self.instance.R]] + [self.instance.jobs[j].r for j in self.instance.N] + [[self.instance.R]]

        # Create model
        model_name = '{}_benders_master'.format(self.instance.name)
        model = gp.Model(model_name)
        model.setParam('OutputFlag', print_log)

        # Create variables
        eta = model.addVar(name="eta")  # variable to represent objective value
        y = model.addVars([(i, j) for i in V for j in V], name="y",
                          vtype=GRB.BINARY)  # y_ij = 1 if i -> j in transitive project precedence network
        x = model.addVars([(i, m) for i in V for m in M[i]], name="x",
                          vtype=GRB.BINARY)  # x_im = 1 if i is executed in mode m
        f = model.addVars([(i, j, k) for i in V for j in V for k in self.instance.K_renew], name="f",
                          vtype=GRB.CONTINUOUS, lb=0)  # f_ijk = flow from i to j of renewable resource k
        s = model.addVars([i for i in V], name="s", vtype=GRB.INTEGER)  # start time of each job
        # longest path variables
        phi = model.addVars([(i, j) for i in V for j in V], name="phi", vtype=GRB.CONTINUOUS, lb=0)

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
        # Non-renewable resource constraints
        model.addConstrs(gp.quicksum(r[i][m][k] * x[i, m] for m in M[i]) <= self.instance.R[k] for i in V for k in
                         self.instance.K_nonrenew)

        # UB on min makespan, i.e. largest possible gap between s[j] and s[i].
        N = sum(max(d[i][m] + d_bar[i][m] for m in M[i]) for i in V)
        model.addConstrs(
            s[j] - s[i] >= d[i][m] - N * (1 - y[i, j]) - N * (1 - x[i, m]) for i in V for j in V if i != j for m in
            M[i])
        model.addConstr(eta >= s[n])

        # longest path with minimal job durations
        model.addConstr(eta >= gp.quicksum(min(d[i][m] for m in M[i]) * phi[i, j] for i in V for j in V))
        model.addConstrs(
            gp.quicksum(phi[j, i] for j in V) == gp.quicksum(phi[i, j] for j in V) for i in V if i != 0 if i != n + 1)
        model.addConstr(gp.quicksum(phi[0, j] for j in V) == 1)
        model.addConstr(gp.quicksum(phi[j, n + 1] for j in V) == 1)
        model.addConstrs(phi[i, j] <= y[i, j] for i in V for j in V if i != j)

        # update model and variable parameters
        self.master_model = model
        self.master_eta = eta
        self.master_y = y
        self.master_x = x
        self.master_f = f
        self.master_s = s
        self.master_phi = phi

    def add_cut(self, t):

        """
        Adds optimality cut to master problem based on solution of the sub-problem.

        :param t Iteration number
        :type t: int
        """
        # instance data
        V = self.instance.V
        M = [self.instance.jobs[j].M for j in V]

        # get transitive network from master problem
        T_E = [(i, j) for i in V for j in V if self.master_y[i, j].X > 0.99]
        # get mode selection from master problem
        master_M = [int(sum(m * self.master_x[i, m].X for m in M[i])) for i in V]

        # some number larger than number of extra precedences in optimal master solution
        N = self.instance.n * self.instance.n
        cut = self.master_model.addConstr(self.master_eta >= (self.sub_objval - self.LB) * gp.quicksum(
            1 / 3 * (self.master_y[e[0], e[1]] + self.master_x[e[0], master_M[e[0]]] +
                     self.master_x[e[1], master_M[e[1]]]) -
            N * (3 - self.master_y[e[0], e[1]] - self.master_x[e[0], master_M[e[0]]] -
                 self.master_x[e[1], master_M[e[1]]]) for e in T_E)
                                          - (self.sub_objval - self.LB) * (len(T_E) - 1) + self.LB,
                                          name="cut_{}".format(t))
        # add cut to list of cuts
        self.cuts.append(cut)

    def get_subproblem(self, t, print_log=False):
        """
        Takes job processing modes and resource allocations selected by the master problem and builds model to find
        longest path through the project network, considering uncertain job durations. When solved, this model produces
        a feasible solution to the original robust MRCPSP problem and therefore gives an UB.

        :param t: Iteration number
        :type t: int
        :param print_log: Indicates whether or not to print Gurobi solve log to terminal when the model is solved.
            Defaults to False.
        :type print_log: bool
        """
        # get instance data
        V = self.instance.V
        n = self.instance.n
        d = [self.instance.jobs[j].d for j in V]
        d_bar = [self.instance.jobs[j].d_bar for j in V]
        M = [self.instance.jobs[j].M for j in V]

        # get transitive network from master problem
        T_E = [(i, j) for i in V for j in V if self.master_y[i, j].X > 0.99]
        # get mode selection from master problem
        master_M = [int(sum(m * self.master_x[i, m].X for m in M[i])) for i in V]

        # Create model
        model_name = '{}_benders_sub_{}'.format(self.instance.name, t)
        model = gp.Model(model_name)
        model.setParam('OutputFlag', print_log)

        # Add variables to model
        alpha = model.addVars([(i, j) for i in V for j in V], name="alpha", vtype=GRB.BINARY)
        w = model.addVars([(i, j) for i in V for j in V], name="w", vtype=GRB.CONTINUOUS, lb=0)
        xi = model.addVars([i for i in V], name="xi", vtype=GRB.CONTINUOUS, lb=0, ub=1)

        # Add constraints to model
        model.addConstr(gp.quicksum(alpha[i, n + 1] for i in V if (i, n + 1) in T_E) == 1)
        model.addConstr(gp.quicksum(alpha[0, i] for i in V if (0, i) in T_E) == 1)
        model.addConstrs(gp.quicksum(alpha[e[0], e[1]] for e in T_E if e[0] == i) - gp.quicksum(
            alpha[e[0], e[1]] for e in T_E if e[1] == i) == 0 for i in V if i != 0 if i != n + 1)
        model.addConstrs(w[e[0], e[1]] <= xi[e[0]] for e in T_E)
        model.addConstrs(w[e[0], e[1]] <= alpha[e[0], e[1]] for e in T_E)
        model.addConstr(gp.quicksum(xi[i] for i in V) <= self.Gamma)

        # set objective function
        model.setObjective(gp.quicksum(
            d[e[0]][master_M[e[0]]] * alpha[e[0], e[1]] + d_bar[e[0]][master_M[e[0]]] * w[e[0], e[1]] for e in T_E),
            GRB.MAXIMIZE)

        # update sub-problem model parameter
        self.sub_model = model
