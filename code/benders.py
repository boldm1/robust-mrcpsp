import time
import gurobipy as gp
from gurobipy import GRB


class Benders:
    """
    Class to solve uncertain MRCPSP instance using Bender' style decomposition.
    """

    def __init__(self, instance, Gamma, time_limit):
        """
        Initialises Benders object with instance to solve and algorithm parameters.

        :param instance: Robust MRCPSP instance to solve
        :type instance: Instance
        :param Gamma: Robustness parameter, i.e. max number of jobs that can achieve worst-case durations
            simultaneously. Takes an integer value from 0 to n (= number of non-dummy jobs in instance).
        :type Gamma: int
        :param time_limit: Max running time of algorithm in seconds. Current UB is returned if time_limit is reached.
        :type time_limit: int
        """
        self.instance = instance
        self.Gamma = Gamma
        self.time_limit = time_limit

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
        self.master_S = None

        # Sub-problem model, objective value and parameters. Updated each time sub-problem is generated and solved.
        self.sub_model = None
        self.sub_objval = None
        self.sub_alpha = None
        self.sub_w = None
        self.sub_xi = None

    def solve(self, num_threads=4, print_log=False, master_type='compact'):
        """
        Runs Benders' decomposition algorithm to solve robust MRCPSP instance. Terminates upon finding an optimal
        solution or reaching the specified time limit.

        :param num_threads: Number of threads to use when solving. Defaults to 4.
        :type num_threads: int
        :param print_log: Indicates whether or not to print Gurobi solve log to terminal. Defaults to False.
        :type print_log: bool
        :param master_type: Master problem formulation to use. Either 'balouka' or 'compact'. Defaults to 'compact'.
        :type str
        :return: Dictionary containing solution information.
        :rtype: dict
        """
        # check master type is recognised
        if master_type not in ['balouka', 'compact']:
            raise Exception("Specified master problem formulation type is not recognised. Must be either 'balouka' or "
                            "'compact'. Please correct this and try again!")

        # Benders' decomposition algorithm
        t = 1  # iteration number
        iteration_times = []  # List to store time of each iteration. Average iteration time is reported in solution.
        best_sol = {'modes': {}, 'network': [], 'flows': {}}  # dict to store x, y, f variable values of best solution
        start_time = time.perf_counter()
        if master_type == 'compact':
            self.create_compact_master_problem(num_threads, print_log)
        elif master_type == 'balouka':
            self.create_balouka_master_problem(num_threads, print_log)
        while self.LB < self.UB:
            iteration_start = time.perf_counter()  # start timing iteration
            # set master problem time limit to remaining algorithm time limit to prevent algorithm getting stuck in the
            # master problem
            time_remaining = self.time_limit - (iteration_start - start_time)
            # if time limit reached, terminate and return solution
            if time_remaining <= 0:
                average_iteration_time = sum(iteration_times) / t
                return {'objval': self.UB, 'objbound': self.LB, 'runtime': self.time_limit, 'n_iterations': t,
                        'avg_iteration': average_iteration_time,
                        'modes': best_sol['modes'], 'network': best_sol['network'], 'flows': best_sol['flows']}
            else:
                if print_log is True:
                    print("Iteration {}".format(t))
                    print("------------")
                    print("Solving master problem...")
                self.master_model.setParam('TimeLimit', time_remaining)
                # solve master problem
                self.master_model.optimize()
                if print_log is True:
                    self.print_master_sol()
                if self.master_model.status == 9:  # if time_limit reached whilst in master problem, terminate
                    # compute average iteration time, setting to +inf if terminated whilst in first iteration
                    if len(iteration_times) == 0:
                        average_iteration_time = float('inf')
                    else:
                        average_iteration_time = sum(iteration_times) / t
                    return {'objval': self.UB, 'objbound': self.LB, 'runtime': self.time_limit, 'n_iterations': t - 1,
                            'avg_iteration': average_iteration_time,
                            'modes': best_sol['modes'], 'network': best_sol['network'], 'flows': best_sol['flows']}
                else:
                    # if master problem provides tighter LB, update
                    master_objval = self.master_model.ObjVal
                    if master_objval > self.LB:
                        self.LB = self.master_model.ObjVal
                    # solve sub-problem and update UB
                    self.get_subproblem(t, num_threads, print_log)
                    if print_log is True:
                        print("\nSolving sub-problem...")
                    self.sub_model.optimize()
                    self.sub_objval = self.sub_model.ObjVal
                    if print_log is True:
                        self.print_subproblem_sol()
                    # if sub-problem provides tighter UB, update and record solution
                    if self.sub_objval < self.UB:
                        self.UB = self.sub_objval
                        # get solution info and save to dict
                        modes = {i: int(sum(m * self.master_x[i, m].X for m in self.instance.jobs[i].M)) for i in
                                 self.instance.V}
                        T_E = [(i, j) for i in self.instance.V for j in self.instance.V if i != j if
                               self.master_y[i, j].X > 0.99]
                        flows = {(i, j): [self.master_f[i, j, k].X for k in self.instance.K_renew] for i in
                                 self.instance.V for j in self.instance.V if i != self.instance.n + 1 if j != 0 if
                                 i != j}
                        best_sol = {'modes': modes, 'network': T_E, 'flows': flows}
                    if print_log is True:
                        print("LB = {}, UB = {}".format(self.LB, self.UB))
                        print("-----------------------------------------------------------------------------\n")
                    # check if LB == UB
                    if self.UB - self.LB < 1e-6:  # allow for numerical imprecision
                        # optimal solution
                        runtime = time.perf_counter() - start_time
                        iteration_times.append(time.perf_counter() - iteration_start)  # add final iteration time
                        average_iteration_time = sum(iteration_times) / t
                        return {'objval': self.UB, 'objbound': self.LB, 'runtime': runtime, 'n_iterations': t,
                                'avg_iteration': average_iteration_time,
                                'modes': best_sol['modes'], 'network': best_sol['network'], 'flows': best_sol['flows']}
                    # add cut and go to next iteration
                    else:
                        self.add_cut(t)
                        iteration_times.append(time.perf_counter() - iteration_start)
                        t += 1

    def print_master_sol(self):
        """
        Prints mode selection, transitive network and resource flows from master solution to the terminal.
        """
        modes = {i: int(sum(m * self.master_x[i, m].X for m in self.instance.jobs[i].M)) for i in
                 self.instance.V}
        T_E = [(i, j) for i in self.instance.V for j in self.instance.V if i != j if
               self.master_y[i, j].X > 0.99]
        flows = {(i, j): [self.master_f[i, j, k].X for k in self.instance.K_renew] for i in
                 self.instance.V for j in self.instance.V if i != self.instance.n + 1 if j != 0 if
                 i != j}
        print("master objval:", self.master_model.ObjVal)
        print("modes:", modes)
        print("network:", T_E)
        print("resource flows:", flows)

    def print_subproblem_sol(self):
        """
        Prints longest path found in subproblem to the terminal.
        """
        T_E = [(i, j) for i in self.instance.V for j in self.instance.V if i != j if self.master_y[i, j].X > 0.99]
        longest_path = [(e[0], e[1]) for e in T_E if self.sub_alpha[e[0], e[1]].X > 0.99]
        print("longest path:", longest_path)

    def create_compact_master_problem(self, num_threads=4, print_log=False):
        """
        Creates Gurobi model to represent basic master problem without optimality cuts. The master problem selects job
        processing modes and resource allocations to minimise project duration under nominal job processing times. The
        solution to the master problem provides a LB to optimal robust solution. Uses compact reformulation model
        without uncertainty.

        :param num_threads: Number of threads to use when solving. Defaults to 4.
        :type num_threads: int
        :param print_log: Indicates whether or not to print Gurobi solve log to terminal when the model is solved.
            Defaults to False.
        :type print_log: bool
        """
        # get instance data
        V = self.instance.V
        n = self.instance.n
        E = self.instance.E
        d = [self.instance.jobs[j].d for j in V]  # nominal durations
        d_bar = [self.instance.jobs[j].d_bar for j in V]  # max durational deviation
        M = [self.instance.jobs[j].M for j in V]  # available processing modes
        # Dummy source and sink require all available resource
        r = [[self.instance.R]] + [self.instance.jobs[j].r for j in self.instance.N] + [[self.instance.R]]

        # Create model
        model_name = '{}_benders_master'.format(self.instance.name)
        model = gp.Model(model_name)
        model.setParam('OutputFlag', print_log)
        model.setParam('Threads', num_threads)

        # Create variables
        eta = model.addVar(name="eta")  # variable to represent objective value
        S = model.addVars([i for i in V], name='S', vtype=GRB.CONTINUOUS, lb=0)
        y = model.addVars([(i, j) for i in V for j in V], name="y",
                          vtype=GRB.BINARY)  # y_ij = 1 if i -> j in transitive project precedence network
        x = model.addVars([(i, m) for i in V for m in M[i]], name="x",
                          vtype=GRB.BINARY)  # x_im = 1 if i is executed in mode m
        f = model.addVars([(i, j, k) for i in V for j in V for k in self.instance.K_renew], name="f",
                          vtype=GRB.CONTINUOUS, lb=0)  # f_ijk = flow from i to j of renewable resource k

        # set model objective
        model.setObjective(eta, GRB.MINIMIZE)

        # Precedence constraints
        model.addConstr(S[0] == 0)
        BigM = sum(max(d[i][m] + d_bar[i][m] for m in M[i]) for i in V)  # largest possible makespan
        model.addConstrs(
            S[j] - S[i] >= d[i][m] * x[i, m] - BigM * (1 - y[i, j]) for i in V for j in V for m in M[i])
        model.addConstr(eta >= S[n + 1])
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

        # Transitivity constraints to improve model performance.
        model.addConstrs(y[i, j] + y[j, i] <= 1 for i in V for j in V if i != n + 1 if j != n + 1)
        model.addConstrs(y[i, j] >= y[i, l] + y[l, j] - 1 for i in V for j in V for l in V)

        # update model and variable parameters
        self.master_model = model
        self.master_eta = eta
        self.master_S = S
        self.master_y = y
        self.master_x = x
        self.master_f = f

    def create_balouka_master_problem(self, num_threads=4, print_log=False):
        """
        Creates Gurobi model to represent basic master problem without optimality cuts. The master problem selects job
        processing modes and resource allocations to minimise project duration under nominal job processing times. The
        solution to the master problem provides a LB to optimal robust solution. Uses formulation presented in Balouka &
        Cohen (2021).

        :param num_threads: Number of threads to use when solving. Defaults to 4.
        :type num_threads: int
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
        model.setParam('Threads', num_threads)

        # Create variables
        eta = model.addVar(name="eta")  # variable to represent objective value
        y = model.addVars([(i, j) for i in V for j in V if i != j], name="y",
                          vtype=GRB.BINARY)  # y_ij = 1 if i -> j in transitive project precedence network
        x = model.addVars([(i, m) for i in V for m in M[i]], name="x",
                          vtype=GRB.BINARY)  # x_im = 1 if i is executed in mode m
        f = model.addVars(
            [(i, j, k) for i in V for j in V if i != j if i != n + 1 if j != 0 for k in self.instance.K_renew],
            name="f", vtype=GRB.CONTINUOUS, lb=0)  # f_ijk = flow from i to j of renewable resource k
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
            V if i != n + 1 for j in V if j != 0 if i != j for m_i in M[i] for m_j in M[j] for k in
            self.instance.K_renew)
        model.addConstrs(
            gp.quicksum(f[i, j, k] for j in V if j != i if j != 0) == gp.quicksum(r[i][m][k] * x[i, m] for m in M[i])
            for i in V if i != n + 1 for k in self.instance.K_renew)
        model.addConstrs(
            gp.quicksum(f[i, j, k] for i in V if i != j if i != n + 1) == gp.quicksum(
                r[j][m][k] * x[j, m] for m in M[j]) for j in V if j != 0 for k in self.instance.K_renew)
        model.addConstrs(gp.quicksum(x[i, m] for m in M[i]) == 1 for i in V)
        # Non-renewable resource constraints
        model.addConstrs(gp.quicksum(r[i][m][k] * x[i, m] for m in M[i]) <= self.instance.R[k] for i in V for k in
                         self.instance.K_nonrenew)

        # UB on min makespan, i.e. largest possible gap between s[j] and s[i].
        N = sum(max(d[i][m] + d_bar[i][m] for m in M[i]) for i in V)
        model.addConstrs(
            s[j] - s[i] >= d[i][m] - N * (1 - y[i, j]) - N * (1 - x[i, m]) for i in V for j in V if i != j for m in
            M[i])
        model.addConstr(eta >= s[n + 1])

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

        # transitive network from master problem
        T_E = [(i, j) for i in V for j in V if i != j if self.master_y[i, j].X > 0.99]

        # get mode selection from master problem
        master_M = [int(sum(m * self.master_x[i, m].X for m in M[i])) for i in V]
        # get longest path from sub-problem
        sub_path = [(e[0], e[1]) for e in T_E if self.sub_alpha[e[0], e[1]].X > 0.99]

        # some number larger than number of extra precedences in optimal master solution
        cut = self.master_model.addConstr(self.master_eta >= (self.sub_objval - self.LB) * gp.quicksum(
            1 / 3 * (self.master_y[e[0], e[1]] + self.master_x[e[0], master_M[e[0]]] +
                     self.master_x[e[1], master_M[e[1]]]) -
            (3 - self.master_y[e[0], e[1]] - self.master_x[e[0], master_M[e[0]]] - self.master_x[e[1], master_M[e[1]]])
            for e in sub_path) - (self.sub_objval - self.LB) * (
                                                  len(sub_path) - 1) + self.LB, name="cut_{}".format(t))
        # add cut to list of cuts
        self.cuts.append(cut)

    def get_subproblem(self, t, num_threads=4, print_log=False):
        """
        Takes job processing modes and resource allocations selected by the master problem and builds model to find
        longest path through the project network, considering uncertain job durations. When solved, this model produces
        a feasible solution to the original robust MRCPSP problem and therefore gives an UB.

        :param t: Iteration number
        :type t: int
        :param num_threads: Number of threads to use when solving. Defaults to 4.
        :type num_threads: int
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
        T_E = [(i, j) for i in V for j in V if i != j if self.master_y[i, j].X > 0.99]
        # get mode selection from master problem
        master_M = [int(sum(m * self.master_x[i, m].X for m in M[i])) for i in V]

        # Create model
        model_name = '{}_benders_sub_{}'.format(self.instance.name, t)
        model = gp.Model(model_name)
        model.setParam('OutputFlag', print_log)
        model.setParam('Threads', num_threads)

        # Add variables to model
        alpha = model.addVars([(e[0], e[1]) for e in T_E], name="alpha", vtype=GRB.BINARY)
        w = model.addVars([(e[0], e[1]) for e in T_E], name="w", vtype=GRB.CONTINUOUS, lb=0)
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

        # update sub-problem model and variable parameters
        self.sub_model = model
        self.sub_alpha = alpha
        self.sub_w = w
        self.sub_xi = xi
