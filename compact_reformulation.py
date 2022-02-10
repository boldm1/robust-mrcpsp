import gurobipy as gp
from gurobipy import GRB


def compact_reformulation(instance, Gamma, time_limit, print_log=False):
    """
    Solves given robust MRCPSP instance using compact mixed-integer programming reformulation.

    :param instance: Robust MRCPSP instance to solve
    :type instance: Instance
    :param Gamma: Robustness parameter, i.e. max number of jobs that can achieve worst-case durations
            simultaneously. Takes an integer value from 0 to n.
    :type Gamma: int
    :param time_limit: Time limit for Gurobi, in seconds
    :type time_limit: int
    :param print_log: Indicates whether or not to print Gurobi solve log to terminal. Defaults to False.
    :type print_log: bool
    :return: Dictionary containing solution information.
    :rtype: dict
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
    model.setParam('OutputFlag', print_log)
    model.setParam('TimeLimit', time_limit)
    model.setParam('Threads', 4)  # always limit number of threads to 4

    G = range(Gamma + 1)  # levels in auxiliary graph

    # Create variables
    S = model.addVars([(i, g) for i in instance.V for g in G], name='S', vtype=GRB.CONTINUOUS, lb=0)
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
        gp.quicksum(f[i, j, k] for j in instance.V if j != 0) == gp.quicksum(r[i][m][k] * x[i, m] for m in M[i]) for i
        in instance.V if i != instance.n + 1 for k in instance.K_renew)
    model.addConstrs(
        gp.quicksum(f[i, j, k] for i in instance.V if i != instance.n + 1) == gp.quicksum(
            r[j][m][k] * x[j, m] for m in M[j]) for j in instance.V if j != 0 for k in instance.K_renew)
    model.addConstrs(gp.quicksum(x[i, m] for m in M[i]) == 1 for i in instance.V)
    # Non-renewable resource constraints
    model.addConstrs(gp.quicksum(r[i][m][k] * x[i, m] for m in M[i]) <= instance.R[k] for i in instance.V for k in
                     instance.K_nonrenew)

    # solve model
    model.optimize()

    # dict to store solution info
    modes = {i: int(sum(m * x[i, m].X for m in M[i])) for i in instance.V}
    T_E = [(i, j) for i in instance.V for j in instance.V if i != j if y[i, j].X > 0.99]
    flows = {(i, j): [f[i, j, k].X for k in instance.K_renew] for i in instance.V for j in instance.V if
             i != instance.n + 1 if j != 0 if i != j}
    return {'objval': model.ObjVal, 'objbound': model.ObjBound, 'runtime': model.Runtime, 'modes': modes,
            'network': T_E, 'flows': flows}
