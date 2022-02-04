import gurobipy as gp
from gurobipy import GRB


def compact_reformulation(instance, Gamma, time_limit):
    """
    Solves given robust MRCPSP instance using compact mixed-integer programming reformulation.

    :param instance: Robust MRCPSP instance to solve
    :type instance: Instance
    :param Gamma: Value of Gamma to protect against
    :type Gamma: int
    :return: Solution object
    :rtype: Solution
    """
    # get instance data
    d = [instance.jobs[j].d for j in instance.V]  # d[j][m] = duration of j in mode m
    d_bar = [instance.jobs[j].d_bar for j in instance.V]  # d_bar[j][m] = max durational deviation of j in mode m
    M = [instance.jobs[j].M for j in instance.V]  # M[j] = list of modes of j
    # r[j][m][k] = amount of resource k required by j in mode m. Dummy source and sink all available resource.
    r = [instance.R] + [instance.jobs[j].r for j in instance.N] + [instance.R]

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
    x = model.addVars([(i, m) for i in instance.V for m in instance.jobs[i].modes], name="x",
                      vtype=GRB.BINARY)  # x_im = 1 if i is executed in mode m
    f = model.addVars([(i, j, k) for i in instance.V for j in instance.V for k in instance.K_renew], name="f",
                      vtype=GRB.CONTINUOUS, lb=0)  # f_ijk = flow from i to j of renewable resource k

    # set model objective
    model.setObjective(S[instance.n + 1, Gamma], GRB.MINIMIZE)

    # add constraints
    BigM = sum(max(d[i][m] + d_bar[i][m] for m in M[i]) for i in instance.V)  # largest possible makespan
    BigR = [[[min(max(r[i][m][k] for m in M[i]), max(r[j][m][k] for m in M[j])) for k in instance.K_renew] for j in
             instance.V] for i in instance.V]

    model.addConstr(S[0, 0] == 0)
    model.addConstrs(
        S[j, g] - S[i, g] >= d[i][m] * x[i, m] - BigM * (1 - y[i, j]) for i in instance.V for j in instance.V for m in
        M[i] for g in G)
    model.addConstrs(
        S[j, g + 1] - S[i, g] >= (d[i][m] + d_bar[i][m]) * x[i, m] - BigM * (1 - y[i, j]) for i in instance.V for j in
        instance.V for m in
        M[i] for g in G[:-1])

    model.addConstrs(y[e[0], e[1]] == 1 for e in instance.E)
    model.addConstr(y[instance.n + 1, instance.n + 1] == 1)

    model.addConstrs(
        f[i, j, k] <= BigR[i][j][k] * y[i, j] for i in instance.V for j in instance.V if i != instance.n + 1 if j != 0
        for k in instance.K_renew)
    model.addConstrs(
        gp.quicksum(f[i, j, k] for j in instance.V) == gp.quicksum(r[i][m][k] * x[i, m] for m in M[i]) for i in
        instance.V if i != instance.n + 1 for k in instance.K)
    model.addConstrs(
        gp.quicksum(f[i, j, k] for i in instance.V) == gp.quicksum(r[j][m][k] * x[j, m] for m in M[j]) for j in
        instance.V if j != 0 for k in instance.K)
    model.addConstrs(gp.quicksum(x[i, m] for m in M[i]) == 1 for i in instance.V)

    # solve model
    model.optimize()
    model.write('{}.sol'.format(model_name))
