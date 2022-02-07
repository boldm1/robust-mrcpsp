
import mrcpsp
from compact_reformulation import compact_reformulation
from benders import Benders

# read nominal data
instance = mrcpsp.load_nominal_mrcpsp('/home/boldm1/OneDrive/robust-mrcpsp/instances/toy_example_data.txt')
# set d_bar values explicitly
instance.set_dbar_explicitly({0: [0], 1: [0, 3], 2: [8, 2], 3: [1], 4: [3], 5: [0], 6: [0]})
#instance.set_dbar_uncertainty_level(0.7)

Gamma = 3
time_limit = 20

print("Compact reformulation")
print("---------------------")
objval, solve_time = compact_reformulation(instance, Gamma, time_limit, write_sol=False)
print("objval:", objval)
print("solve_time:", solve_time)

print("\nBenders'")
print("--------")
benders_sol = Benders(instance, Gamma).solve(time_limit, print_log=True)
print("objval:", benders_sol['objval'])
print("solve_time:", benders_sol['solve_time'])