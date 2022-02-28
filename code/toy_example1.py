
import mrcpsp
from compact_reformulation import CompactRefomulation
from benders import Benders

"""
Testing code by solving toy example from Balouka and Cohen (2021).
"""

# read nominal data
file_path = '/home/boldm1/OneDrive/robust-mrcpsp/instances/mrcpsp_toy_example_data.txt'
instance = mrcpsp.load_nominal_mrcpsp(file_path)
# set d_bar values explicitly
instance.set_dbar_explicitly({0: [0], 1: [0, 3], 2: [8, 2], 3: [1], 4: [3], 5: [0], 6: [0]})

# solve parameters
Gamma = 3
time_limit = 20

print(f"\nSolving toy example given in {file_path}...")
print("\nJob data [r, d, d_bar]:")
print("-------------------------")
for i in instance.V:
    print("job {}:".format(i),
          [[instance.jobs[i].r[m][0], instance.jobs[i].d[m], instance.jobs[i].d_bar[m]] for m in instance.jobs[i].M])

print("\nCompact reformulation")
print("---------------------")
compact_sol = CompactRefomulation(instance, Gamma, time_limit).solve(print_log=False)
print("objval:", compact_sol['objval'])
print("runtime:", compact_sol['runtime'])
print("modes:", compact_sol['modes'])
print("network:", compact_sol['network'])
print("resource flows:", compact_sol['flows'])

print("\nBenders'")
print("--------")
benders_sol = Benders(instance, Gamma, time_limit).solve(print_log=False)
print("objval:", benders_sol['objval'])
print("runtime:", benders_sol['runtime'])
print("modes:", benders_sol['modes'])
print("network:", compact_sol['network'])
print("resource flows:", benders_sol['flows'])
