
import mrcpsp


# read nominal data
instance = mrcpsp.load_nominal_mrcpsp('/home/boldm1/OneDrive/robust-mrcpsp/instances/toy_example_data.txt')
# set d_bar values explicitly
instance.set_dbar_explicitly({0: [0], 1: [0, 3], 2: [8, 2], 3: [1], 4: [3], 5: [0], 6: [0]})


print(instance.jobs[4].d_bar)
