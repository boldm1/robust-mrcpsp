
from mrcpsp import load_mrcpsp


instance = load_mrcpsp('/home/boldm1/OneDrive/robust-mrcpsp/instances/j10.mm/j102_2.mm')
print(instance.jobs[1].LS)
