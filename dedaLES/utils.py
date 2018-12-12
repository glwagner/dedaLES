import numpy as np
from mpi4py import MPI

def random_noise(domain):
    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand = np.random.RandomState(seed=23)
    return rand.standard_normal(gshape)[slices]

def mpiprint(msg):
    if MPI.COMM_WORLD.Get_rank() is 0:
        print(msg)
