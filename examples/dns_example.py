import sys; sys.path.append("..")

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

import dedaLES

# Basic stuff
n = 128
model = dedaLES.BoussinesqChannelFlow(Lx=2*pi, Ly=2*pi, Lz=1, nx=n, ny=4, nz=n, ν=1, κ=0.1, Fb=0, closure=None)

# Boundary conditions
model.set_bc("nopenetration", "top", "bottom")
model.set_bc("noslip", "top", "bottom")
model.set_bc("noflux", "top", "bottom")

#model.problem.add_bc("left(bz) = 0")
#model.problem.add_bc("right(bz) = Fb")

model.build_solver()
