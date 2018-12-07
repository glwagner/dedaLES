import sys; sys.path.append("..")

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

import dedaLES

# Basic stuff
closure = dedaLES.ConstantSmagorinsky()
model = dedaLES.BoussinesqChannelFlow(Lx=2*pi, Ly=2*pi, Lz=1, nx=32, ny=32, nz=32, ν=1, κ=0.1, closure=closure)
model.set_default_bcs()
model.build_solver()

# Initial condition
Re = 100
H = model.Lz/10
B = Re**2 / H
z0 = -model.Lz/2
τ = np.sqrt(H/B)

b0 = -B * np.exp(-(model.z-z0)**2/(2*H**2)) # unstable blob
model.set_b(b0)

model.run(initial_dt=1e-2*τ, iterations=10, logcadence=1)
