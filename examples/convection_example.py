import sys; sys.path.append("..")

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

import dedaLES

# Constant
second = 1.0
minute = 60*second
hour = 60*minute
day = 24*hour

# Domain parameters
nx = 64
ny = 64
nz = 32
Lx = 200
Ly = 200
Lz = 100

# Physical parameters
Q  = -0.1
N2inf = 9.5e-3
h0 = 50
d = 10

# Physical constants
α  = 2.5e-4
g  = 9.81
ρ0 = 1028.1
cP = 3993.0
κ  = 1.43e-7
ν  = 1.05e-7

# Calculated parameters
bz0 = Q*α*g / (cP*ρ0*κ)

# Construct model
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, ν=ν, κ=κ, 
                                      N2inf=N2inf, bz0=bz0, tb=day/2, closure=None)

def smoothstep(z, d): 0.5*(1 + np.tanh(z/d))

# Boundary conditions
model.set_bc("nopenetration", "top", "bottom")
model.set_bc("noslip", "top", "bottom")
model.set_tracer_gradient_bc("b", "top", gradient="bz0*tanh((t-tb)/tb)")
model.set_tracer_gradient_bc("b", "bottom", gradient="Ninf")

model.build_solver()

b0 = N2inf * (model.z + h0) * smoothstep(-model.z-h0, d)
model.set_b(b0)

model.run(initial_dt=10*minute, sim_time=0.5*days, logcadence=100)
