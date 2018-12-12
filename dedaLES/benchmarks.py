import sys; sys.path.append("..")

import time, logging
import numpy as np
import matplotlib.pyplot as plt

from numpy import pi
from dedalus.extras import flow_tools

import dedaLES


def build_rayleigh_benard_benchmark(nx, ny, nz, closure=None):
    Lx = 25.0       # Domain horizontal extent
    Ly = 5.0        # Domain horizontal extent
    Lz = 1.0        # Domain vertical extent

    # Parameters
    Pr = 1.0        # Prandtl number
    f  = 0.0        # Coriolis parameter
    κ  = 1.0        # Thermal diffusivity 
    ε  = 0.8        # Perturbation above criticality
    a  = 1e-3       # Noise amplitude for initial condition

    # Constants
    Ra_critical = 1707.762
    Ra = Ra_critical + ε
    ν  = Pr*κ                   # Viscosity 
    Bz = -Ra*Pr                 # Unstable buoyancy gradient

    # Construct model
    model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, 
                                          nx=nx, ny=ny, nz=nz, 
                                          ν=ν, κ=κ, B0=Bz*Lz, closure=closure)

    model.set_bc("nopenetration", "top", "bottom")
    model.set_bc("freeslip", "top", "bottom")
    model.problem.add_bc("right(b) = B0")
    model.problem.add_bc("left(b) = 0")

    model.build_solver()

    # Random perturbations, initialized globally for same results in parallel
    noise = dedaLES.random_noise(model.domain)

    # Linear background + perturbations damped at walls
    z = model.z
    pert = a * noise * z * (Lz - z)
    b0 = Bz*(z - pert)
    model.set_b(b0)

    return model


def run_rayleigh_benard_benchmark(model, iterations=100):

    # Integration parameters
    model.solver.stop_sim_time = 100
    model.solver.stop_wall_time = np.Inf
    model.solver.stop_iteration = iterations

    dt = 1e-4
    while model.solver.ok:
        model.solver.step(dt)
