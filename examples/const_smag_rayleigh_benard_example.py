import sys; sys.path.append("..")

import time, logging
import numpy as np
import matplotlib.pyplot as plt

from numpy import pi
from dedalus.extras import flow_tools

import dedaLES

logger = logging.getLogger(__name__)

# Domain parameters
nx = 64         # Horizontal resolution
ny = 16         # Horizontal resolution
nz = 16         # Vertical resolution
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
                                      ν=ν, κ=κ, B0=Bz*Lz, closure=dedaLES.ConstantSmagorinsky())

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

# Integration parameters
model.solver.stop_sim_time = 100
model.solver.stop_wall_time = 60 * 60.
model.solver.stop_iteration = np.inf

# Analysis
snap = model.solver.evaluator.add_file_handler('snapshots', sim_dt=0.2,
                                               max_writes=10)
snap.add_task("interp(b, z=0)", scales=1, name='b midplane')
snap.add_task("interp(u, z=0)", scales=1, name='u midplane')
snap.add_task("interp(v, z=0)", scales=1, name='v midplane')
snap.add_task("interp(w, z=0)", scales=1, name='w midplane')

# CFL
CFL = flow_tools.CFL(model.solver, initial_dt=1e-4, cadence=5,
                     safety=1.5, max_change=1.5, min_change=0.5, max_dt=0.05)
CFL.add_velocities(('u', 'v', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(model.solver, cadence=10)
flow.add_property("sqrt(u*u + v*v + w*w) / ν", name='Re')

# Main loop
try:
    logger.info('Starting loop')
    start_run_time = time.time()
    while model.solver.ok:
        dt = CFL.compute_dt()
        model.solver.step(dt)
        if (model.solver.iteration-1) % 100 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(
                        model.solver.iteration, model.solver.sim_time, dt))
            logger.info('Max Re = %f' %flow.max('Re'))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %model.solver.iteration)
    logger.info('Sim end time: %f' %model.solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/60/60*model.domain.dist.comm_cart.size))
