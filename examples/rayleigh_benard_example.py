"""
This script reproduces results from

Robert M Kerr, "Rayleigh number scaling in numerical convection", 
Journal of Fluid Mechanics (1996)
"""

import sys; sys.path.append("..")

import time, logging
import numpy as np
import matplotlib.pyplot as plt

from numpy import pi
from dedalus.extras import flow_tools

import dedaLES

logger = logging.getLogger(__name__)

resolutions = {
    1: {nh:  64, nz: 32},
    2: {nh:  96, nz: 48},
    3: {nh: 128, nz: 48},
    4: {nh: 192, nz: 64},
    5: {nh: 288, nz: 96} 
}

experiments = {
    5000: { **resolutions[1],
            Ra : 5e4,
            t0 : 24.0,
            tf : 44.0 },
    10000: { **resolutions[1],
            Ra : 1e5,
            t0 : 24.0,
            tf : 44.0 },
    20000: { **resolutions[1]
            Ra : 2e5,
            t0 : 24.0,
            tf : 44.0 },
    40000: { **resolutions[3],
            Ra : 4e5,
            t0 : 26.0,
            tf : 34.0 },
    50000: { **resolutions[2],
            Ra : 4e5,
            t0 : 27.0,
            tf : 40.0 },
} 
            
    
# Domain parameters
nx = 64         # Horizontal resolution
ny = 64         # Horizontal resolution
nz = 16         # Vertical resolution
Lx = 25.0       # Domain horizontal extent
Ly = 5.0        # Domain horizontal extent
Lz = 1.0        # Domain vertical extent

# Parameters
Pr = 1.0        # Prandtl number
Ra = 2e4        # Rayleigh number
f  = 0.0        # Coriolis parameter
a  = 1e-3       # Noise amplitude for initial condition

κ  = Pr         # Thermal diffusivity 
Bz = -Ra*Pr     # Unstable buoyancy gradient
ν  = 1/Re       # Viscosity 

# Construct model
closure = None #dedaLES.AnisotropicMinimumDissipation()
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, 
                                      ν=ν, κ=κ, B0=Bz*Lz, closure=closure)

model.set_bc("no penetration", "top", "bottom")
model.set_bc("free slip", "top", "bottom")
model.problem.add_bc("right(b) = B0")
model.problem.add_bc("left(b) = 0")

model.build_solver()

# Iniital condition: unstable buoyancy grad + random perturbations
noise = dedaLES.random_noise(model.domain)
z = model.z
pert = a * noise * z * (Lz - z)
b0 = Bz*(z - pert)
model.set_fields(b=b0)

# Integration parameters
model.solver.stop_sim_time = 100
model.solver.stop_wall_time = 60 * 60.
model.solver.stop_iteration = np.inf

# Analysis
if closure is None:
    closure_name = 'DNS'
else:
    closure_name = closure.__class__.__name__

snap = model.solver.evaluator.add_file_handler(
        "snapshots_rayleigh_benard_{:s}".format(closure_name), sim_dt=0.2, max_writes=10)
snap.add_task("interp(b, z=0)", scales=1, name='b midplane')
snap.add_task("interp(u, z=0)", scales=1, name='u midplane')
snap.add_task("interp(v, z=0)", scales=1, name='v midplane')
snap.add_task("interp(w, z=0)", scales=1, name='w midplane')
snap.add_task("integ(b, 'z')", name='b integral x4', scales=4)

# CFL
CFL = flow_tools.CFL(model.solver, initial_dt=1e-6, cadence=5, 
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
        if (model.solver.iteration-1) % 10 == 0:
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
