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

# From Table 2 in Kerr (1996):
resolutions = {
    1: {'nh':  64, 'nz': 32},
    2: {'nh':  96, 'nz': 48},
    3: {'nh': 128, 'nz': 48},
    4: {'nh': 192, 'nz': 64},
    5: {'nh': 288, 'nz': 96} 
}

kerr_parameters = {
    5000: { 'nh' : 64,
            'nz' : 32,
            't0' : 24.0,
            'tf' : 44.0 },
   10000: { 'nh' : 64,
            'nz' : 32,
            't0' : 24.0,
            'tf' : 44.0 },
   20000: { 'nh' : 64,
            'nz' : 32,
            't0' : 24.0,
            'tf' : 44.0 }}
#   40000: { **resolutions[3],
#            't0' : 26.0,
#            'tf' : 34.0 },
#   50000: { **resolutions[2],
#            't0' : 27.0,
#            'tf' : 40.0 },
#  100000: { **resolutions[3],
#            't0' : 26.0,
#            'tf' : 1000.0 },
#  250000: { **resolutions[3],
#            't0' : 26.0,
#            'tf' : 36.0 },
#  500000: { **resolutions[4],
#            't0' : 49.0,
#            'tf' : 64.0 },
# 1000000: { **resolutions[4],
#            't0' : 28.0,
#            'tf' : 37.0 },
# 2000000: { **resolutions[5],
#            't0' : 24.0,
#            'tf' : 140.0 }
#} 
#
    
# Rayleigh number
Ra = 10000

# Fixed parameters
nx = ny = kerr_parameters[Ra]['nh']
nz = kerr_parameters[Ra]['nz']
Lx = Ly = 12.0          # Horizonal extent
Lz = 2.0                # Vertical extent
Pr = 0.7                # Prandtl number
f  = 0.0                # Coriolis parameter
a  = 1e-3               # Noise amplitude for initial condition

# Calculated parameters
κ  = Pr                 # Thermal diffusivity 
Bz = -Ra*Pr             # Unstable buoyancy gradient
ν  = -Lz**4*Bz/(Ra*κ)    # Viscosity 

# Construct model
closure = None #dedaLES.AnisotropicMinimumDissipation()
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, zbottom=-Lz/2, nx=nx, ny=ny, nz=nz, 
                                      ν=ν, κ=κ, B0=Bz*Lz, closure=closure)

model.set_bc("no penetration", "top", "bottom")
model.set_bc("no slip", "top", "bottom")
model.problem.add_bc("right(b) = B0")
model.problem.add_bc("left(b) = 0")

model.build_solver()

# Initial condition: unstable buoyancy grad + random perturbations
noise = dedaLES.random_noise(model.domain)
z = model.z
pert = a * noise * z * (Lz - z)
b0 = Bz*(z - pert)
model.set_fields(b=b0)

model.stop_at(iteration=1000) #sim_time=kerr_parameters[Ra]['tf'])

# Analysis
if closure is None: closure_name = 'DNS'
else:               closure_name = closure.__class__.__name__
    
analysis = model.solver.evaluator.add_file_handler(
    "rayleigh_benard_snapshots_{:s}".format(closure_name), iter=100, max_writes=100)

analysis.add_system(model.solver.state, layout='g')

# CFL
CFL = flow_tools.CFL(
    model.solver, initial_dt=1e-4, cadence=20, safety=1.5, max_change=1.5, min_change=0.5, max_dt=0.05)
CFL.add_velocities(('u', 'v', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(model.solver, cadence=10)
flow.add_property("sqrt(u*u + v*v + w*w) / ν", name='Re')
flow.add_property("ux*ux + uy*uy + uz*uz + vx*vx + vy*vy + vz*vz + wx*wx + wy*wy + wz*wz", name="epsilon")

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
            logger.info('Average epsilon = %f' %flow.volume_average('epsilon'))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %model.solver.iteration)
    logger.info('Sim end time: %f' %model.solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/60/60*model.domain.dist.comm_cart.size))
