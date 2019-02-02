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

# Parameters
Ra = 1e5    # Rayleigh number. Ra = Δb*L^3 / ν*κ = Δb*L^3*Pr / ν^2
Pr = 0.7    # Prandtl number
a  = 1e-1   # Noise amplitude for initial condition
Δb = 1.0    # Buoyancy difference

# Calculated parameters
ν = np.sqrt(Pr/Ra) # Viscosity. ν = sqrt(Pr/Ra) with Lz=Δb=1
κ = ν/Pr           # Thermal diffusivity 

# Construct model
closure = None #dedaLES.AnisotropicMinimumDissipation() #dedaLES.ConstantSmagorinsky()
model = dedaLES.BoussinesqChannelFlow(Lx=1, Ly=1, Lz=6, nx=64, ny=64, nz=32, ν=ν, κ=κ, Δb=1, closure=closure)

model.set_bc("no penetration", "top", "bottom")
model.set_bc("no slip", "top", "bottom")
model.problem.add_bc("right(b) = 0")
model.problem.add_bc("left(b) = Δb")

model.build_solver()

# Set initial condition: unstable buoyancy grad + random perturbations
noise = a * dedaLES.random_noise(model.domain) * model.z * (model.Lz - model.z) / model.Lz**2,
model.set_fields(
    u = noise,
    v = noise,
    w = noise,
    b = model.z/model.Lz - noise
)
model.stop_at(sim_time=100)

# Analysis
if closure is None: closure_name = 'DNS'
else:               closure_name = closure.__class__.__name__

simulation_name = "rayleigh_benard_snapshots_{:s}".format(closure_name)
analysis = model.solver.evaluator.add_file_handler(simulation_name, iter=100, max_writes=100)
analysis.add_system(model.solver.state, layout='g')

stats = flow_tools.GlobalFlowProperty(model.solver, cadence=10)
stats.add_property("sqrt(u*u + v*v + w*w) / ν", name="Re")
stats.add_property("1 + w*b/κ", name="Nu_wb")
stats.add_property("1 + (ε + ε_sgs)/κ", name="Nu_ε")
stats.add_property("(χ + χ_sgs)/κ", name="Nu_χ")

# CFL
CFL = flow_tools.CFL(model.solver, initial_dt=1e-4, cadence=10, safety=0.5, max_change=1.5)
CFL.add_velocities(('u', 'v', 'w'))

# Main loop
try:
    logger.info('Starting loop')
    start_run_time = time.time()
    log_time = time.time()  
    while model.solver.ok:
        dt = CFL.compute_dt()
        model.solver.step(dt)
        if (model.solver.iteration-1) % 10 == 0:
            compute_time = time.time() - log_time
            log_time = time.time()
            logger.info(
                "i: {:d}, t: {:.2e}, t_wall: {:.1f} s, dt: {:.2e}, max Re: {:.0f}, Nu_wb: {:.2f}, Nu_ε: {:.2f}, Nu_χ: {:.2f}".format(
                model.solver.iteration, model.solver.sim_time, compute_time, dt, stats.max("Re"), 
                stats.volume_average("Nu_wb"), stats.volume_average("Nu_ε"), stats.volume_average("Nu_χ")))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise

finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %model.solver.iteration)
    logger.info('Sim end time: %f' %model.solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/60/60*model.domain.dist.comm_cart.size))
