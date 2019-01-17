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
    
# Re = U*L/ν = 
# Rayleigh number. Ra = Δb*L^3 / ν*κ = Δb*L^3*Pr / ν^2
Ra = 2000 #for testing

# Fixed parameters
#nx = ny = kerr_parameters[Ra]['nh']
#nz = kerr_parameters[Ra]['nz']
nx = ny = 8     # for testing
nz = 8          # for testing
Lx = Ly = 6.0               # Horizonal extent
Lz = 1.0                    # Vertical extent
Pr = 0.7                    # Prandtl number
f  = 0.0                    # Coriolis parameter
a  = 1e-3                   # Noise amplitude for initial condition
Δb = 1.0                    # Buoyancy difference

# Calculated parameters
ν = np.sqrt(Δb*Pr*Lz**3/Ra) # Viscosity. ν = sqrt(Pr/Ra) with Lz=Δb=1
κ = ν/Pr                      # Thermal diffusivity 

# Construct model
closure = dedaLES.AnisotropicMinimumDissipation()
#closure = dedaLES.ConstantSmagorinsky()
#closure = None
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, 
                                      ν=ν, κ=κ, Δb=Δb, closure=closure, nu=ν,V=Lx*Ly*Lz)

model.set_bc("no penetration", "top", "bottom")
model.set_bc("no slip", "top", "bottom")
model.problem.add_bc("right(b) = 0")
model.problem.add_bc("left(b) = Δb")


model.build_solver()

# Initial condition: unstable buoyancy grad + random perturbations
noise = dedaLES.random_noise(model.domain)
z = model.z
pert = a * noise * z * (Lz - z)
b0 = (z - pert) / Lz
model.set_fields(b=b0)

model.stop_at(iteration=24000) #sim_time=kerr_parameters[Ra]['tf'])

#Load previous timestepping
#model.solver.load_state('rayleigh_benard_snapshots_DNS/rayleigh_benard_snapshots_DNS_s1/rayleigh_benard_snapshots_DNS_s1_p0.h5', -1)
model.solver.load_state('rayleigh_benard_snapshots_AnisotropicMinimumDissipation/rayleigh_benard_snapshots_AnisotropicMinimumDissipation_s1/rayleigh_benard_snapshots_AnisotropicMinimumDissipation_s1_p0.h5')

# Analysis: Some Timesteps
if closure is None: closure_name = 'DNS'
else:               closure_name = closure.__class__.__name__

analysis = model.solver.evaluator.add_file_handler("rayleigh_benard_snapshots_{:s}".format(closure_name), iter=100, max_writes=100)
analysis.add_system(model.solver.state, layout='g')
# we should consider saving the total state only at the final timestep for memory purposes


# Analysis: All Timesteps
# because Δb = 1, Lz = 1, to obtain the non-dimensional Nusselt number all we must do is divide the vertical heat flux by kappa
nusselt = model.solver.evaluator.add_file_handler('nusselt_AMD', iter=1)
nusselt.add_task("integ(integ(integ(b*w, 'z'), 'x'), 'y') / V / κ", layout='g', name='Nu1')
nusselt.add_task("integ(integ(integ(ε, 'z'), 'x'), 'y')/V / κ", layout='g', name='Nu2')
nusselt.add_task("(integ(integ(integ(bz*bz + dx(b)*dx(b) + dy(b)*dy(b), 'z'), 'x'), 'y')/V-1)", layout='g', name='Nu3')

# CFL
CFL = flow_tools.CFL(
    model.solver, initial_dt=1e-4, cadence=20, safety=1.5, max_change=1.5, min_change=0.5, max_dt=0.05)
CFL.add_velocities(('u', 'v', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(model.solver, cadence=10)
flow.add_property("sqrt(u*u + v*v + w*w) / nu", name='Re_domain')
flow.add_property("sqrt(ε) / nu", name='Re_dissipation')
flow.add_property("ε", name="dissipation")

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
            logger.info('Max domain Re = %f' %flow.max("Re_domain"))
            logger.info('Max dissipation Re = %f' %flow.max("Re_dissipation"))
            logger.info('Average epsilon = %f' %flow.volume_average("dissipation"))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %model.solver.iteration)
    logger.info('Sim end time: %f' %model.solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/60/60*model.domain.dist.comm_cart.size))
