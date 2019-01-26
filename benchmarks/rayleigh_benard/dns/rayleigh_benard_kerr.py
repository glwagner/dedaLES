"""
This script reproduces results from

Robert M Kerr, "Rayleigh number scaling in numerical convection", 
Journal of Fluid Mechanics (1996)
"""

import os, sys; sys.path.append(os.path.join("..", "..", ".."))
import time, logging
import numpy as np
import matplotlib.pyplot as plt

from numpy import pi
from dedalus.extras import flow_tools

import dedaLES

logger = logging.getLogger(__name__)

# From Table 2 in Kerr (1996):
kerr_parameters = {
    50000: {  'nh' : 64,
              'nz' : 32,
              't0' : 24.0,
              'tf' : 44.0,
              'dt' : 0.025,
            'iters': 40000,
            },
   100000: {  'nh' : 64,
              'nz' : 32,
              't0' : 24.0,
              'tf' : 44.0,
              'dt' : 0.025,
            'iters': 40000,
            },
   200000: {  'nh' : 64,
              'nz' : 32,
              't0' : 24.0,
              'tf' : 44.0,
              'dt' : 0.025,
            'iters': 40000,
            },
   400000: {  'nh' : 96,
              'nz' : 48,
              't0' : 26.0,
              'tf' : 44.0,
              'dt' : 0.0025,
            'iters': 40000,
            },
   500000: {  'nh' : 96,
              'nz' : 48,
              't0' : 27.0,
              'tf' : 40.0,
              'dt' : 0.0025,
            'iters': 40000,
            },
  1000000: {  'nh' : 128,
              'nz' : 48,
              't0' : 26.0,
              'tf' : 1000.0, 
              'dt' : 0.0025,
            'iters': 40000,
            },
  2500000: {  'nh' : 128,
              'nz' : 48,
              't0' : 26.0,
              'tf' : 36.0,
              'dt' : 0.0025,
            'iters': 40000,
            },
  5000000: {  'nh' : 192,
              'nz' : 64,
              't0' : 49.0,
              'tf' : 64.0,
              'dt' : 0.0025,
            'iters': 40000,
            },
}

def identifier(model):
    return "nh{:d}_nz{:d}_dt{:.0f}_Ra{:d}".format(model.nx, model.nz, 10000*model.dt, model.Ra)
        
# Rayleigh number. Ra = Δb*L^3 / ν*κ = Δb*L^3*Pr / ν^2
Ra = 500000
load_previous = False # load from previous simulation

# Parameters
nx = ny = kerr_parameters[Ra]['nh']
nz = kerr_parameters[Ra]['nz']
dt = kerr_parameters[Ra]['dt']
#nx = ny = nz = 8
Lx = Ly = 6.0                 # Horizonal extent
Lz = 1.0                      # Vertical extent
Pr = 0.7                      # Prandtl number
f  = 0.0                      # Coriolis parameter
a  = 1e-3                     # Noise amplitude for initial condition
Δb = 1.0                      # Buoyancy difference
ν  = np.sqrt(Δb*Pr*Lz**3/Ra)  # Viscosity. ν = sqrt(Pr/Ra) with Lz=Δb=1
κ  = ν/Pr                     # Thermal diffusivity 

# Construct model
closure = None
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, ν=ν, κ=κ, Δb=Δb, 
                                      closure=closure, nu=ν, V=Lx*Ly*Lz, H=Lz, Ra=Ra, dt=dt)
model.set_bc("no penetration", "top", "bottom")
model.set_bc("no slip", "top", "bottom")
model.problem.add_bc("right(b) = 0")
model.problem.add_bc("left(b) = Δb")

model.build_solver(timestepper='SBDF3')

logger.info("Simulation identifier: {}\n".format(identifier(model)))

# Analysis: Some Timesteps
if closure is None: closure_name = 'DNS'
else:               closure_name = closure.__class__.__name__

snapshot_directory = "snapshots_{}_{}".format(identifier(model), closure_name)
nusselt_directory = "nusselt_{}_{}".format(identifier(model), closure_name)

if load_previous:
    # Note: we need to move the h5 file if we want to load from a previous state
    filepath = os.path.join(snapshot_directory, "{}_s1.h5".format(snapshot_directory))
    absolute_filepath = os.path.abspath(filepath)
    dedaLES.mpiprint("Loading from {}".format(absolute_filepath))
    model.solver.load_state(absolute_filepath)
    model.stop_at(iteration=solver.iteration + kerr_parameters[Ra]['iters']) # sim_time=kerr_parameters[Ra]['tf'])
else: # Set initial condition: unstable buoyancy grad + random perturbations
    noise = dedaLES.random_noise(model.domain)
    z = model.z
    pert = a * noise * z * (Lz - z)
    u0 = pert / Lz
    b0 = (z - pert) / Lz
    model.set_fields(u=u0, b=b0)
    model.stop_at(iteration=kerr_parameters[Ra]['iters']) # sim_time=kerr_parameters[Ra]['tf'])

snapshots = model.solver.evaluator.add_file_handler(snapshot_directory, iter=10000, max_writes=30)
snapshots.add_system(model.solver.state, layout='g')

# Analysis: All Timesteps
# because Δb = 1, Lz = 1, to obtain the non-dimensional Nusselt number all we must do is divide the vertical heat flux by kappa
nusselt = model.solver.evaluator.add_file_handler(nusselt_directory, iter=1)
nusselt.add_task("integ(integ(integ(w*b, 'z'), 'x'), 'y') / (V*κ)", layout='g', name='Nu1')
nusselt.add_task("integ(integ(integ(ε, 'z'), 'x'), 'y') / (V*κ)", layout='g', name='Nu2')
nusselt.add_task("(integ(integ(integ(bx*bx + by*by + bz*bz, 'z'), 'x'), 'y')/V - 1)", layout='g', name='Nu3')

if closure is not None:
    nusselt.add_task("integ(integ(integ(ε_sgs, 'z'), 'x'), 'y')/V ", layout='g', name='sgs_dis')
    nusselt.add_task("integ(integ(integ(χ_sgs, 'z'), 'x'), 'y')/V ", layout='g', name='sgs_xi')

# Flow properties
flow = flow_tools.GlobalFlowProperty(model.solver, cadence=10)
flow.add_property("sqrt(u*u + v*v + w*w) / nu", name='Re_domain')
flow.add_property("sqrt(ε)", name='Re_dissipation')
flow.add_property("ε", name="dissipation")
flow.add_property("w*b/(κ*V)", name="Nu1")
flow.add_property("ε/(κ*V)", name="Nu2")
flow.add_property("χ/V - 1", name="Nu3")


# Main loop
try:
    logger.info("Starting loop. Ra: {}, closure: {}".format(Ra, closure_name))
    start_run_time = time.time()
    while model.solver.ok:
        model.solver.step(model.dt)
        if (model.solver.iteration-1) % 10 == 0:
            logger.info("\nIteration: {:d}, Time: {:e}\n".format(model.solver.iteration, model.solver.sim_time))

            logger.info("     Max domain Re: {:.4f}".format(flow.max("Re_domain")))
            logger.info("Max dissipation Re: {:.4f}".format(flow.max("Re_dissipation")))
            logger.info("   Average epsilon: {:.4f}".format(flow.volume_average("dissipation")))
            logger.info("               Nu1: {:.6f}".format(flow.volume_average("Nu1")))
            logger.info("               Nu2: {:.6f}".format(flow.volume_average("Nu2")))
            logger.info("               Nu3: {:.6f}".format(flow.volume_average("Nu3")))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise

finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %model.solver.iteration)
    logger.info('Sim end time: %f' %model.solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/3600 * model.domain.dist.comm_cart.size))

logger.info("\n\nRa: {}, closure: {}".format(Ra, closure_name))
