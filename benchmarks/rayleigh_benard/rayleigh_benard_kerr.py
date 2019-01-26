"""
This script reproduces results from

Robert M Kerr, "Rayleigh number scaling in numerical convection", 
Journal of Fluid Mechanics (1996)
"""

import os, sys; sys.path.append(os.path.join("..", "..",))
import time, logging
import numpy as np
import matplotlib.pyplot as plt

from mpi4py import MPI
from numpy import pi
from dedalus.extras import flow_tools

import dedaLES

debug = False
logger = logging.getLogger(__name__)

# From Table 2 in Kerr (1996):
kerr_key = {
    '1': 50000,
    '2': 100000,
    '3': 200000,
    '4': 400000,
    '5': 500000,
    '6': 1000000,
    '7': 2500000,
    '8': 5000000
}

kerr_parameters = {
    50000: {          'nh' : 64,
                      'nz' : 32,
                      't0' : 24.0,
                      'tf' : 44.0,
                      'dt' : 0.0025,
             'spinup_time' : 100,
                   'iters' : 40000,
            },
   100000: {          'nh' : 64,
                      'nz' : 32,
                      't0' : 24.0,
                      'tf' : 44.0,
                      'dt' : 0.0025,
             'spinup_time' : 100,
                   'iters' : 40000,
            },
   200000: {          'nh' : 64,
                      'nz' : 32,
                      't0' : 24.0,
                      'tf' : 44.0,
                      'dt' : 0.0025,
             'spinup_time' : 100,
                   'iters' : 40000,
            },
   400000: {          'nh' : 96,
                      'nz' : 48,
                      't0' : 26.0,
                      'tf' : 44.0,
                      'dt' : 0.0025,
             'spinup_time' : 100,
                   'iters' : 40000,
            },
   500000: {          'nh' : 96,
                      'nz' : 48,
                      't0' : 27.0,
                      'tf' : 40.0,
                      'dt' : 0.0025,
             'spinup_time' : 100,
                   'iters' : 40000,
            },
  1000000: {          'nh' : 128,
                      'nz' : 48,
                      't0' : 26.0,
                      'tf' : 1000.0, 
                      'dt' : 0.0025,
             'spinup_time' : 100,
                   'iters' : 40000,
            },
  2500000: {          'nh' : 128,
                      'nz' : 48,
                      't0' : 26.0,
                      'tf' : 36.0,
                      'dt' : 0.0025,
             'spinup_time' : 100,
                   'iters' : 40000,
            },
  5000000: {          'nh' : 192,
                      'nz' : 64,
                      't0' : 49.0,
                      'tf' : 64.0,
                      'dt' : 0.0025,
             'spinup_time' : 100,
                   'iters' : 40000,
            },
}


def identifier(model):
    return "nh{:d}_nz{:d}_dt{:.0f}_Ra{:d}".format(model.nx, model.nz, 
        10000*kerr_parameters[model.Ra]['dt'], model.Ra)
        
# Rayleigh number. Ra = Δb*L^3 / ν*κ = Δb*L^3*Pr / ν^2
Ra = kerr_key[sys.argv[1]]

# Parameters
dt = kerr_parameters[Ra]['dt']
Lx = Ly = 6.0                 # Horizonal extent
Lz = 1.0                      # Vertical extent
Pr = 0.7                      # Prandtl number
f  = 0.0                      # Coriolis parameter
a  = 1e-3                     # Noise amplitude for initial condition
Δb = 1.0                      # Buoyancy difference
ν  = np.sqrt(Δb*Pr*Lz**3/Ra)  # Viscosity. ν = sqrt(Pr/Ra) with Lz=Δb=1
κ  = ν/Pr                     # Thermal diffusivity 
log_cadence = 100

if debug:
    nx = ny = nz = 8
    spinup_time = 10*dt
else:
    nx = ny = kerr_parameters[Ra]['nh']
    nz = kerr_parameters[Ra]['nz']
    spinup_time = kerr_parameters[Ra]['spinup_time']

# Construct model
closure = None
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, ν=ν, κ=κ, Δb=Δb, 
                                      closure=closure, nu=ν, V=Lx*Ly*Lz, H=Lz, Ra=Ra)
model.set_bc("no penetration", "top", "bottom")
model.set_bc("no slip", "top", "bottom")
model.problem.add_bc("right(b) = 0")
model.problem.add_bc("left(b) = Δb")

model.build_solver(timestepper='SBDF3')

logger.info("Simulation identifier: {}\n".format(identifier(model)))

# Analysis: Some Timesteps
if closure is None: closure_name = 'DNS'
else:               closure_name = closure.__class__.__name__

# Set initial condition: unstable buoyancy grad + random perturbations
noise = dedaLES.random_noise(model.domain)
z = model.z
pert = a * noise * z * (Lz - z)
u0 = pert / Lz
b0 = (z - pert) / Lz
model.set_fields(u=u0, b=b0)
model.stop_at(sim_time=spinup_time)

# CFL
CFL = flow_tools.CFL(model.solver, initial_dt=dt, cadence=10, safety=0.5)
CFL.add_velocities(('u', 'v', 'w'))

# Flow properties
spinup = flow_tools.GlobalFlowProperty(model.solver, cadence=log_cadence)
spinup.add_property("sqrt(u*u + v*v + w*w) / nu", name='Re_domain')

# Three Nusselt numbers
spinup.add_property("1 + w*b/κ", name="Nu1")
spinup.add_property("1 + ε/κ", name="Nu2")
spinup.add_property("χ/κ", name="Nu3")

# Main loop
try:
    logger.info("Starting loop. Ra: {}, closure: {}".format(Ra, closure_name))
    start_run_time = time.time()
    while model.solver.ok:
        dt = CFL.compute_dt()
        model.solver.step(dt)
        if (model.solver.iteration-1) % log_cadence == 0:
            logger.info("iter: {: 7d}, dt: {:e}, t: {:.2f}, max Re: {:.0f}, Nu1: {:.2f}, Nu2: {:.2f}, Nu3: {:.2f}".format(
                            model.solver.iteration, dt, model.solver.sim_time, spinup.max("Re_domain"), 
                            spinup.volume_average("Nu1"), spinup.volume_average("Nu2"), spinup.volume_average("Nu3")))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise

finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %model.solver.iteration)
    logger.info('Sim end time: %f' %model.solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/3600 * model.domain.dist.comm_cart.size))

logger.info("\n\n Spinup complete for Ra: {}, closure: {}".format(Ra, closure_name))

# Reset
model.solver.iteration = 0
model.stop_at(iteration=kerr_parameters[Ra]['iters'])
dt = kerr_parameters[Ra]['dt']

if debug: model.stop_at(iteration=10)

# Flow properties
stats = flow_tools.GlobalFlowProperty(model.solver, cadence=1)

# Three Nusselt numbers
stats.add_property("1 + w*b/κ", name="Nu1")
stats.add_property("1 + ε/κ", name="Nu2")
stats.add_property("χ/κ", name="Nu3")

t = [0.0]
Nu1 = [spinup.volume_average("Nu1")] 
Nu2 = [spinup.volume_average("Nu2")]
Nu3 = [spinup.volume_average("Nu3")]

nusselt_filename = "nusselt_{}_{}".format(identifier(model), closure_name)

def save_arrays(filename, **kwargs):
    for k, a in kwargs.items():
        kwargs[k] = np.array(a)
    np.savez(nusselt_filename, **kwargs)

try: 
    logger.info("Gathering statistics with dt = {:e}...".format(dt))
    while model.solver.ok:
        model.solver.step(dt)

        Nu1.append(stats.volume_average("Nu1"))
        Nu2.append(stats.volume_average("Nu2"))
        Nu3.append(stats.volume_average("Nu3"))
        t.append(model.solver.sim_time)
    
        logger.info("iter: {: 7d}, Nu1: {:.2f}, Nu2: {:.2f}, Nu3: {:.2f}".format(model.solver.iteration, Nu1[-1], Nu2[-1], Nu3[-1]))
                            
except:
    if MPI.COMM_WORLD.Get_rank() is 0: save_arrays(nusselt_filename, Nu1=Nu1, Nu2=Nu2, Nu3=Nu3, t=t) 
    logger.error('Exception raised, triggering end of main loop.')
    raise

finally:
    if MPI.COMM_WORLD.Get_rank() is 0: save_arrays(nusselt_filename, Nu1=Nu1, Nu2=Nu2, Nu3=Nu3, t=t) 
    logger.info("Statistics gathering complete.")
