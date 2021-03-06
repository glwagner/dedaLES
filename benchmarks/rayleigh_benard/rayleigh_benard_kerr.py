"""
This script reproduces results from

Robert M Kerr, "Rayleigh number scaling in numerical convection", 
Journal of Fluid Mechanics (1996)
"""

import os, sys; sys.path.append(os.getenv('DEDALES', os.path.join("..", "..")))

import time, logging
import numpy as np
import matplotlib.pyplot as plt

from mpi4py import MPI
from numpy import pi
from dedalus.extras import flow_tools

import dedaLES

# Partly from Table 2 in Kerr (1996), partly new:
kerr_parameters = {
     '1': {'Ra':    50000, 'nh':  64, 'nz':  32, 'spinup_time': 100, 'equil_time': 40, 'stats_time': 200},  
     '2': {'Ra':   100000, 'nh':  64, 'nz':  32, 'spinup_time': 100, 'equil_time': 40, 'stats_time': 200},
     '3': {'Ra':   200000, 'nh':  64, 'nz':  32, 'spinup_time': 100, 'equil_time': 40, 'stats_time': 200},
     '4': {'Ra':   400000, 'nh':  96, 'nz':  48, 'spinup_time': 100, 'equil_time': 20, 'stats_time': 100},
     '5': {'Ra':   500000, 'nh':  96, 'nz':  48, 'spinup_time': 100, 'equil_time': 20, 'stats_time':  20},
     '6': {'Ra':  1000000, 'nh': 128, 'nz':  48, 'spinup_time': 100, 'equil_time': 10, 'stats_time':  10},
     '7': {'Ra':  2500000, 'nh': 128, 'nz':  48, 'spinup_time': 100, 'equil_time': 10, 'stats_time':  10},
     '8': {'Ra':  5000000, 'nh': 192, 'nz':  64, 'spinup_time': 100, 'equil_time': 10, 'stats_time':  10},
     '9': {'Ra': 10000000, 'nh': 192, 'nz':  64, 'spinup_time': 100, 'equil_time': 10, 'stats_time':  10},
    '10': {'Ra': 20000000, 'nh': 288, 'nz':  96, 'spinup_time':  30, 'equil_time': 20, 'stats_time':  10},
    '11': {'Ra': 20000000, 'nh': 576, 'nz': 192, 'spinup_time':  30, 'equil_time': 10, 'stats_time':  10},
}


def identifier(m): return "nh{:d}_nz{:d}_Ra{:d}".format(m.nx, m.nz, m.Ra)

def set_rayleigh_benard_bcs(model):
    model.set_bc("no penetration", "top", "bottom")
    model.set_bc("no slip", "top", "bottom")
    model.problem.add_bc("right(b) = 0")
    model.problem.add_bc("left(b) = Δb")

def add_reynolds_number(flow):
    # Assumes L=1
    flow.add_property("sqrt(u*u + v*v + w*w) / ν", name='Re')

def add_nusselt_numbers(flow):
    flow.add_property("1 + w*b/κ", name="Nu_wb")
    flow.add_property("1 + (ε + ε_sgs)/κ", name="Nu_ε")
    flow.add_property("(χ + χ_sgs)/κ", name="Nu_χ")

logger = logging.getLogger(__name__)

# Setup
if len(sys.argv) == 1 or sys.argv[1] is 'debug':
    debug = True
    run = '1'
else:
    debug = False
    run = sys.argv[1]

if len(sys.argv) > 2:
    closure = getattr(dedaLES, sys.argv[2])()
    closure_name = sys.argv[2]
else:
    closure = None
    closure_name = 'DNS'

# Parameters
Ra = kerr_parameters[run]['Ra']         # Rayleigh number. Ra = Δb*L^3 / ν*κ = Δb*L^3*Pr / ν^2
nx = ny = kerr_parameters[run]['nh']    # Horizontal resolution
nz = kerr_parameters[run]['nz']         # Vertical resolution
Lx = Ly = 6.0                           # Horizonal extent
Lz = 1.0                                # Vertical extent
Pr = 0.7                                # Prandtl number
Δb = 1.0                                # Buoyancy difference
a  = 1e-1                               # Noise amplitude for initial condition
dt0 = 1e-4                              # Initial time-step for spinup and equilibration: 
                                        #   small frac of non-dim time, t=1.
ν  = np.sqrt(Δb*Pr*Lz**3/Ra)            # Viscosity. ν = sqrt(Pr/Ra) with Lz=Δb=1
κ  = ν/Pr                               # Thermal diffusivity 

CFL_cadence = 10
CFL_safety = 0.5

spinup_time = kerr_parameters[run]['spinup_time']
spinup_coarsening_factor = 2
spinup_log_cadence = 100

equil_time = kerr_parameters[run]['equil_time']
equil_log_cadence = 100

stats_dt_frac = 0.25 # fraction by which to reduce time-step for stats-gathering
stats_time = kerr_parameters[run]['stats_time']

if debug: # Overwrite a few things for a test run
    nx = ny = nz = 8
    CFL_cadence = np.inf
    CFL_safety = 0.1
    spinup_time = 10*dt0
    spinup_log_cadence = 1
    equil_log_cadence = 1
    equil_time = 10*dt0
    stats_time = 10*dt0
    
# Messaging
logger.info("""\n\n
    *** Rayleigh-Benard convection ***\n
    Ra: {}
    nh: {}
    nz: {}
    closure: {}\n\n""".format(Ra, nx, nz, closure_name)
)

## 
## Spinup!
##

nx_spinup = int(nx/spinup_coarsening_factor)
ny_spinup = int(ny/spinup_coarsening_factor)
nz_spinup = int(nz/spinup_coarsening_factor)

spinup_model = dedaLES.BoussinesqChannelFlow(
    Lx=Lx, Ly=Ly, Lz=Lz, nx=nx_spinup, ny=ny_spinup, nz=nz_spinup, ν=ν, κ=κ, Δb=Δb, closure=closure, Ra=Ra)
set_rayleigh_benard_bcs(spinup_model)
spinup_model.build_solver(timestepper='SBDF3')

# Set initial condition: unstable buoyancy grad + random perturbations
noise = a * dedaLES.random_noise(spinup_model.domain) * spinup_model.z * (Lz - spinup_model.z) / Lz**2,
spinup_model.set_fields(
    u = noise,
    v = noise,
    w = noise,
    b = spinup_model.z/Lz - noise
)
spinup_model.stop_at(sim_time=spinup_time)

# Make CFL
spinup_CFL = flow_tools.CFL(spinup_model.solver, initial_dt=dt0, cadence=10, max_change=1.5, safety=0.5)
spinup_CFL.add_velocities(('u', 'v', 'w'))

# Flow properties: Reynolds number and three Nusselt numbers
spinup_stats = flow_tools.GlobalFlowProperty(spinup_model.solver, cadence=spinup_log_cadence)
add_reynolds_number(spinup_stats)
add_nusselt_numbers(spinup_stats)

# Spinup loop with coarse model
try:
    logger.info("Starting coarse spinup. Ra: {}, closure: {}".format(Ra, closure_name))
    start_run_time = time.time()
    log_time = time.time()  
    while spinup_model.solver.ok:
        dt = spinup_CFL.compute_dt()
        spinup_model.solver.step(dt)
        if (spinup_model.solver.iteration-1) % spinup_log_cadence == 0:
            compute_time = time.time() - log_time
            log_time = time.time()
            logger.info(
                "i: {:d}, t: {:.2e}, t_wall: {:.1f} s, dt: {:.2e}, max Re: {:.0f}, Nu_wb: {:.2f}, Nu_ε: {:.2f}, Nu_χ: {:.2f}".format(
                spinup_model.solver.iteration, spinup_model.solver.sim_time, compute_time, dt, spinup_stats.max("Re"), 
                spinup_stats.volume_average("Nu_wb"), spinup_stats.volume_average("Nu_ε"), spinup_stats.volume_average("Nu_χ")))
except:
    logger.error('Exception raised, triggering end of spinup loop.')
    raise

spinup_wall_time = time.time() - start_run_time
logger.info("""\n
    *** Spinup complete. Starting equilibration. ***\n
    Spinup time: {:f} sec ({:.2f} hr)
    Spinup time: {:f} cpu-hr\n\n""".format(
        spinup_wall_time, spinup_wall_time/3500, spinup_wall_time/3600 * spinup_model.domain.dist.comm_cart.size)
)

##
## Equilibration!
##

# Store final state of coarse model
spinup_model.u.set_scales(spinup_coarsening_factor)
spinup_model.v.set_scales(spinup_coarsening_factor)
spinup_model.w.set_scales(spinup_coarsening_factor)
spinup_model.b.set_scales(spinup_coarsening_factor)

u0 = np.copy(spinup_model.u['g']) 
v0 = np.copy(spinup_model.v['g'])
w0 = np.copy(spinup_model.w['g'])
b0 = np.copy(spinup_model.b['g'])

# Deallocate!
spinup_model = None
spinup_stats = None
spinup_CFL   = None

# Set up full-resolution model
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, ν=ν, κ=κ, Δb=Δb, closure=closure, Ra=Ra)
set_rayleigh_benard_bcs(model)
model.build_solver(timestepper='SBDF3')

model.set_fields(u=u0, v=v0, w=w0, b=b0)
model.stop_at(sim_time=equil_time)

equil_CFL = flow_tools.CFL(model.solver, initial_dt=dt0, cadence=CFL_cadence, safety=CFL_safety)
equil_CFL.add_velocities(('u', 'v', 'w'))

# Flow properties
equil_stats = flow_tools.GlobalFlowProperty(model.solver, cadence=equil_log_cadence)
add_reynolds_number(equil_stats)
add_nusselt_numbers(equil_stats)

# Equilibration run
try:
    logger.info("Starting equilibration. Ra: {}, closure: {}".format(Ra, closure_name))
    start_run_time = time.time()
    log_time = time.time()  
    while model.solver.ok:
        dt = equil_CFL.compute_dt()
        model.solver.step(dt)

        if (model.solver.iteration-1) % equil_log_cadence == 0:
            compute_time = time.time() - log_time
            log_time = time.time()
            logger.info(
                "i: {:d}, t: {:.2e}, t_wall: {:.1f} s, dt: {:.2e}, max Re: {:.0f}, Nu_wb: {:.2f}, Nu_ε: {:.2f}, Nu_χ: {:.2f}".format(
                model.solver.iteration, model.solver.sim_time, compute_time, dt, equil_stats.max("Re"), 
                equil_stats.volume_average("Nu_wb"), equil_stats.volume_average("Nu_ε"), equil_stats.volume_average("Nu_χ")))

except:
    logger.error('Exception raised, triggering end of equilibration loop.')
    raise

equil_wall_time = time.time() - start_run_time
logger.info("""\n
    *** Equilibration complete. Starting statistics gathering.\n
    Equilibration time: {:f} sec ({:.2f} hr)
    Equilibration + spinup time: {:f} sec ({:.2f} hr)
    Equilibration time: {:f} cpu-hr\n\n""".format(equil_wall_time, equil_wall_time/3600, 
        spinup_wall_time+equil_wall_time, (spinup_wall_time+equil_wall_time)/3600, 
        equil_wall_time/3600 * model.domain.dist.comm_cart.size)
)


##
## Statistics!
##

# Reset
model.stop_at(sim_time=model.solver.sim_time+stats_time)
dt *= stats_dt_frac # Fix time-step at safe value

# Flow properties
stats = flow_tools.GlobalFlowProperty(model.solver, cadence=1)
add_nusselt_numbers(stats)

t, Nu_wb, Nu_ε, Nu_χ = [], [], [], []
nusselt_filename = "nusselt_{}_{}".format(identifier(model), closure_name)

def save_arrays(filename, **kwargs):
    if MPI.COMM_WORLD.Get_rank() is 0: 
        for k, a in kwargs.items():
            kwargs[k] = np.array(a)
        np.savez(filename, **kwargs)

try: 
    logger.info("Gathering statistics with dt = {:e}...".format(dt))
    start_run_time = time.time()
    while model.solver.ok:
        model.solver.step(dt)

        Nu_wb.append(stats.volume_average("Nu_wb"))
        Nu_ε.append(stats.volume_average("Nu_ε"))
        Nu_χ.append(stats.volume_average("Nu_χ"))
        t.append(model.solver.sim_time)
    
        logger.info("i: {:d}, Nu_wb: {:.4f}, Nu_ε: {:.4f}, Nu_χ: {:.4f}".format(model.solver.iteration, Nu_wb[-1], Nu_ε[-1], Nu_χ[-1]))
                            
except:
    logger.error("Exception raised, triggering end of statistics-gathering loop.")
    raise

finally:
    save_arrays(nusselt_filename, Nu_wb=Nu_wb, Nu_ε=Nu_ε, Nu_χ=Nu_χ, t=t) 

stats_wall_time = time.time() - start_run_time
logger.info("""\n
    *** Statistics gathering complete.\n
    Statistics-gathering time: {:f} sec ({:.2f} hr)
    Total time: {:f} sec ({:.2f} hr)
    Run time: {:f} cpu-hr\n\n""".format(stats_wall_time, stats_wall_time/3600, spinup_wall_time+equil_wall_time+stats_wall_time,
        (spinup_wall_time+equil_wall_time+stats_wall_time)/3600, stats_time/3600 * model.domain.dist.comm_cart.size)
)
