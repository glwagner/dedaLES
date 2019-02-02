import os, sys; sys.path.append(os.getenv('DEDALES', os.path.join("..", "..")))

import time, logging
import numpy as np

from mpi4py import MPI
from numpy import pi
from dedalus.extras import flow_tools

import dedaLES

logger = logging.getLogger(__name__)

# Some convenient constants
second = 1.0
minute = 60*second
hour   = 60*minute
day    = 24*hour

try:
    if sys.argv[1] is 'debug':
        debug = True
except:
    debug = False

def identifier(model, closure=None): 
    if closure is None: closure_name = 'DNS'
    else:               closure_name = closure.__class__.__name__
    return "freeconvection_nh{:d}_nz{:d}_Q{:.0f}_bfreq{:.0f}_{:s}".format(
            model.nx, model.nz, -model.Q, 1/sqrt(initial_N2), closure_name)

# Domain parameters
nx = ny = 256   # Horizontal resolution
nz = 128        # Vertical resolution
Lx = Ly = 128.0 # Domain horizontal extent [m]
Lz = 150.0      # Domain horizontal extent [m]

# Physical parameters
Q = -75.0  # Cooling rate [W m⁻²]
initial_T_surface = 20.0 # Deep buoyancy gradient [s⁻²]
initial_Tz = 1e-2 # Deep buoyancy gradient [s⁻²]

# Physical constants
a  = 1e-2       # Noise amplitude
α  = 2.0e-4     # Thermal expansion coefficient [K⁻¹]
β  = 8.0e-4     # Thermal expansion coefficient [K⁻¹]
g  = 9.81       # Graviational acceleration [m s⁻²]
ρ0 = 1028.1     # Reference density [kg m⁻³]
cP = 3993.0     # Specific heat of oceanic water [J kg⁻¹ K⁻¹]
κ  = 1.43e-7    # Thermal diffusivity [m² s⁻¹]
ν  = 1.05e-6    # Viscosity [m² s⁻¹]
dt = 10*second

# Numerical parameters
CFL_cadence = 10
stats_cadence = 10
run_time = 4*hour

if debug:
    nx = ny = 16
    nz = 8
    CFL_cadence = np.inf

# Calculated parameters
surface_flux = Q*α*g / (cP*ρ0*κ) # [s⁻²]
initial_N2 = g*α*initial_Tz / ρ0

# Construct model
closure = None
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, ν=ν, κ=κ, closure=closure,
                                      surface_flux=surface_flux, initial_N2=initial_N2)

# Boundary conditions
model.set_bc("no penetration", "top", "bottom")
model.set_bc("free slip", "top", "bottom")
model.set_tracer_gradient_bc("b", "top", gradient="surface_flux") #*tanh(t/tb)")
model.set_tracer_gradient_bc("b", "bottom", gradient="initial_N2")

model.build_solver(timestepper='SBDF3')

# Initial condition
noise = a * dedaLES.random_noise(model.domain) * model.z * (Lz - model.z) / Lz**2,
initial_b = α * g * (initial_T + initial_Tz*model.z)
model.set_fields(
    u = noise,
    v = noise,
    w = noise,
    b = initial_b * (1 + noise)
)

model.stop_at(sim_time=run_time)

CFL = flow_tools.CFL(model.solver, initial_dt=dt, cadence=CFL_cadence, max_change=1.5, safety=0.5)
CFL.add_velocities(('u', 'v', 'w'))

stats = flow_tools.GlobalFlowProperty(model.solver, cadence=stats_cadence)

stats.add_property("w*b", name="buoyancyflux")
stats.add_property("ε + ε_sgs", name="epsilon")
stats.add_property("χ + χ_sgs", name="chi")
stats.add_property("w**2", name="wsquared")
stats.add_property("sqrt(u*u + v*v + w*w) / ν", name='Re')


analysis = model.solver.evaluator.add_file_handler(identifier(model, closure=closure), iter=1000, max_writes=100)
analysis.add_system(model.solver.state, layout='g')

# Main loop
try:
    logger.info('Starting loop')
    start_run_time = time.time()
    log_time = time.time()  

    while model.solver.ok:
        model.solver.step(dt)
        dt = CFL.compute_dt()

        if model.time_to_log(stats_cadence): 
            compute_time = time.time() - log_time
            log_time = time.time()

            logger.info(
                "i: {:d}, t: {:.1f} hr, twall: {:.1s} s, dt: {:.1f} s, max Re {:.0f}, max w^2: {:.2f}".format(
                        model.solver.iteration, model.solver.sim_time/hour, log_time, dt, 
                        stats.max("Re"), stats.max("wsquared")
            ))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise

finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %model.solver.iteration)
    logger.info('Sim end time: %f' %model.solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time) / hour * model.domain.dist.comm_cart.size))
