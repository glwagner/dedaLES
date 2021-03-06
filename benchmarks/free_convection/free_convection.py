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

debug = False

def identifier(model, closure=None): 
    if closure is None: closure_name = 'DNS'
    else:               closure_name = closure.__class__.__name__
    return "freeconvection_nh{:d}_nz{:d}_Q{:.0f}_Ninv{:.0f}_{:s}".format(
            model.nx, model.nz, -model.Q, 1/np.sqrt(initial_N2), closure_name)

# Input:
# python3 free_convection.py n L |Q| τ (closure)
#
# Args
# ----
#       n : resolution
#       L : domain size
#     |Q| : absolute value of surface heat flux
#       τ : buoyancy time-scale (1/N)
# closure : LES closure name (defaults to DNS if not supplied)

# Domain parameters
nx = ny = nz = int(sys.argv[1])     # Horizontal resolution
Lx = Ly = Lz = float(sys.argv[2])   # Domain extent [m]
Q = -float(sys.argv[3])             # Cooling rate [W m⁻²]
initial_N = 1/float(sys.argv[4])    # Buoyancy frequency [s⁻¹]

if len(sys.argv) > 5:
    closure_name = sys.argv[5]
    closure = getattr(dedaLES, closure_name)(stratified=True)
else:
    closure_name = 'DNS'
    closure = None

# Physical constants
a  = 1.0e-4     # Noise amplitude [m s⁻¹]
α  = 2.0e-4     # Thermal expansion coefficient for temperature [m s⁻² C⁻¹]
β  = 8.0e-4     # Thermal expansion coefficient for salinity [m s⁻² PSU⁻¹]
g  = 9.81       # Graviational acceleration [m s⁻²]
ρ0 = 1028.1     # Reference density [kg m⁻³]
cP = 3993.0     # Specific heat of oceanic water [J kg⁻¹ K⁻¹]
κ  = 1.43e-7    # Thermal diffusivity [m² s⁻¹]
ν  = 1.05e-6    # Viscosity [m² s⁻¹]

# Physical parameters
initial_T_surface = 20.0 # Initial surface temperature [C]
initial_Tz = initial_N**2 * ρ0 / (g*α) # Deep temperature gradient [C m⁻¹]
surface_flux = Q*α*g / (cP*ρ0*κ) # [s⁻²]
initial_N2 = initial_N**2

dt = 0.01 * np.sqrt(-1/surface_flux) # Initial time-step

# Numerical parameters
CFL_cadence = 100
stats_cadence = 10
analysis_cadence = 10
run_time = 2*hour

if debug:
    nx = ny = nz = 8
    CFL_cadence = np.inf
    dt = 1e-16
    run_time = 10*dt
    stats_cadence = analysis_cadence = 1

# Construct model
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, ν=ν, κ=κ, closure=closure,
                                      Q=Q, surface_flux=surface_flux, initial_N2=initial_N2)

# Boundary conditions
model.set_bc("no penetration", "top", "bottom")
model.set_bc("free slip", "top", "bottom")
model.set_tracer_gradient_bc("b", "top", gradient="surface_flux") #*tanh(t/tb)")
model.set_tracer_gradient_bc("b", "bottom", gradient="initial_N2")

model.build_solver(timestepper='SBDF3')

# Initial condition
noise = a * dedaLES.random_noise(model.domain) * model.z * (Lz - model.z) / Lz**2
initial_b = α * g * (initial_T_surface + initial_Tz*model.z)
model.set_fields(
    u = noise,
    v = noise,
    w = noise,
    b = initial_b + np.sqrt(initial_N2) * noise
)

model.stop_at(sim_time=run_time)

CFL = flow_tools.CFL(model.solver, initial_dt=dt, cadence=CFL_cadence, max_change=1.5, safety=0.5)
CFL.add_velocities(('u', 'v', 'w'))

stats = flow_tools.GlobalFlowProperty(model.solver, cadence=stats_cadence)

stats.add_property("w*b", name="buoyancyflux")
stats.add_property("ε + ε_sgs", name="epsilon")
stats.add_property("χ + χ_sgs", name="chi")
stats.add_property("w*w", name="wsquared")
stats.add_property("sqrt(u*u + v*v + w*w) / ν", name='Re')

analysis = model.solver.evaluator.add_file_handler(identifier(model, closure=closure), iter=analysis_cadence, max_writes=100)
analysis.add_system(model.solver.state, layout='g')
analysis.add_task("interp(b, y=0)", scales=1, name='b midplane')
analysis.add_task("interp(u, y=0)", scales=1, name='u midplane')
analysis.add_task("interp(v, y=0)", scales=1, name='v midplane')
analysis.add_task("interp(w, y=0)", scales=1, name='w midplane')

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
                "i: {:d}, t: {:.1f} hr, twall: {:.1f} s, dt: {:.1f} s, max Re {:.0f}, max sqrt(w^2): {:e}".format(
                        model.solver.iteration, model.solver.sim_time/hour, compute_time, dt, 
                        stats.max("Re"), np.sqrt(stats.max("wsquared"))
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
