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

if len(sys.argv) > 1:
    debug = True
else:
    debug = False

def identifier(model, closure=None): 
    if closure is None: closure_name = 'DNS'
    else:               closure_name = closure.__class__.__name__
    return "free_convection_nx{:d}_ny{:d}_nz{:d}_100Q{:.0f}_Ninv{:.0f}_{:s}".format(
            model.nx, model.ny, model.nz, -100*model.Q, 1/np.sqrt(initial_N2), closure_name)

# Main parameters
nx = nz = 128  # Horizontal resolution
Lx = Lz = 1.0  # Domain extent [m]
ny = 16        # Horizontal resolution
Ly = 0.1       # Domain extent [m]
Q = -0.1       # Cooling rate [W m⁻²]
initial_N2 = (1/1000.0)**2 # Initial buoyancy gradient [s⁻²]

# Physical constants
α  = 2.0e-4     # Thermal expansion coefficient [K⁻¹]
β  = 8.0e-4     # Thermal expansion coefficient [K⁻¹]
g  = 9.81       # Graviational acceleration [m s⁻²]
ρ0 = 1028.1     # Reference density [kg m⁻³]
cP = 3993.0     # Specific heat of oceanic water [J kg⁻¹ K⁻¹]
κ  = 1.43e-7    # Thermal diffusivity [m² s⁻¹]
ν  = 1.05e-6    # Viscosity [m² s⁻¹]

initial_T0   = 20.0                              # Surface temperature [C]
initial_Tz   = initial_N2 * ρ0 / (g*α)           # Temperature gradient [C m⁻¹]
surface_flux = Q*α*g / (cP*ρ0*κ)                 # Surface buoyancy flux [s⁻²]
initial_dt   = 0.01 / np.sqrt(-surface_flux)     # Initial time step
w_turb       = (-Lz*surface_flux)**(1/3)         # Domain turbulent velocity scale [m s⁻¹]
noise_amp    = 0.01*w_turb                       # Noise amplitude [m s⁻¹]
l_kolmo       = (-ν**3/surface_flux)**(1/4)      # Kolmogorov length scale
t_erosion    = -initial_N2*Lz**2/surface_flux    # Time-scale for stratification erosion

# Numerical parameters
CFL_cadence = 10
stats_cadence = 100
analysis_cadence = 100
run_time = 2*hour
max_writes = 1000

if debug:
    nx = ny = nz = 8
    CFL_cadence = np.inf
    initial_dt = 1e-16
    run_time = 10*initial_dt
    stats_cadence = analysis_cadence = 1


# Construct model
closure = None
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, ν=ν, κ=κ, closure=closure,
                                      Q=Q, surface_flux=surface_flux, initial_N2=initial_N2)

Δx_min, Δx_max = dedaLES.grid_stats(model, 0)
Δy_min, Δy_max = dedaLES.grid_stats(model, 1)
Δz_min, Δz_max = dedaLES.grid_stats(model, 2)

logger.info("""\n
    *** Convection into a linearly stratified fluid ***

                       Simulation info
                       ---------------

                            Q : {:.2e} m C s⁻¹
        surface buoyancy flux : {:.2e} m² s⁻³
                          1/N : {:.2e} s
                   initial dt : {:.2e} s
                     run time : {:.2e} s

                           Lx : {:.1f} m
                           Ly : {:.1f} m
                           Lz : {:.1f} m

                           nx : {:d}
                           ny : {:d}
                           nz : {:d}

            turbulent w-scale : {:.2e} m s⁻¹
           erosion time-scale : {:.2e} s
             Kolmogorov scale : {:.2e} m
           x-spacing min, max : {:.2e} m, {:.2e} m 
           y-spacing min, max : {:.2e} m, {:.2e} m
           z-spacing min, max : {:.2e} m, {:.2e} m

    """.format(Q, surface_flux, np.sqrt(1/initial_N2), initial_dt, run_time, 
               Lx, Ly, Lz, nx, ny, nz,
               w_turb, t_erosion, l_kolmo, Δx_min, Δx_max, Δy_min, Δy_max, Δz_min, Δz_max)
)


# Boundary conditions
model.set_bc("no penetration", "top", "bottom")
model.set_bc("free slip", "top", "bottom")
model.set_tracer_gradient_bc("b", "top", gradient="surface_flux") #*tanh(t/tb)")
model.set_tracer_gradient_bc("b", "bottom", gradient="initial_N2")

model.build_solver(timestepper='SBDF3')

# Initial condition
noise = noise_amp * dedaLES.random_noise(model.domain) * model.z / Lz #* (Lz - model.z) / Lz**2
initial_b = α*g*initial_Tz*model.z
model.set_fields(
    u = noise,
    v = noise,
    b = initial_b # + np.sqrt(initial_N2) * noise
    #w = noise,
)

model.stop_at(sim_time=run_time)

CFL = flow_tools.CFL(model.solver, initial_dt=initial_dt, cadence=CFL_cadence, max_change=1.5, safety=0.5)
CFL.add_velocities(('u', 'v', 'w'))

stats = flow_tools.GlobalFlowProperty(model.solver, cadence=stats_cadence)

stats.add_property("w*b", name="buoyancyflux")
stats.add_property("ε + ε_sgs", name="epsilon")
stats.add_property("χ + χ_sgs", name="chi")
stats.add_property("w*w", name="wsquared")
stats.add_property("sqrt(u*u + v*v + w*w) / ν", name='Re')

analysis = model.solver.evaluator.add_file_handler(identifier(model, closure=closure), iter=analysis_cadence, max_writes=max_writes)
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
        dt = CFL.compute_dt()
        model.solver.step(dt)

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
