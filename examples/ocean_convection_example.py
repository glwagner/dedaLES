import sys; sys.path.append("..")

import time, logging
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi

logger = logging.getLogger(__name__)

import dedaLES

# Some convenient constants
second = 1.0
minute = 60*second
hour   = 60*minute
day    = 24*hour
siderealyear = 365*day + 6*hour + 9*min + 9.76*second
omega = 2*pi / sideralyear

# Physical constants
α  = 1.19e-4                # Thermal expansion coefficient for water at 10ᵒC [K⁻¹]
g  = 9.81                   # Graviational acceleration [m s⁻²]
ρ0 = 1027.62                # Reference density [kg m⁻³]
cP = 4003.0                 # Specific heat of seawater at 10ᵒC [J kg⁻¹ K⁻¹]
κ  = 1.43e-7                # Thermal diffusivity of seawater [m² s⁻¹]
ν  = 1.05e-6                # Viscosity of seawater [m² s⁻¹]

# Denbo and Skyllinstad (1995) experiment parameters
Q = -300.0                  # Surface cooling rate [W m⁻²]
Tz_deep = 3.2e-4            # Deep temperature gradient [ᵒC/m]
bz_deep = g*α/ρ0*Tz_deep    # Deep buoyancy gradient
h0 = 450.0                  # Initial mixed layer depth
d = h0/10                   # Transition thickness from deep stratified to mixed layers
f = 1.4e-4                  # Coriolis parameters [s⁻¹]

nx = ny = 128               # Horizontal resolution
nz = 34                     # Vertical resolution
Lx = Ly = 3840.0            # Domain horizontal extent [m]
Lz = 1020                   # Domain vertical extent [m]
a = 1e-2                    # Non-dimensional initial noise amplitude
dt = 90.0                   # Timestep

# Calculated parameters
bz_surf = Q*α*g / (cP*ρ0*κ) # Unstable surface buoyancy gradient [s⁻²]
b_init = (Lz-h0)*bz_deep    # Initial buoyancy at bottom of domain
b_noise = a*b_init          # Noise amplitude for initial buoyancy condition

# Construct model
closure = dedaLES.ConstantSmagorinsky()
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, ν=ν, κ=κ, zbottom=-Lz
                                      bz_deep=bz_deep, bz_surf=bz_surf, tb=day/2, closure=closure, H=Lz)
    
# Boundary conditions
model.set_bc("no penetration", "top", "bottom")
model.set_bc("free slip", "top", "bottom")
model.set_tracer_gradient_bc("b", "top", gradient="bz_surf") #*tanh(t/tb)")
model.set_tracer_gradient_bc("b", "bottom", gradient="bz_deep")

model.build_solver()

# Initial condition
def smoothstep(z, d): 
    return 0.5*(1 + np.tanh(z/d))

b0 = bz_deep * (model.z + h0) * smoothstep(-model.z-h0, d)

# Add noise ... ?
noise = dedaLES.random_noise(model.domain)
b0 += b_noise * noise * model.z * (model.z + Lz) / Lz**2 # max std dev b_noise/4
model.set_fields(b=b0)

# Flow properties
flow = flow_tools.GlobalFlowProperty(model.solver, cadence=100)
flow.add_property("sqrt(u*u + v*v + w*w)*H / nu", name='Re_domain')
flow.add_property("sqrt(ε) / nu", name='Re_dissipation')
flow.add_property("ε", name="dissipation")
flow.add_property("w*w", name="w_square")

model.stop_at(iteration=4800)

# Main loop
try:
    logger.info('Starting loop')
    start_run_time = time.time()

    while model.solver.ok:
        model.solver.step(dt)

        if model.time_to_log(100): 
            logger.info('Iter: {}, Time: {:.2f}, dt: {:.1f}'.format(
                        model.solver.iteration, model.solver.sim_time/hour, dt))
            logger.info("     Max domain Re = {:.6f}".format(flow.max("Re_domain")))
            logger.info("Max dissipation Re = {:.6f}".format(flow.max("Re_dissipation")))
            logger.info("   Average epsilon = {:.6f}".format(flow.volume_average("dissipation")))
            logger.info("     Average rms w = {:.6f}".format(np.sqrt(flow.volume_average("w_square"))))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise

finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %model.solver.iteration)
    logger.info('Sim end time: %f' %model.solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time) / hour * model.domain.dist.comm_cart.size))
        










































