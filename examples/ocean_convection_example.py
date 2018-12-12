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

# Domain parameters
nx = 128        # Horizontal resolution
ny = 16         # Horizontal resolution
nz = 64         # Vertical resolution
Lx = 200        # Domain horizontal extent [m]
Ly = 20         # Domain horizontal extent [m]
Lz = 100        # Domain vertical extent [m]

# Physical parameters
Q     = -100.0  # Cooling rate [W m⁻²]
N2inf = 9.5e-3  # Deep buoyancy gradient [s⁻²]
h0    = 50      # Initial mixed layer depth [m]
d     = 10      # Mixed layer - interior transition scale [m]

# Physical constants
α  = 2.5e-4     # Thermal expansion coefficient [K⁻¹]
g  = 9.81       # Graviational acceleration [m s⁻²]
ρ0 = 1028.1     # Reference density [kg m⁻³]
cP = 3993.0     # Specific heat of oceanic water [J kg⁻¹ K⁻¹]
κ  = 1.43e-7    # Thermal diffusivity [m² s⁻¹]
ν  = 1.05e-6    # Viscosity [m² s⁻¹]

# Calculated parameters
bz0 = Q*α*g / (cP*ρ0*κ) # [s⁻²]

# Construct model
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, ν=ν, κ=κ, 
                                      N2inf=N2inf, bz0=bz0, tb=day/2, closure=None)

def smoothstep(z, d): 
    return 0.5*(1 + np.tanh(z/d))

# Boundary conditions
model.set_bc("nopenetration", "top", "bottom")
model.set_bc("freeslip", "top", "bottom")
model.set_tracer_gradient_bc("b", "top", gradient="bz0") #*tanh(t/tb)")
model.set_tracer_gradient_bc("b", "bottom", gradient="N2inf")

model.build_solver()

# Initial condition
b0 = N2inf * (model.z + h0) * smoothstep(-model.z-h0, d)
model.set_b(b0)

def get_KE(model): return model.flow.volume_average('KE')
def get_variance(model): return model.flow.volume_average('variance')

#model.add_log_task("Average KE", get_KE)
#model.add_log_task("Buoyancy variance", get_variance)
#model.run(dt=0.1*cooling_scale, sim_time=4*day, log_cadence=10)

cooling_scale = 1/np.sqrt(-bz0)
print("Cooling time scale: {} s".format(cooling_scale))

dt = 10*cooling_scale
log_cadence = 100
model.solver.stop_wall_time = np.inf
model.solver.stop_iteration = np.inf
model.solver.stop_sim_time = 4*day

# Main loop
try:
    logger.info('Starting loop')
    start_run_time = time.time()

    while model.solver.ok:
        model.solver.step(dt)

        if model.time_to_log(log_cadence): 

            logger.info('Iter: {}, Time: {:.2f}, dt: {:.1f}'.format(
                        model.solver.iteration, model.solver.sim_time/hour, dt))

            logger.info("Average KE = {}".format(get_KE(model)))
            logger.info("Average buoyancy variance = {}".format(get_variance(model)))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise

finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %model.solver.iteration)
    logger.info('Sim end time: %f' %model.solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time) / hour * model.domain.dist.comm_cart.size))
        










































