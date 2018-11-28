"""
Dedalus script for 3D Rayleigh-Benard convection.

This script uses parity-bases in the x and y directions to mimick no-slip,
insulating sidewalls.  The equations are scaled in units of the thermal
diffusion time (Pe = 1).

This script should be ran in parallel, and would be most efficient using a
2D process mesh.  It uses the built-in analysis framework to save 2D data slices
in HDF5 files.  The `merge.py` script in this folder can be used to merge
distributed analysis sets from parallel runs, and the `plot_2d_series.py` script
can be used to plot the slices.

To run, merge, and plot using 4 processes, for instance, you could use:
    $ mpiexec -n 4 python3 rayleigh_benard.py
    $ mpiexec -n 4 python3 merge.py snapshots
    $ mpiexec -n 4 python3 plot_2d_series.py snapshots/*.h5

The simulation should take roughly 400 process-minutes to run, but will
automatically stop after an hour.

"""

import numpy as np
from mpi4py import MPI
import time

from dedalus import public as de
from dedalus.extras import flow_tools

import logging
logger = logging.getLogger(__name__)


# Parameters
nx, ny, nz = (64, 64, 16)
Lx, Ly, Lz = (400., 400., 100.)

f = 1e-4 # s⁻¹
κ = 1.43e-7 # b-m²/s
ν = 1e-6 # m²/s
Bz = 1e-6 # s⁻²

# Create bases and domain
start_init_time = time.time()
xbasis = de.Fourier('x', nx, interval=(-Lx/2, Lx/2), dealias=3/2)
ybasis = de.Fourier('y', ny, interval=(-Ly/2, Ly/2), dealias=3/2)
zbasis = de.Chebyshev('z', nz, interval=(-Lz, 0), dealias=3/2)
domain = de.Domain([xbasis, ybasis, zbasis], grid_dtype=np.float64)

# 3D rotating Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['p', 'b', 'u', 'v', 'w', 'bz', 'uz', 'vz', 'wz'], time='t')

problem.parameters['κ'] = κ    # diffusivity
problem.parameters['ν'] = ν    # kinematic viscosity
problem.parameters['f'] = f    # Coriolis parameter
problem.parameters['Bz'] = Bz  # background (unstable) buoyancy gradient

# Momentum equations
problem.add_equation("dt(u) - f*v - ν*(dx(dx(u)) + dy(dy(u)) + dz(uz)) + dx(p) = - u*dx(u) - v*dy(u) - w*uz")
problem.add_equation("dt(v) + f*u - ν*(dx(dx(v)) + dy(dy(v)) + dz(vz)) + dy(p) = - u*dx(v) - v*dy(v) - w*vz")
problem.add_equation("dt(w) - b   - ν*(dx(dx(w)) + dy(dy(w)) + dz(wz)) + dz(p) = - u*dx(w) - v*dy(w) - w*wz")

problem.add_equation("dt(b) - κ*(dx(dx(b)) + dy(dy(b)) + dz(bz)) = - u*dx(b) - v*dy(b) - w*bz")

problem.add_equation("dx(u) + dy(v) + wz = 0")

problem.add_equation("bz - dz(b) = 0")
problem.add_equation("uz - dz(u) = 0")
problem.add_equation("vz - dz(v) = 0")
problem.add_equation("wz - dz(w) = 0")

problem.add_bc("left(b) = -left(Bz*z)")
problem.add_bc("left(u) = 0")
problem.add_bc("left(v) = 0")
problem.add_bc("left(w) = 0")
problem.add_bc("right(b) = -right(Bz*z)")
problem.add_bc("right(u) = 0")
problem.add_bc("right(v) = 0")
problem.add_bc("right(w) = 0", condition="(nx != 0) or (ny != 0)")
problem.add_bc("integ_z(p) = 0", condition="(nx == 0) and (ny == 0)")

# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

# Initial conditions
x = domain.grid(0)
y = domain.grid(1)
z = domain.grid(2)
b = solver.state['b']
bz = solver.state['bz']

# Random perturbations, initialized globally for same results in parallel
gshape = domain.dist.grid_layout.global_shape(scales=1)
slices = domain.dist.grid_layout.slices(scales=1)
rand = np.random.RandomState(seed=23)
noise = rand.standard_normal(gshape)[slices]

# Linear background + perturbations damped at walls
zb, zt = zbasis.interval
pert =  1e-3 * noise * (zt - z) * (z - zb) / Lz
b['g'] = -Bz*(z - pert)
b.differentiate('z', out=bz)

# Integration parameters
solver.stop_sim_time = 100
solver.stop_wall_time = 60 * 60.
solver.stop_iteration = np.inf

# CFL
CFL = flow_tools.CFL(solver, initial_dt=1e-4, cadence=5, safety=1.5,
                     max_change=1.5, min_change=0.5, max_dt=0.05)
CFL.add_velocities(('u', 'v', 'w'))

flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("0.5*(u*u + v*v + w*w)", name='E')

# Main loop
end_init_time = time.time()
logger.info('Initialization time: %f' %(end_init_time-start_init_time))
try:
    logger.info('Starting loop')
    start_run_time = time.time()
    while solver.ok:
        dt = CFL.compute_dt()
        solver.step(dt)
        if (solver.iteration-1) % 100 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('Max KE = %f' %flow.max('E'))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/60/60*domain.dist.comm_cart.size))

