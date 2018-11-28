"""
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

# A few globals
minute = 60
hour = 60*minute

def tanhstep(z, z0, d):
    return 0.5*(np.tanh((z-z0)/d) + 1)

def mixedlayerN(z, dmix, Nmix, Ndeep):
    return Nmix + (Ndeep-Nmix)*tanhstep(z, -dmix, -0.2*dmix)

def addparams(problem, **params):
    for k, v in params.items():
        problem.parameters[k] = v
        
class OceanModel():
    def __init__(self,
        nx = 128, 
        ny = 128,
        nz = 16, 
        Lx = 800.0,     # Domain length [m] 
        Ly = 800.0,     # Domain length [m] 
        Lz = 100.0,     # Domain height [m]
        f = 1e-4,       # Coriolis parameters [s⁻¹]
        κ = 1.43e-7,    # Diffusivity [m²/s]
        ν = 1e-6,       # Kinematic viscosity [m²/s]
        **params):

        # Create bases and domain
        xbasis = de.Fourier('x', nx, interval=(-Lx/2, Lx/2), dealias=3/2)
        ybasis = de.Fourier('y', ny, interval=(-Ly/2, Ly/2), dealias=3/2)
        zbasis = de.Chebyshev('z', nz, interval=(-Lz, 0), dealias=3/2)
        domain = de.Domain([xbasis, ybasis, zbasis], grid_dtype=np.float64)

        self.xbasis = xbasis
        self.ybasis = ybasis
        self.zbasis = zbasis
        self.domain = domain

        self.x = domain.grid(0)
        self.y = domain.grid(1)
        self.z = domain.grid(2)

        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz

        self.f = f
        self.κ = κ
        self.ν = ν

        # 3D rotating Boussinesq hydrodynamics
        problem = de.IVP(domain, variables=['p', 'b', 'u', 'v', 'w', 'bz', 'uz', 'vz', 'wz'], time='t')

        addparams(problem, κ=κ, ν=ν, f=f, **params)
        for k, v in params.items():
            setattr(self, k, v)

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

        self.problem = problem

    def build_solver(self, timestepper=de.timesteppers.RK443):

        start_build_time = time.time()
        solver = self.problem.build_solver(timestepper)
        logger.info('Solver built. (t = %f) ' %(time.time()-start_init_time)))

        self.solver = solver

        self.u = solver.state['u']
        self.v = solver.state['v']
        self.w = solver.state['w']
        self.b = solver.state['b']

        self.uz = solver.state['uz']
        self.vz = solver.state['vz']
        self.wz = solver.state['wz']
        self.bz = solver.state['bz']


    def run(self, initial_dt=1e-4, sim_time=100): 

        # Integration parameters
        self.solver.stop_sim_time = self.solver.sim_time + sim_time
        self.solver.stop_wall_time = np.inf # 60 * 60.
        self.solver.stop_iteration = np.inf

        # CFL
        CFL = flow_tools.CFL(self.solver, initial_dt=initial_dt, cadence=5, safety=1.5,
                             max_change=1.5, min_change=0.5, max_dt=0.05)
        CFL.add_velocities(('u', 'v', 'w'))

        flow = flow_tools.GlobalFlowProperty(self.solver, cadence=10)
        flow.add_property("0.5*(u*u + v*v + w*w)", name='KE')

        self.CFL = CFL
        self.flow = flow

        try:
            logger.info('Starting loop')
            start_run_time = time.time()
            while self.solver.ok:
                dt = CFL.compute_dt()
                self.solver.step(dt)
                if (self.solver.iteration-1) % 100 == 0:
                    logger.info('Iteration: %i, Time: %e, dt: %e' %(self.solver.iteration, self.solver.sim_time, dt))
                    logger.info('Average KE = %e' %flow.volume_average('KE'))
        except:
            logger.error('Exception raised, triggering end of main loop.')
            raise
        finally:
            end_run_time = time.time()
            logger.info('Iterations: %i' %self.solver.iteration)
            logger.info('Sim end time: %f' %self.solver.sim_time)
            logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
            logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time) / hour * self.domain.dist.comm_cart.size))


class DeepConvectionModel(OceanModel):
    def __init__(self,
        nx = 128, 
        ny = 128,
        nz = 16, 
        Lx = 800.0,     # Domain length [m] 
        Ly = 800.0,     # Domain length [m] 
        Lz = 100.0,     # Domain height [m]
        f = 1e-4,       # Coriolis parameters [s⁻¹]
        κ = 1.43e-7,    # Diffusivity [m²/s]
        ν = 1e-6,       # Kinematic viscosity [m²/s]
        bflux = 0.0,      # Buoyancy gradient [s⁻²]
        ):

        OceanModel.__init__(self, nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz, f=f, κ=κ, ν=ν, bflux=bflux)

        self.problem.add_bc("left(bz) = bflux")
        self.problem.add_bc("left(u) = 0")
        self.problem.add_bc("left(v) = 0")
        self.problem.add_bc("left(w) = 0")
        self.problem.add_bc("right(b) = 0")
        self.problem.add_bc("right(u) = 0")
        self.problem.add_bc("right(v) = 0")
        self.problem.add_bc("right(w) = 0", condition="(nx != 0) or (ny != 0)")
        self.problem.add_bc("integ_z(p) = 0", condition="(nx == 0) and (ny == 0)")

        self.build_solver()


class RayleighBernardConvection(OceanModel):
    def __init__(self,
        nx = 128, 
        ny = 128,
        nz = 16, 
        Lx = 800.0,     # Domain length [m] 
        Ly = 800.0,     # Domain length [m] 
        Lz = 100.0,     # Domain height [m]
        f = 1e-4,       # Coriolis parameters [s⁻¹]
        κ = 1.43e-7,    # Diffusivity [m²/s]
        ν = 1e-6,       # Kinematic viscosity [m²/s]
        Bz = 1e-6,      # Buoyancy gradient [s⁻²]
        ):

        OceanModel.__init__(self, nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz, f=f, κ=κ, ν=ν, Bz=Bz)

        self.problem.add_bc("left(b) = -left(Bz*z)")
        self.problem.add_bc("left(u) = 0")
        self.problem.add_bc("left(v) = 0")
        self.problem.add_bc("left(w) = 0")
        self.problem.add_bc("right(b) = -right(Bz*z)")
        self.problem.add_bc("right(u) = 0")
        self.problem.add_bc("right(v) = 0")
        self.problem.add_bc("right(w) = 0", condition="(nx != 0) or (ny != 0)")
        self.problem.add_bc("integ_z(p) = 0", condition="(nx == 0) and (ny == 0)")

        self.build_solver()
        

    def set_unstable_ic(self, magnitude=1e-3):

        # Random perturbations, initialized globally for same results in parallel
        gshape = self.domain.dist.grid_layout.global_shape(scales=1)
        slices = self.domain.dist.grid_layout.slices(scales=1)
        rand = np.random.RandomState(seed=23)
        noise = rand.standard_normal(gshape)[slices]

        # Linear background + perturbations damped at walls
        zb, zt = self.zbasis.interval
        pert =  magnitude * noise * (zt - self.z) * (self.z - zb) / self.Lz
        self.b['g'] = -self.Bz*(self.z - pert)
        self.b.differentiate('z', out=self.bz)

