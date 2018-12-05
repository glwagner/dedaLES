import numpy as np
import time, logging
from mpi4py import MPI

from dedalus import public as de
from dedalus.extras import flow_tools

logger = logging.getLogger(__name__)


# A few globals
minute = 60
hour = 60*minute


def addparams(problem, **params):
    for k, v in params.items():
        problem.parameters[k] = v


class ChannelFlow():
    def __init__(self):
        pass

    # Boundary conditions:
    def set_nopenetration_top(self):
        self.problem.add_bc("right(w) = 0", condition="(nx != 0) or (ny != 0)")

    def set_nopenetration_bottom(self):
        self.problem.add_bc("left(w) = 0")

    def set_nopenetration_topandbottom(self):
        self.set_nopenetration_top()
        self.set_nopenetration_bottom()

    def set_noslip_top(self):
        self.problem.add_bc("right(u) = 0")
        self.problem.add_bc("right(v) = 0")

    def set_noslip_bottom(self):
        self.problem.add_bc("left(u) = 0")
        self.problem.add_bc("left(v) = 0")

    def set_noslip_topandbottom(self):
        self.set_noslip_top()
        self.set_noslip_bottom()

    def set_freeslip_top(self):
        self.problem.add_bc("right(uz) = 0")
        self.problem.add_bc("right(vz) = 0")

    def set_freeslip_bottom(self):
        self.problem.add_bc("left(uz) = 0")
        self.problem.add_bc("left(vz) = 0")

    def set_freeslip_topandbottom(self):
        self.set_freeslip_top()
        self.set_freeslip_bottom()
    
    def set_wall_velocity_top(self, utop=0, vtop=0):
        self.parameters['utop'] = utop
        self.parameters['vtop'] = vtop

        self.problem.add_bc("right(u) = utop")
        self.problem.add_bc("right(v) = vtop")

    def set_wall_velocity_bottom(self, ubottom=0, vbottom=0):
        self.parameters['ubottom'] = ubottom
        self.parameters['vbottom'] = vbottom

        self.problem.add_bc("left(u) = ubottom")
        self.problem.add_bc("left(v) = vbottom")

    def build_solver(self, timestepper='RK443'):

        detimestepper = getattr(de.timesteppers, timestepper)
        start_build_time = time.time()
        solver = self.problem.build_solver(detimestepper)
        logger.info('Solver built. (t = %f) ' %(time.time()-start_build_time))

        self.solver = solver

    def log(self, logger, dt):
        logger.info('Iteration: %i, Time: %e, dt: %e' %(self.solver.iteration, self.solver.sim_time, dt))
        logger.info('Average KE = %e' %self.flow.volume_average('KE'))

    def time_to_log(self, logcadence):
        return (self.solver.iteration-1) % logcadence == 0

    def add_flow_properties(self, flow):
        flow.add_property("0.5*(u*u + v*v + w*w)", name='KE')

    def run(self, initial_dt=1e-16, sim_time=None, iterations=100, wall_time=None, logcadence=100):

        # Integration parameters
        if wall_time is not None:
            self.solver.stop_wall_time = wall_time
        else:
            self.solver.stop_wall_time = np.inf

        if sim_time is not None:
            self.solver.stop_sim_time = self.solver.sim_time + sim_time
            self.solver.stop_iteration = np.inf
        else:
            self.solver.stop_sim_time = np.inf
            self.solver.stop_iteration = self.solver.iteration + iterations

        # CFL
        CFL = flow_tools.CFL(self.solver, initial_dt=initial_dt, cadence=5, safety=1.5, 
                             max_change=1.5, min_change=0.5, max_dt=0.05)
        CFL.add_velocities(('u', 'v', 'w'))

        flow = flow_tools.GlobalFlowProperty(self.solver, cadence=10)
        self.add_flow_properties(flow)

        self.CFL = CFL
        self.flow = flow

        try:
            logger.info('Starting loop')
            start_run_time = time.time()
            while self.solver.ok:
                dt = CFL.compute_dt()
                self.solver.step(dt)
                if self.time_to_log(logcadence): self.log(logger, dt)
                    
        except:
            logger.error('Exception raised, triggering end of main loop.')
            raise
        finally:
            end_run_time = time.time()
            logger.info('Iterations: %i' %self.solver.iteration)
            logger.info('Sim end time: %f' %self.solver.sim_time)
            logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
            logger.info(
                'Run time: %f cpu-hr' %((end_run_time-start_run_time) / hour * self.domain.dist.comm_cart.size))

    
                            


class BoussinesqChannelFlow(ChannelFlow):
    def __init__(self,
        nx = 64,
        ny = 64,
        nz = 64,
        Lx = 1.0,        # [m]
        Ly = 1.0,        # [m]
        Lz = 1.0,        # [m]
        f = 0.0,         # Coriolis parameter [s⁻¹]
        κ = 1.43e-7,     # Diffusivity [m²/s]
        ν = 1e-6,        # Kinematic viscosity [m²/s]
        Nsq = 0.0,       # Background buoyancy gradient [s⁻²]
        closure = None,  # Subgrid closure
        **params         # Additional parameters to be added to dedalus problem
        ):

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
        self.problem = problem = de.IVP(domain, 
                                        variables=['p', 'b', 'u', 'v', 'w', 'bz', 'uz', 'vz', 'wz'], time='t')

        # Add additional parameters passed to constructor to the dedalus Problem.
        addparams(problem, κ=κ, ν=ν, f=f, Nsq=Nsq, **params)
        for name, value in params.items():
            setattr(self, name, value)

        # First-order substitutions
        problem.substitutions['ux'] = "dx(u)"
        problem.substitutions['vx'] = "dx(v)"
        problem.substitutions['wx'] = "dx(w)"
        problem.substitutions['bx'] = "dx(b)"
        problem.substitutions['uy'] = "dy(u)"
        problem.substitutions['vy'] = "dy(v)"
        problem.substitutions['wy'] = "dy(w)"
        problem.substitutions['by'] = "dy(b)"

        # Closure substitutions
        if closure in (None, "DNS"):
            # No parameterized subgrid fluxes
            problem.substitutions['Fx_sgs'] = "0"
            problem.substitutions['Fy_sgs'] = "0"
            problem.substitutions['Fz_sgs'] = "0"
            problem.substitutions['Fb_sgs'] = "0"
        else:
            closure.substitutions(problem, tracers=["b"])

        # Momentum equations
        problem.add_equation("dt(u) - f*v   - ν*(dx(ux) + dy(uy) + dz(uz)) + dx(p) = - u*ux - v*uy - w*uz + Fx_sgs")
        problem.add_equation("dt(v) + f*u   - ν*(dx(vx) + dy(vy) + dz(vz)) + dy(p) = - u*vx - v*vy - w*vz + Fy_sgs")
        problem.add_equation("dt(w) - b     - ν*(dx(wx) + dy(wy) + dz(wz)) + dz(p) = - u*wx - v*wy - w*wz + Fz_sgs")
        problem.add_equation("dt(b) + Nsq*w - κ*(dx(bx) + dy(by) + dz(bz))         = - u*bx - v*by - w*bz + Fb_sgs")
        problem.add_equation("ux + vy + wz = 0")

        problem.add_equation("bz - dz(b) = 0")
        problem.add_equation("uz - dz(u) = 0")
        problem.add_equation("vz - dz(v) = 0")
        problem.add_equation("wz - dz(w) = 0")

        problem.add_bc("right(p) = 0", condition="(nx == 0) and (ny == 0)")

    def set_noflux_bottom(self):
        self.problem.add_bc("left(bz) = 0")

    def set_noflux_top(self):
        self.problem.add_bc("right(bz) = 0")

    def set_noflux_topandbottom(self):
        self.set_noflux_top()
        self.set_noflux_bottom()

    def set_default_bcs(self):
        self.set_nopenetration_topandbottom()
        self.set_noslip_topandbottom()
        self.set_noflux_topandbottom()

    def build_solver(self, timestepper='RK443'):

        ChannelFlow.build_solver(self, timestepper=timestepper)

        self.u = self.solver.state['u']
        self.v = self.solver.state['v']
        self.w = self.solver.state['w']
        self.b = self.solver.state['b']

        self.uz = self.solver.state['uz']
        self.vz = self.solver.state['vz']
        self.wz = self.solver.state['wz']
        self.bz = self.solver.state['bz']

    def set_b(self, b0):
        self.b['g'] = b0
        self.b.differentiate('z', out=self.bz)
