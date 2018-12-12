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

    # Boundary condition API.
    def set_nopenetration_bc_top(self):
        self.problem.add_bc("right(w) = 0", condition="(nx != 0) or (ny != 0)")

    def set_nopenetration_bc_bottom(self):
        self.problem.add_bc("left(w) = 0")

    def set_noslip_bc_top(self):
        self.problem.add_bc("right(u) = 0")
        self.problem.add_bc("right(v) = 0")

    def set_noslip_bc_bottom(self):
        self.problem.add_bc("left(u) = 0")
        self.problem.add_bc("left(v) = 0")
    
    def set_freeslip_bc_top(self):
        self.problem.add_bc("right(uz) = 0")
        self.problem.add_bc("right(vz) = 0")

    def set_freeslip_bc_bottom(self):
        self.problem.add_bc("left(uz) = 0")
        self.problem.add_bc("left(vz) = 0")

    def set_velocity_bc_top(self, u=0, v=0):
        self.problem.add_bc("right(u) = {}".format(u))
        self.problem.add_bc("right(v) = {}".format(v))

    def set_velocity_bc_bottom(self, u=0, v=0):
        self.problem.add_bc("left(u) = {}".format(u))
        self.problem.add_bc("left(v) = {}".format(v))

    def set_tracer_gradient_bc_bottom(self, tracers, gradient=0):
        if isinstance(tracers, str):
            self.problem.add_bc("left({}z) = {}".format(tracers, gradient))
        else:
            for tracer in tracers:
                self.problem.add_bc("left({}z) = {}".format(tracer, gradient))

    def set_tracer_gradient_bc_top(self, tracers, gradient=0):
        if isinstance(tracers, str):
            self.problem.add_bc("right({}z) = {}".format(tracers, gradient))
        else:
            for tracer in tracers:
                self.problem.add_bc("right({}z) = {}".format(tracer, gradient))

    def set_tracer_gradient_bc(self, tracers, *walls, gradient=0):
        for wall in walls:
            method = getattr(self, "set_tracer_gradient_bc_%s" %wall)
            method(tracers, gradient=gradient)

    def set_tracer_noflux_bc(self, tracers, *walls):
        self.set_tracer_gradient_bc(tracers, *walls, gradient=0)

    def set_bc(self, bctype, *walls, **kwargs):
        """
        Set boundary conditions for the velocity field in a channel model.

        Args
        ----
            bctype : A string that indicates the type of boundary condition being specified.

            walls : The walls on which the boundary condition will be specified.
                    (either "top" or "bottom").

        Keyword args
        ------------
            u, v, etc: Keyword arguments associated with the specified boundary condition.

        Notes
        -----
            This method sets boundary condition on the momentum variables :math:`u`, :math:`v`,
            and :math:`w`.

            Boundary condition types that are currently implemented are
                * "nopenetration" or :math:`w=0`
                * "noslip" or :math:`u=v=0`
                * "freeslip" or :math:`uz=vz=0`
                * "velocity" or :math:`u=u_0` and :math:`v=v_0`

            The velocity boundary condition accepts keyword arguments ``u`` and ``v``.
        """
        for wall in walls:
            method = getattr(self, 'set_%s_bc_%s' %(bctype, wall))
            method(**kwargs)

    def set_tracer_bc(self, bctype, tracers, *walls, **kwargs):
        """
        Set boundary conditions for a tracer in a channel model.

        Args
        ----
            bctype : A string that indicates the type of boundary condition being specified.

            tracers : A string or list of strings that indicate the names of the tracer fields
                      that the boundary condition applies to.

            walls : The walls on which the boundary condition will be specified
                    (either "top" or "bottom").

            kwargs : Keyword arguments associated with the boundary condition.

        Notes
        -----
            This method sets boundary condition on the tracer variables indicated by
            ``tracer``.

            Currently implemented boundary condition types for a tracer :math:`c` are
                * "gradient" or :math:`cz=G`
        """
        for wall in walls:
            method = getattr(self, 'set_tracer_%s_bc_%s' %(bctype, wall))
            method(tracers, **kwargs)

    def set_field(self, phi, gridvalue):
        """
        Set the state field `phi` as `gridvalue`. Calculate the derivative `phiz`.

        Args
        ----
            phi : The name of the field to be set.

            gridvalue : The grid values of phi.
        """

        field = getattr(self, phi)
        field['g'] = gridvalue

        fieldz = getattr(self, phi+'z')
        field.differentiate('z', out=fieldz)

    def set_fields(self, **fields):
        for name, value in fields.items():
            self.set_field(name, value)

    def build_solver(self, timestepper='RK443'):
        """
        Build a dedalus solver for the model with `timestepper`.

        Args
        ----
            timestepper : The name of the timestepper to be used by the solver. 
        """
        detimestepper = getattr(de.timesteppers, timestepper)

        start_build_time = time.time()
        solver = self.problem.build_solver(detimestepper)
        logger.info('Solver built. (t = %f) ' %(time.time()-start_build_time))

        self.solver = solver

    def log(self, logger, dt):
        pass

    def time_to_log(self, log_cadence):
        return (self.solver.iteration-1) % log_cadence == 0

    def init_flow(self, cadence=10):
        self.flow = flow_tools.GlobalFlowProperty(self.solver, cadence=cadence)

    def add_flow_KE(self):
        self.flow.add_property("0.5*(u*u + v*v + w*w)", name="KE")

    """
        run(self, **kwargs)

    Run a ChannelFlow model (work in progress).
    """
    def run(self, 
                     dt = 1e-16, 
             iterations = 100, 
             sim_time   = None, 
             wall_time  = None, 
            log_cadence = 100
        ):

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

        try:
            logger.info('Starting loop')
            start_run_time = time.time()
            while self.solver.ok:
                self.solver.step(dt)
                if self.time_to_log(log_cadence): self.log(logger, dt)
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
    """
    Boussinesq flow in a channel (with optional rotation).

    Parameters
    ----------
        nx       : grid resolution in :math:`x`
        ny       : grid resolution in :math:`y`
        nz       : grid resolution in :math:`z`
        Lx       : domain extent in :math:`x`
        Ly       : domain extent in :math:`y`
        Lz       : domain extent in :math:`z`
        f        : Coriolis parameter
        ν        : 'Molecular' viscosity
        κ        : 'Molecular' diffusivity for buoyancy
        Nsq      : Background buoyancy gradient
        closure  : turbulent closure for Large Eddy Simulation
        **params : additional parameters to be added to the dedalus problem.
    """
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
        self.problem = problem = de.IVP(domain, variables=['p', 'b', 'u', 'v', 'w', 'bz', 'uz', 'vz', 'wz'], time='t')

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
        if closure is None:
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

        self.log_tasks = {}
        
    def add_log_task(self, name, task):
        self.log_tasks[name] = task

    def add_flow_variance(self):
        self.flow.add_property("0.5*b*b", name="variance")

    def set_noflux_bc_top(self):
        ChannelFlow.set_tracer_noflux_bc(self, "b", "top")

    def set_noflux_bc_bottom(self):
        ChannelFlow.set_tracer_noflux_bc(self, "b", "bottom")

    def set_default_bcs(self):
        self.set_bc("nopenetration", "top", "bottom")
        self.set_bc("noslip", "top", "bottom")
        self.set_bc("noflux", "top", "bottom")

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

        self.init_flow()
        self.add_flow_KE()
        self.add_flow_variance()


    def log(self, logger, dt):
        logger.info('Iteration: %i, Time: %e, dt: %e' %(self.solver.iteration, self.solver.sim_time, dt))

        for name, task in self.log_tasks.items(): 
            logger.info("{} = {}".format(name, task(self)))
