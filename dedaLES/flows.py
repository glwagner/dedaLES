import numpy as np
import time, logging

from mpi4py import MPI
from dedalus import public as de

logger = logging.getLogger(__name__)


# A few globals
minute = 60
hour = 60*minute

def add_parameters(problem, **params):
    """
    Add parameters to a dedalus problem programmatically.
    """
    for name, value in params.items():
        problem.parameters[name] = value

def add_derivative_substitutions(problem, dim, vars=[]):
    for var in vars:
        problem.substitutions[var+dim] = f"d{dim}({var})"

def add_closure(problem, closure, tracers=[]):
    """
    Add closure substitutions and equations to the problem.
    """
    if closure is None:
        # No parameterized subgrid fluxes
        problem.substitutions['Fx_sgs'] = "0"
        problem.substitutions['Fy_sgs'] = "0"
        problem.substitutions['Fz_sgs'] = "0"
        problem.substitutions['Fb_sgs'] = "0"
    else:
        closure.substitutions(problem, tracers=tracers)

        
def bind_parameters(obj, **params):
    """
    Bind the name, value pairs in `params` as members of the class `obj`.
    """
    for name, value in params.items():
        setattr(obj, name, value)


class ChannelFlow():
    """
    A parent class for fluid models in channels.
    """
    def __init__(self, nx, ny, nz, Lx, Ly, Lz, xleft, yleft, zbottom): 

        # Default origin: (-Lx/2, -Ly/2, -Lz)
        if xleft is None: xleft = -Lx/2
        if yleft is None: yleft = -Ly/2
        if zbottom is None: zbottom = 0.0

        xright = xleft + Lx
        yright = yleft + Ly
        ztop = zbottom + Lz

        self.xlimits = (xleft, xright)
        self.ylimits = (yleft, yright)
        self.zlimits = (zbottom, ztop)

        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz

    # Boundary condition API.
    def set_nopenetration_bc_top(self):
        # Forbid domain-averaged vertical velocity.
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

    def set_tracer_gradient_bc_bottom(self, tracers=[], gradient=0):
        if isinstance(tracers, str):
            self.problem.add_bc("left({}z) = {}".format(tracers, gradient))
        else:
            for tracer in tracers:
                self.problem.add_bc("left({}z) = {}".format(tracer, gradient))

    def set_tracer_gradient_bc_top(self, tracers=[], gradient=0):
        if isinstance(tracers, str):
            self.problem.add_bc("right({}z) = {}".format(tracers, gradient))
        else:
            for tracer in tracers:
                self.problem.add_bc("right({}z) = {}".format(tracer, gradient))

    def set_tracer_gradient_bc(self, tracers=[], *walls, gradient=0):
        for wall in walls:
            method = getattr(self, "set_tracer_gradient_bc_%s" %wall)
            method(tracers, gradient=gradient)

    def set_tracer_noflux_bc(self, *walls, tracers=[]):
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

    def set_tracer_bc(self, bctype, *walls, tracers=[], **kwargs):
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

    def add_log_task(self, name, task):
        try:
            self.log_tasks[name] = task
        except AttributeError:
            self.log_tasks = {name: task}

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
        nx  = 32,
        ny  = 32,
        nz  = 32,
        Lx  = 1.0,       # [m]
        Ly  = 1.0,       # [m]
        Lz  = 1.0,       # [m]
        f   = 0.0,       # Coriolis parameter [s⁻¹]
        κ   = 1.43e-7,   # Diffusivity [m²/s]
        ν   = 1.05e-6,   # Kinematic viscosity [m²/s]
        Nsq = 0.0,       # Background buoyancy gradient [s⁻²]
        closure = None,  # Subgrid closure
        xleft   = None,  # x-origin of domain (default: -Lx/2)
        yleft   = None,  # y-origin of domain (default: -Ly/2)
        zbottom = None,  # z-coordinate of channel bottom (default is -Lz)
        **params         # Additional parameters to be added to dedalus problem
        ):

        ChannelFlow.__init__(self, nx, ny, nz, Lx, Ly, Lz, xleft, yleft, zbottom)

        # Create bases and domain
        self.xbasis = xbasis = de.Fourier('x', nx, interval=self.xlimits, dealias=3/2)
        self.ybasis = ybasis = de.Fourier('y', ny, interval=self.ylimits, dealias=3/2)
        self.zbasis = zbasis = de.Chebyshev('z', nz, interval=self.zlimits, dealias=3/2)
        self.domain = domain = de.Domain([xbasis, ybasis, zbasis], grid_dtype=np.float64)

        self.x = domain.grid(0)
        self.y = domain.grid(1)
        self.z = domain.grid(2)

        # 3D rotating Boussinesq hydrodynamics
        self.problem = problem = de.IVP(domain, variables=['p', 'b', 'u', 'v', 'w', 'bz', 'uz', 'vz', 'wz'], time='t')

        # Add additional parameters passed to constructor to the dedalus Problem.
        add_parameters(problem, κ=κ, ν=ν, f=f, Nsq=Nsq, **params)
        bind_parameters(self, f=f, κ=κ, ν=ν, **params)

        # First-order substitutions
        add_derivative_substitutions(problem, 'x', vars=['u', 'v', 'w', 'b'])
        add_derivative_substitutions(problem, 'y', vars=['u', 'v', 'w', 'b'])

        # LES Closure
        add_closure(problem, closure, tracers=['b'])
                   
        # Momentum equations
        problem.add_equation("dt(u) - ν*(dx(ux) + dy(uy) + dz(uz)) + dx(p) - f*v = - u*ux - v*uy - w*uz + Fx_sgs")
        problem.add_equation("dt(v) - ν*(dx(vx) + dy(vy) + dz(vz)) + dy(p) + f*u = - u*vx - v*vy - w*vz + Fy_sgs")
        problem.add_equation("dt(w) - ν*(dx(wx) + dy(wy) + dz(wz)) + dz(p) - b   = - u*wx - v*wy - w*wz + Fz_sgs")

        # Buoyancy equation
        problem.add_equation("dt(b) - κ*(dx(bx) + dy(by) + dz(bz)) + Nsq*w       = - u*bx - v*by - w*bz + Fb_sgs")

        # Continuity equation
        problem.add_equation("ux + vy + wz = 0")

        # First-order equivalencies
        problem.add_equation("bz - dz(b) = 0")
        problem.add_equation("uz - dz(u) = 0")
        problem.add_equation("vz - dz(v) = 0")
        problem.add_equation("wz - dz(w) = 0")

        # Pressure gauge condition
        problem.add_bc("right(p) = 0", condition="(nx == 0) and (ny == 0)")
        
    def set_noflux_bc_top(self):
        """
        A convenience method to set a no flux condition on buoyancy
        on the top wall of the channel.
        """
        ChannelFlow.set_tracer_noflux_bc(self, "top", tracers="b")

    def set_noflux_bc_bottom(self):
        """
        A convenience method to set a no flux condition on buoyancy
        on the bottom wall of the channel.
        """
        ChannelFlow.set_tracer_noflux_bc(self, "bottom", tracers="b")

    def build_solver(self, timestepper='RK443'):
        """
        Build a BoussinesqChanelFlow solver.
        """

        ChannelFlow.build_solver(self, timestepper=timestepper)

        self.u = self.solver.state['u']
        self.v = self.solver.state['v']
        self.w = self.solver.state['w']
        self.b = self.solver.state['b']

        self.uz = self.solver.state['uz']
        self.vz = self.solver.state['vz']
        self.wz = self.solver.state['wz']
        self.bz = self.solver.state['bz']


    def log(self, logger, dt):
        """
        Print messages.
        """
        logger.info('Iteration: %i, Time: %e, dt: %e' %(self.solver.iteration, self.solver.sim_time, dt))

        for name, task in self.log_tasks.items(): 
            logger.info("{} = {}".format(name, task(self)))
