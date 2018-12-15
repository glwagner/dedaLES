import numpy as np
import time, logging

from mpi4py import MPI
from dedalus import public as de

logger = logging.getLogger(__name__)

# A few globals
minute = 60
hour = 60*minute


class Flow():
    """
    A model for fluid flow.
    """
    def __init__(self, nx, ny, nz, Lx, Ly, Lz):
        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz


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

    Run a flow model (work in progress).
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

    

class ChannelFlow(Flow):
    """
    A models for flow in a channel, between two rigid walls.
    """
    def __init__(self, nx, ny, nz, Lx, Ly, Lz, xleft, yleft, zbottom): 

        Flow.__init__(self, nx, ny, nz, Lx, Ly, Lz)

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
