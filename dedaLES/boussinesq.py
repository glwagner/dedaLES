import numpy as np

from dedalus import public as de

from .flows import ChannelFlow
from .utils import add_parameters, bind_parameters, add_substitutions, add_first_derivative_substitutions
from .closures import add_closure_substitutions, add_closure_variables, add_closure_equations

default_substitutions = {
        'ε' : "ν*(ux*ux + uy*uy + uz*uz + vx*vx + vy*vy + vz*vz + wx*wx + wy*wy + wz*wz)",
    'ε_sgs' : "-u*(Nu_sgs+Lu_sgs) - v*(Nv_sgs+Lv_sgs) - w*(Nw_sgs+Lw_sgs)",
        'χ' : "κ*(bx*bx + by*by + bz*bz)",
    'χ_sgs' : "-b*(Nb_sgs+Lb_sgs)"
}

class BoussinesqChannelFlow(ChannelFlow):
    """
    Boussinesq flow in a channel with optional rotation.

    Args
    ----
        nx : (int)
            Grid resolution in :math:`x`

        ny : (int)
            Grid resolution in :math:`y`

        nz : (int)
            Grid resolution in :math:`z`

        Lx : (float)
            Domain extent in :math:`x`

        Ly : (float)
            Domain extent in :math:`y`

        Lz : (float)
            Domain extent in :math:`z`

        f : (float)
            Coriolis parameter

        ν : (float)
            'Molecular' viscosity

        κ : (float)
            'Molecular' diffusivity for buoyancy

        Nsq : (float)
            Background buoyancy gradient

        closure : (None or closure.EddyViscosityClosure)
            Turbulent closure for Large Eddy Simulation

        xleft : (float)
            Domain x-origin

        yleft : (float)
            Domain y-origin

        zbottom : (float)
            Domain z-origin

        **params : (any)
            Additional parameters to be added to the dedalus problem.
    """
    def __init__(self,
        nx = 32,
        ny = 32,
        nz = 32,
        Lx = 1.0,
        Ly = 1.0,
        Lz = 1.0,
        f = 0.0,
        κ = 1.43e-7,
        ν = 1.05e-6,
        Nsq = 0.0,
        closure = None,
        xleft = None,
        yleft = None,
        zbottom = None,
        substitutions = default_substitutions,
        **params
        ):

        ChannelFlow.__init__(self, nx, ny, nz, Lx, Ly, Lz, xleft, yleft, zbottom)

        # Create bases and domain
        self.xbasis = xbasis = de.Fourier('x', nx, interval=self.xlimits, dealias=3/2)
        self.ybasis = ybasis = de.Fourier('y', ny, interval=self.ylimits, dealias=3/2)
        self.zbasis = zbasis = de.Chebyshev('z', nz, interval=self.zlimits, dealias=3/2)
        self.domain = domain = de.Domain([xbasis, ybasis, zbasis], grid_dtype=np.float64)

        variables = ['p', 'b', 'u', 'v', 'w', 'bz', 'uz', 'vz', 'wz']
        add_closure_variables(variables, closure)

        self.problem = problem = de.IVP(domain, variables=variables, time='t')

        add_parameters(problem, f=f, ν=ν, κ=κ, Nsq=Nsq, **params)
        bind_parameters(self, f=f, ν=ν, κ=κ, Nsq=Nsq, **params)

        add_first_derivative_substitutions(problem, ['u', 'v', 'w', 'b'], ['x', 'y'])

        # LES Closure
        add_closure_substitutions(problem, closure, tracers=['b'])
        add_closure_equations(problem, closure, tracers=['b'])

        # Custom substitutions
        add_substitutions(problem, substitutions)

        # Equations
        problem.add_equation(f"dt(u) - ν*(dx(ux) + dy(uy) + dz(uz)) + dx(p) - f*v - Lu_sgs = - u*ux - v*uy - w*uz + Nu_sgs")
        problem.add_equation(f"dt(v) - ν*(dx(vx) + dy(vy) + dz(vz)) + dy(p) + f*u - Lv_sgs = - u*vx - v*vy - w*vz + Nv_sgs")
        problem.add_equation(f"dt(w) - ν*(dx(wx) + dy(wy) + dz(wz)) + dz(p) - b   - Lw_sgs = - u*wx - v*wy - w*wz + Nw_sgs")
        problem.add_equation(f"dt(b) - κ*(dx(bx) + dy(by) + dz(bz)) + Nsq*w       - Lb_sgs = - u*bx - v*by - w*bz + Nb_sgs")
        problem.add_equation("ux + vy + wz = 0")

        # First-order equivalencies
        problem.add_equation("bz - dz(b) = 0")
        problem.add_equation("uz - dz(u) = 0")
        problem.add_equation("vz - dz(v) = 0")
        problem.add_equation("wz - dz(w) = 0")

        # Pressure gauge condition
        problem.add_bc("right(p) = 0", condition="(nx == 0) and (ny == 0)")

        self.x = domain.grid(0)
        self.y = domain.grid(1)
        self.z = domain.grid(2)

    def set_noflux_bc_top(self):
        """Set a no flux condition on buoyancy on the top wall."""
        ChannelFlow.set_tracer_noflux_bc(self, "top", tracers="b")

    def set_noflux_bc_bottom(self):
        """Set a no flux condition on buoyancy on the bottom wall."""
        ChannelFlow.set_tracer_noflux_bc(self, "bottom", tracers="b")

    def build_solver(self, timestepper='RK443'):
        """Build a BoussinesqChanelFlow solver."""
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
        """Print messages."""
        logger.info('Iteration: %i, Time: %e, dt: %e' %(self.solver.iteration, self.solver.sim_time, dt))

        for name, task in self.log_tasks.items():
            logger.info("{} = {}".format(name, task(self)))
