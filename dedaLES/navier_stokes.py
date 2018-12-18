import numpy as np

from numpy import pi
from dedalus import public as de

from .flows import Flow
from .utils import add_parameters, bind_parameters, add_first_derivative_substitutions
from .closures import add_closure_substitutions, add_closure_variables, add_closure_equations


class NavierStokesTriplyPeriodicFlow(Flow):
    """
    Flow in a triply-periodic box with optional rotation.

    Parameters
    ----------
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

        closure : (None or closure.EddyViscosityClosure)
            Turbulent closure for Large Eddy Simulation 

        **params : (any)
            Additional parameters to be added to the dedalus problem.
    """
    def __init__(self,
        nx = 32,
        ny = 32,
        nz = 32,
        Lx = 2*pi,      # [m]
        Ly = 2*pi,      # [m]
        Lz = 2*pi,      # [m]
        ν = 1.05e-6,    # Kinematic viscosity [m²/s]

        # Background fields:
        ub = "0",
        vb = "0",
        wb = "0",
        pb = "0",

        closure = None,  # Subgrid closure
        **params         # Additional parameters to be added to dedalus problem
        ):

        Flow.__init__(self, nx, ny, nz, Lx, Ly, Lz)

        self.xlimits = (-Lx/2, Lx/2)
        self.ylimits = (-Ly/2, Ly/2)
        self.zlimits = (-Lz/2, Lz/2)

        # Create bases and domain
        self.xbasis = xbasis = de.Fourier('x', nx, interval=self.xlimits, dealias=3/2)
        self.ybasis = ybasis = de.Fourier('y', ny, interval=self.ylimits, dealias=3/2)
        self.zbasis = zbasis = de.Fourier('z', nz, interval=self.zlimits, dealias=3/2)
        self.domain = domain = de.Domain([xbasis, ybasis, zbasis], grid_dtype=np.float64)

        self.x = domain.grid(0)
        self.y = domain.grid(1)
        self.z = domain.grid(2)

        bind_parameters(self, ν=ν, **params)

        # Problem set-up
        variables = ['p', 'u', 'v', 'w']
        add_closure_variables(variables, closure)

        self.problem = problem = de.IVP(domain, variables=variables, time='t')

        add_parameters(problem, ν=ν, **params)
        bind_parameters(self, ν=ν, **params)

        problem.substitutions['ub'] = ub
        problem.substitutions['vb'] = vb
        problem.substitutions['wb'] = wb
        problem.substitutions['pb'] = pb

        problem.substitutions['U'] = "u + ub"
        problem.substitutions['V'] = "v + vb"
        problem.substitutions['W'] = "w + wb"

        u = ['u', 'v', 'w', 'p']
        add_first_derivative_substitutions(problem, coordinate='x', variables=u)
        add_first_derivative_substitutions(problem, coordinate='y', variables=u)
        add_first_derivative_substitutions(problem, coordinate='z', variables=u)

        U = ['U', 'V', 'W']
        add_first_derivative_substitutions(problem, coordinate='x', variables=U)
        add_first_derivative_substitutions(problem, coordinate='y', variables=U)
        add_first_derivative_substitutions(problem, coordinate='z', variables=U)

        ub = ['ub', 'vb', 'wb', 'pb']
        add_first_derivative_substitutions(problem, coordinate='x', variables=ub)
        add_first_derivative_substitutions(problem, coordinate='y', variables=ub)
        add_first_derivative_substitutions(problem, coordinate='z', variables=ub)

        add_closure_substitutions(problem, closure)
        add_closure_equations(problem, closure)

        problem.substitutions['div(f1, f2, f3)'] = "dx(f1) + dy(f2) + dz(f3)"

        ## Momentum equations:
        # Primary linear terms
        linear_x = "px - ν*div(ux, uy, uz)" 
        linear_y = "py - ν*div(vx, vy, vz)" 
        linear_z = "pz - ν*div(wx, wy, wz)" 
        
        # Linear background terms
        linear_background_x = "pbx - ν*div(ubx, uby, ubz)" 
        linear_background_y = "pby - ν*div(vbx, vby, vbz)" 
        linear_background_z = "pbz - ν*div(wbx, wby, wbz)" 

        xmom = f"dt(u) + {linear_x} = - {linear_background_x} - U*Ux - V*Uy - W*Uz + Fx_sgs"
        ymom = f"dt(v) + {linear_y} = - {linear_background_y} - U*Vx - V*Vy - W*Vz + Fy_sgs"
        zmom = f"dt(w) + {linear_z} = - {linear_background_z} - U*Wx - V*Wy - W*Wz + Fz_sgs"
                   
        problem.add_equation(xmom)
        problem.add_equation(ymom)
        problem.add_equation(zmom)

        # Continuity equation
        problem.add_equation("ux + vy + wz = - ubx - vby - wbz",
                             condition="(nx != 0) or (ny != 0) or (nz != 0)")

        problem.add_equation("p = 0", condition="(nx == 0) and (ny == 0) and (nz == 0)")

    def build_solver(self, timestepper='RK443'):
        """
        Build a NavierStokesTriplyPeriodicFlow solver.
        """

        Flow.build_solver(self, timestepper=timestepper)

        self.u = self.solver.state['u']
        self.v = self.solver.state['v']
        self.w = self.solver.state['w']
