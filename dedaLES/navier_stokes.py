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

        # "Background" (non-evolving) fields:
        u_bg = "0",
        v_bg = "0",
        w_bg = "0",
        p_bg = "0",

        include_linear_bg = False, # Include linear background terms?

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

        problem.substitutions['u_bg'] = u_bg 
        problem.substitutions['v_bg'] = v_bg
        problem.substitutions['w_bg'] = w_bg
        problem.substitutions['p_bg'] = p_bg

        problem.substitutions['U'] = "u + u_bg"
        problem.substitutions['V'] = "v + v_bg"
        problem.substitutions['W'] = "w + w_bg"

        add_first_derivative_substitutions(problem, 
            ['u', 'v', 'w', 'p', 'U', 'V', 'W', 'u_bg', 'v_bg', 'w_bg', 'p_bg'], 
            ['x', 'y', 'z']
        )

        add_closure_substitutions(problem, closure, u='U', v='V', w='W')
        add_closure_equations(problem, closure, u='U', v='V', w='W')

        problem.substitutions['div(f1, f2, f3)'] = "dx(f1) + dy(f2) + dz(f3)"

        ## Momentum equations:
        # Primary linear terms
        linear_x = "px - ν*div(ux, uy, uz)" 
        linear_y = "py - ν*div(vx, vy, vz)" 
        linear_z = "pz - ν*div(wx, wy, wz)" 
        
        # Linear background terms
        if include_linear_bg_terms:
            linear_bg_x = "p_bgx - ν*div(u_bgx, u_bgy, u_bgz)" 
            linear_bg_y = "p_bgy - ν*div(v_bgx, v_bgy, v_bgz)" 
            linear_bg_z = "p_bgz - ν*div(w_bgx, w_bgy, w_bgz)" 
        else:
            linear_bg_x = "0"
            linear_bg_y = "0"
            linear_bg_z = "0"

        xmom = f"dt(u) + {linear_x} = - {linear_bg_x} - U*Ux - V*Uy - W*Uz + Fx_sgs"
        ymom = f"dt(v) + {linear_y} = - {linear_bg_y} - U*Vx - V*Vy - W*Vz + Fy_sgs"
        zmom = f"dt(w) + {linear_z} = - {linear_bg_z} - U*Wx - V*Wy - W*Wz + Fz_sgs"
                   
        problem.add_equation(xmom)
        problem.add_equation(ymom)
        problem.add_equation(zmom)

        # Continuity equation
        pro_bglem.add_equation("ux + vy + wz = - u_bgx - v_bgy - w_bgz",
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
