"""
closures.py

Classes and functions for implementing turbulence closures
in dedaLES.

Turbulence closures are generalized to consist of 'nonlinear' and
'linear' momentum stress and tracer flux divergences.
For u-momentum, for example, the nonlinear stress divergence is
denoted Nu_sgs, while the linear stress divergence is denoted Lu_sgs.
For a tracer c, the nonlinear flux divergence is Nc_sgs, and the linear
flux divergence is Lc_sgs.
"""

import numpy as np

from .utils import bind_parameters

tensor_components_3d = ['xx', 'yy', 'zz', 'xy', 'yz', 'zx', 'yx', 'zy', 'xz']

def add_closure_variables(variables, closure):
    # variables = variables + closure.variables
    pass

def add_closure_substitutions(problem, closure, tracers=[], **kwargs):
    """
    Add closure substitutions to the problem.
    """
    if closure is None:
        # No parameterized subgrid fluxes
        problem.substitutions['Nu_sgs'] = "0"
        problem.substitutions['Nv_sgs'] = "0"
        problem.substitutions['Nw_sgs'] = "0"
        problem.substitutions['Lu_sgs'] = "0"
        problem.substitutions['Lv_sgs'] = "0"
        problem.substitutions['Lw_sgs'] = "0"
        # For convenience:
        problem.substitutions['ν_sgs'] = "0"
        for c in tracers:
            problem.substitutions[f'N{c}_sgs'] = "0"
            problem.substitutions[f'L{c}_sgs'] = "0"
            # For convenience:
            problem.substitutions[f'κ{c}_sgs'] = "0"
    else:
        # closure.equations(problem, tracers=tracers)
        closure.add_substitutions(problem, tracers=tracers, **kwargs)

def add_closure_equations(problem, closure, **kwargs):
    """
    Add closure equations to the problem.
    """
    pass


class EddyViscosityClosure():
    """
    Generic, three-dimensional LES closure based on an eddy viscosity and diffusivity.

    Notes
    -----
        * `ν_sgs` denotes the nonlinear, isotropic eddy viscosity.
        * `νij_sgs_const` denotes the constant eddy viscosity tensor with components `ij = ['xx', 'yy', 'xy', ...]`
        * `κc_sgs` denotes the nonlinear, isotropic eddy diffusivity for tracer `c`. 
        * `κij_c_sgs_const` denotes the constant eddy diffusivity tensor for tracer `c` and component `ij`. 
    """
    def __init__(self):
        pass

    def add_substitutions_strain_rate_tensor(self, problem, u='u', v='v', w='w'):
        """
        Defines the strain-rate tensor `Sij` in `problem`.

        For the the indicies of the strain rate tensor,
        we use the shorthand x=1, y=2, z=3 for clarity.

        Args
        ----
        problem : de.Problem
            dedalus problem to which the substitiutions are added

        u : str
            Name of the flow field in the x-direction

        v : str
            Name of the flow field in the y-direction

        w : str
            Name of the flow field in the z-direction
        """
        bind_parameters(self, u='u', v='v', w='w')

        problem.substitutions['Sxx'] = f"{u}x"
        problem.substitutions['Syy'] = f"{v}y"
        problem.substitutions['Szz'] = f"{w}z"
        problem.substitutions['Sxy'] = f"({u}y + {v}x) / 2"
        problem.substitutions['Syz'] = f"({v}z + {w}y) / 2"
        problem.substitutions['Szx'] = f"({w}x + {u}z) / 2"
        problem.substitutions['Syx'] = "Sxy"
        problem.substitutions['Szy'] = "Syz"
        problem.substitutions['Sxz'] = "Szx"
        problem.substitutions['tr_S2'] = (
            "Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz")

    def add_substitutions_subgrid_stress(self, problem):
        """Defines the subgrid stress tensor `τij` and subgrid stress
        divergence `(Nu_sgs, Nv_sgs, Nw_sgs)` in `problem`.
        """
        # Get constant viscosity parameters
        problem.parameters['ν_split'] = getattr(self, 'ν_split', 0)

        for ij in tensor_components_3d:
            ν_component = f'ν{ij}_sgs_const'
            problem.parameters[ν_component] = getattr(self, ν_component, 0)

        # Linear terms
        for ij in tensor_components_3d:
            problem.substitutions[f'τ{ij}_linear'] = f"2 * ν{ij}_sgs_const * S{ij}" # this will fail if Sij involves non-constant terms.

        # Linear subgrid momentum fluxes
        u, v, w = self.u, self.v, self.w
        problem.substitutions['Lu_sgs'] = f"dx(τxx_linear) + dy(τyx_linear) + dz(τzx_linear) + ν_split*(dx({u}x) + dy({u}y) + dz({u}z))" 
        problem.substitutions['Lv_sgs'] = f"dx(τxy_linear) + dy(τyy_linear) + dz(τzy_linear) + ν_split*(dx({v}x) + dy({v}y) + dz({v}z))" 
        problem.substitutions['Lw_sgs'] = f"dx(τxx_linear) + dy(τyx_linear) + dz(τzx_linear) + ν_split*(dx({w}x) + dy({w}y) + dz({w}z))" 

        # Subgrid stress proportional to eddy viscosity
        for ij in tensor_components_3d:
            problem.substitutions[f'τ{ij}'] = f"2 * (ν_sgs - ν_split) * S{ij}" 

        # Nonlinear subgrid momentum fluxes
        problem.substitutions['Nu_sgs'] = "dx(τxx) + dy(τyx) + dz(τzx)"
        problem.substitutions['Nv_sgs'] = "dx(τxy) + dy(τyy) + dz(τzy)"
        problem.substitutions['Nw_sgs'] = "dx(τxz) + dy(τyz) + dz(τzz)"


    def add_substitutions_subgrid_flux(self, problem, c):
        """Defines the subgrid tracer flux for `qcx` for the tracer `c` and
        the subgrid flux divergence `Nc_sgs` in `problem`.
        """
        # Retrieve constant viscosity parameters
        problem.parameters['κ_split'] = getattr(self, 'κ_split', 0)

        for ij in tensor_components_3d:
            κ_component = f'κ{ij}_{c}_sgs_const'
            problem.parameters[κ_component] = getattr(self, κ_component, 0)

        # Linear subgrid tracer fluxes
        split_diffusion = f"κ_split*(dx({c}x) + dy({c}y) + dz({c}z))"
        diagonal_diffusion = f"κxx_{c}_sgs_const * dx({c}x) + κyy_{c}_sgs_const * dy({c}y) + κzz_{c}_sgs_const * dz({c}z)"
        offdiagonal_diffusion = (
            f"   (κxy_{c}_sgs_const + κyx_{c}_sgs_const) * dx({c}y)"
            f" + (κxz_{c}_sgs_const + κzx_{c}_sgs_const) * dx({c}z)"
            f" + (κyz_{c}_sgs_const + κzy_{c}_sgs_const) * dy({c}z)"
        )

        Lc_sgs = f"L{c}_sgs"
        problem.substitutions[Lc_sgs] = f"{split_diffusion} + {diagonal_diffusion} + {offdiagonal_diffusion}"

        # Nonlinear subgrid tracer fluxes
        qx = f"q{c}x"
        qy = f"q{c}y"
        qz = f"q{c}z"

        # Diffusivity for c
        κ_sgs = f"κ{c}_sgs - κ_split"

        problem.substitutions[qx] = f"- {κ_sgs} * dx({c})"
        problem.substitutions[qy] = f"- {κ_sgs} * dy({c})"
        problem.substitutions[qz] = f"- {κ_sgs} * dz({c})"

        Nc_sgs = f"N{c}_sgs"
        problem.substitutions[Nc_sgs] = f"- dx({qx}) - dy({qy}) - dz({qz})"

        
    def add_closure_substitutions(self):
        pass


class ConstantSmagorinsky(EddyViscosityClosure):
    """
    Constant Smagorinsky closure for Large Eddy Simulation.

    Parameters
    ----------
    Δ_const : float
        Constant-filter size

    C : float
        Poincare constant for grid-relative filter

    Sc : float
        Turbulent Schmidt number (Sc = ν_sgs / κ_sgs)

    ν_split : float
        Viscosity splitting parameter. A viscous term with viscosity ν_split 
        is added to the LHS and subtracted from the RHS.

    Notes
    -----
    The subgrid viscosity is calculated with constant-filter and grid-relative
    filter sizes as

        ν_sgs = [Δ_const**2 + (C*Δ_grid)**2] * sqrt(2*tr(S**2))
        Δ_grid = (Δx * Δy * Δz)**(1/3)

    where the grid-based fitler lengths are determined from the grid spacing
    and dealiasing factors D as

        Δx[i] = Dx * 2 * Δx[i] = Dx * (x[i+1] - x[i-1])

    """
    def __init__(self, Δ_const=0, C=0.17, Sc=1, ν_split=0):
        self.Δ_const = Δ_const
        self.C = C
        self.Sc = Sc
        self.ν_split = ν_split
        self.κ_split = ν_split / Sc

    def add_substitutions(self, problem, tracers=[], u='u', v='v', w='w', **kwargs):
        # Construct grid-based filter field
        Δx = problem.domain.bases[0].dealias * 2 * problem.domain.grid_spacing(0)
        Δy = problem.domain.bases[1].dealias * 2 * problem.domain.grid_spacing(1)
        Δz = problem.domain.bases[2].dealias * 2 * problem.domain.grid_spacing(2)

        self.Δ = problem.domain.new_field()
        self.Δ['g'] = (Δx*Δy*Δz)**(1/3)

        # Add subgrid parameters to problem
        problem.parameters['Δ0'] = self.Δ_const
        problem.parameters['Δ'] = self.Δ
        problem.parameters['C_poin'] = self.C
        problem.parameters['Sc_sgs'] = self.Sc

        # Add subgrid substitutions to problem
        self.add_substitutions_strain_rate_tensor(problem, u=u, v=v, w=w)
        problem.substitutions['ν_sgs'] = "(Δ0**2 + (C_poin*Δ)**2) * sqrt(2*tr_S2)"
        self.add_substitutions_subgrid_stress(problem)

        # Add tracer terms to problem
        for c in tracers:
            κ_sgs = f"κ{c}_sgs"
            problem.substitutions[κ_sgs] = "ν_sgs / Sc_sgs"
            self.add_substitutions_subgrid_flux(problem, c)



class AnisotropicMinimumDissipation(EddyViscosityClosure):
    """
    Anisotropic minimum dissipation (AMD) turbulence closure for Large Eddy Simulation.

    Parameters
    ----------
    C : float
        Poincare constant

    stratified : bool
        Set to 'True' to use the stratified version of AMD

    ν_split : float
        Viscosity splitting parameter. A viscous term with viscosity ν_split 
        is added to the LHS and subtracted from the RHS.

    κ_split : float
        Diffusivity splitting parameter. A diffusive term with diffusivity κ_split 
        is added to the LHS and subtracted from the RHS.

    quasi_strain : float
        A 'regularization' parameter that prevents the AMD eddy viscosity from being NaN.
        Only use this with trivial initial condition that lead to an AMD eddy viscosity of ν_sgs = 0/0.
        quasi_strain has units of strain, and should be much smaller than resolved strain in non-quiescent 
        regions of the flow.

    quasi_gradient : float
        A 'regularization' parameter that prevents the AMD eddy diffusivity from being NaN.
        See quasi_strain for caveats.
    """
    def __init__(self, C=1/12, stratified=False, ν_split=0, κ_split=0, quasi_strain=0, quasi_gradient=0):
        self.C = C
        self.stratified = stratified
        self.ν_split = ν_split
        self.κ_split = κ_split
        self.quasi_strain = quasi_strain
        self.quasi_gradient = quasi_gradient

    def add_substitutions(self, problem, u='u', v='v', w='w', b='b', tracers=[], **kwargs):
        """Add substitutions associated with the Anisotropic
        Minimum Dissipation model to the `problem`.

        Args
        ----
            problem : dedalus.IVP
                The dedalus initial-value-problem to add the substitutions to.

            u : str
                The name of the x-velocity in `problem`

            v : str
                The name of the y-velocity in `problem`

            w : str
                The name of the w-velocity in `problem`

            b : str
                The name of the buoyancy field in the stratified `problem`.

            tracers : list
                A list of the names of tracers in `problem`.
        """
        # Construct grid-based filter field
        Δx = problem.domain.bases[0].dealias * 2 * problem.domain.grid_spacing(0)
        Δy = problem.domain.bases[1].dealias * 2 * problem.domain.grid_spacing(1)
        Δz = problem.domain.bases[2].dealias * 2 * problem.domain.grid_spacing(2)

        self.Δx = problem.domain.new_field()
        self.Δy = problem.domain.new_field()
        self.Δz = problem.domain.new_field()
        self.Δx['g'] = Δx
        self.Δy['g'] = Δy
        self.Δz['g'] = Δz

        # Add subgrid parameters to problem
        problem.parameters['Δx'] = self.Δx
        problem.parameters['Δy'] = self.Δy
        problem.parameters['Δz'] = self.Δz
        problem.parameters['C_poin'] = self.C
        problem.parameters['quasi_strain_sq'] = self.quasi_strain**2
        problem.parameters['quasi_gradient_sq'] = self.quasi_gradient**2

        # Add subgrid substitutions to problem
        self.add_substitutions_strain_rate_tensor(problem, u=u, v=v, w=w)

        # Soft max function
        problem.substitutions['zero_max(x)'] = "0.5 * (abs(x) + x)"

        # AMD substitutions
        problem.substitutions['tr_uij'] = (
            f"{u}x*{u}x + {u}y*{u}y + {u}z*{u}z +" +
            f"{v}x*{v}x + {v}y*{v}y + {v}z*{v}z +" +
            f"{w}x*{w}x + {w}y*{w}y + {w}z*{w}z"
        )

        problem.substitutions['uik_ujk_Sij'] = (
            f"Δx**2 * ({u}x*{u}x*Sxx + {v}x*{v}x*Syy + {w}x*{w}x*Szz + 2*{u}x*{v}x*Sxy + 2*{u}x*{w}x*Sxz + 2*{v}x*{w}x*Syz) + " +
            f"Δy**2 * ({u}y*{u}y*Sxx + {v}y*{v}y*Syy + {w}y*{w}y*Szz + 2*{u}y*{v}y*Sxy + 2*{u}y*{w}y*Sxz + 2*{v}y*{w}y*Syz) + " +
            f"Δz**2 * ({u}z*{u}z*Sxx + {v}z*{v}z*Syy + {w}z*{w}z*Szz + 2*{u}z*{v}z*Sxy + 2*{u}z*{w}z*Sxz + 2*{v}z*{w}z*Syz)"
        )

        if self.stratified:
            # This implementation assumes that `w` is aligned with gravity.
            problem.substitutions['wk_bk'] = f"Δx**2 * {w}x*{b}x + Δy**2 * {w}y*{b}y + Δz**2 * {w}z*{b}z"
        else:
            problem.substitutions['wk_bk'] = "0"

        problem.substitutions['ν_sgs'] = "zero_max(C_poin * (wk_bk - uik_ujk_Sij)) / (tr_uij + quasi_strain_sq)"

        self.add_substitutions_subgrid_stress(problem)

        for c in tracers:
            # mod_Dc = |∇c|²
            mod_Dc = f"mod_D{c}"
            problem.substitutions[mod_Dc] = f"{c}x**2 + {c}y**2 + {c}z**2"

            # Dc_dot_ui = ∇c • ∂ᵢu
            Dc_dot_ux = f"D{c}_dot_ux"
            Dc_dot_uy = f"D{c}_dot_uy"
            Dc_dot_uz = f"D{c}_dot_uz"

            problem.substitutions[Dc_dot_ux] = f"{u}x*{c}x + {v}x*{c}y + {w}x*{c}z"
            problem.substitutions[Dc_dot_uy] = f"{u}y*{c}x + {v}y*{c}y + {w}y*{c}z"
            problem.substitutions[Dc_dot_uz] = f"{u}z*{c}x + {v}z*{c}y + {w}z*{c}z"

            # uik_ck_ci = Δₖ² ∂ₖuᵢ ∂ₖc ∂ᵢc = Δₖ² ∂ₖc (∇c • ∂ₖu)
            uik_ck_ci = f"uik_{c}k_{c}i"
            problem.substitutions[uik_ck_ci] = (
                   f"Δx**2 * {c}x * {Dc_dot_ux}" +
                f" + Δy**2 * {c}y * {Dc_dot_uy}" +
                f" + Δz**2 * {c}z * {Dc_dot_uz}")

            # κ_sgs = -C^2 Δₖ² ∂ₖuᵢ ∂ₖc ∂ᵢc / |∇c|²
            κ_sgs = f"κ{c}_sgs"
            problem.substitutions[κ_sgs] = f"zero_max(-C_poin * {uik_ck_ci} / ({mod_Dc} + quasi_gradient_sq))"
            self.add_substitutions_subgrid_flux(problem, c)


class StratifiedAnisotropicMinimumDissipation(AnisotropicMinimumDissipation):
    """
    A stratification-aware Anisotropic minimum dissipation (AMD) turbulence closure for Large Eddy Simulation.

    Parameters
    ----------
    C : float
        Poincare constant

    ν_split : float
        Viscosity splitting parameter. A viscous term with viscosity ν_split 
        is added to the LHS and subtracted from the RHS.

    κ_split : float
        Diffusivity splitting parameter. A diffusive term with diffusivity κ_split 
        is added to the LHS and subtracted from the RHS.

    quasi_strain : float
        A 'regularization' parameter that prevents the AMD eddy viscosity from being NaN.
        Only use this with trivial initial conditions that lead to an AMD eddy viscosity of ν_sgs = 0/0.
        quasi_strain has units of strain, and should be much smaller than resolved strain in non-quiescent 
        regions of the flow.

    quasi_gradient : float
        A 'regularization' parameter that prevents the AMD eddy diffusivity from being NaN.
        See quasi_strain for caveats.
    """
    def __init__(self, C=1/12, ν_split=0, κ_split=0, quasi_strain=0, quasi_gradient=0):
        AnisotropicMinimumDissipation.__init__(self, stratified=True, C=C, ν_split=0, κ_split=0, quasi_strain=quasi_strain, 
                                               quasi_gradient=quasi_gradient)
