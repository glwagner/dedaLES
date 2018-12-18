import numpy as np

def add_closure_variables(variables, closure):
    # variables = variables + closure.variables
    pass

def add_closure_substitutions(problem, closure, tracers=[], **kwargs):
    """
    Add closure substitutions to the problem.
    """
    if closure is None:
        # No parameterized subgrid fluxes
        problem.substitutions['Fx_sgs'] = "0"
        problem.substitutions['Fy_sgs'] = "0"
        problem.substitutions['Fz_sgs'] = "0"
        for c in tracers:
            problem.substitutions[f'F{c}_sgs'] = "0"
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
    Generic LES closure based on an eddy viscosity and diffusivity.
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
        divergence `(Fx_sgs, Fy_sgs, Fz_sgs)` in `problem`.
        """
        # Subgrid stress proportional to eddy viscosity
        problem.substitutions['τxx'] = "2 * ν_sgs * Sxx"
        problem.substitutions['τxy'] = "2 * ν_sgs * Sxy"
        problem.substitutions['τxz'] = "2 * ν_sgs * Sxz"
        problem.substitutions['τyx'] = "2 * ν_sgs * Syx"
        problem.substitutions['τyy'] = "2 * ν_sgs * Syy"
        problem.substitutions['τyz'] = "2 * ν_sgs * Syz"
        problem.substitutions['τzx'] = "2 * ν_sgs * Szx"
        problem.substitutions['τzy'] = "2 * ν_sgs * Szy"
        problem.substitutions['τzz'] = "2 * ν_sgs * Szz"

        # Subgrid momentum fluxes
        problem.substitutions['Fx_sgs'] = "dx(τxx) + dy(τyx) + dz(τzx)"
        problem.substitutions['Fy_sgs'] = "dx(τxy) + dy(τyy) + dz(τzy)"
        problem.substitutions['Fz_sgs'] = "dx(τxz) + dy(τyz) + dz(τzz)"

    def add_substitutions_subgrid_flux(self, problem, c):
        """Defines the subgrid tracer flux for `qcx` for the tracer `c` and
        the subgrid flux divergence `Fc_sgs` in `problem`.
        """
        # Subgrid buoyancy fluxes
        qx = f"q{c}x"
        qy = f"q{c}y"
        qz = f"q{c}z"

        # Diffusivity for c
        κ_sgs = f"κ{c}_sgs"
        Fc_sgs = f"F{c}_sgs"

        problem.substitutions[qx] = f"- {κ_sgs} * dx({c})" 
        problem.substitutions[qy] = f"- {κ_sgs} * dy({c})"
        problem.substitutions[qz] = f"- {κ_sgs} * dz({c})"
        problem.substitutions[Fc_sgs] = f"- dx({qx}) - dy({qy}) - dz({qz})"

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
    def __init__(self, Δ_const=0, C=0.17, Sc=1):
        self.Δ_const = Δ_const
        self.C = C
        self.Sc = Sc

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
            self.add_substitutions_subgrid_flux(problem, tracer)



class AnisotropicMinimumDissipation(EddyViscosityClosure):
    """
    Anisotropic minimum dissipation (AMD) turbulence closure for Large Eddy Simulation.

    Parameters
    ----------
    C : float
        Poincare constant

    stratified : bool
        Set to 'True' to use the stratified version of AMD
    """
    def __init__(self, C=1.0/np.sqrt(12), stratified=False):
        self.C = C
        self.stratified = stratified

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

        problem.substitutions['ν_sgs'] = "zero_max(-C_poin**2 * (uik_ujk_Sij - wk_bk) / tr_uij)"

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
            problem.substitutions[κ_sgs] = f"zero_max(-C_poin**2 * {uik_ck_ci} / {mod_Dc})"
            self.add_substitutions_subgrid_flux(problem, c)
