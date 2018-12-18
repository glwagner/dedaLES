def add_closure_variables(variables, closure):
    # variables = variables + closure.variables
    pass

def add_closure_substitutions(problem, closure, tracers=[]):
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
        closure.add_substitutions(problem, tracers=tracers)


def add_closure_equations(problem, closure, tracers=[]):
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
        """
        Defines the subgrid stress tensor `τij` and subgrid stress
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

    def add_substitutions_subgrid_flux(self, problem, tracer):
        """
        Defines the subgrid tracer flux for `qcx` for `tracer=c` and
        the subgrid flux divergence `Fc_sgs` in `problem`.
        """
        # Subgrid buoyancy fluxes
        qx = f"q{tracer}x"
        qy = f"q{tracer}y"
        qz = f"q{tracer}z"

        # Diffusivity for tracer
        κ_sgs = f"κ{tracer}_sgs"

        problem.substitutions[qx] = f"- {κ_sgs} * dx({tracer})" 
        problem.substitutions[qy] = f"- {κ_sgs} * dy({tracer})"
        problem.substitutions[qz] = f"- {κ_sgs} * dz({tracer})"
        problem.substitutions[f"F{tracer}_sgs"] = f"- dx({qx}) - dy({qy}) - dz({qz})"

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

    def add_substitutions(self, problem, tracers=[]):
        # Construct grid-based filter field
        Δx = problem.domain.bases[0].dealias * 2 * problem.domain.grid_spacing(0)
        Δy = problem.domain.bases[1].dealias * 2 * problem.domain.grid_spacing(1)
        Δz = problem.domain.bases[2].dealias * 2 * problem.domain.grid_spacing(2)

        self.Δ = problem.domain.new_field()
        self.Δ['g'] = (Δx*Δy*Δz)**(1/3)

        # Add subgrid parameters to problem
        problem.parameters['Δ0'] = self.Δ_const
        problem.parameters['Δ'] = self.Δ
        problem.parameters['C_poincare'] = self.C
        problem.parameters['Sc_sgs'] = self.Sc

        # Add subgrid substitutions to problem
        self.add_substitutions_strain_rate_tensor(problem)
        problem.substitutions['ν_sgs'] = "(Δ0**2 + (C_poincare*Δ)**2) * sqrt(2*tr_S2)"
        self.add_substitutions_subgrid_stress(problem)

        # Add tracer terms to problem
        for tracer in tracers:
            κ_sgs = f"κ{tracer}_sgs"
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
    def __init__(self, C=0.2887, stratified=False):
        self.C = C
        self.stratified = stratified

    def add_substitutions(self, problem, tracers=[]):
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
        problem.parameters['C_poincare'] = self.C

        # Add subgrid substitutions to problem
        self.add_substitutions_strain_rate_tensor(problem, u=u, v=v, w=w)

        # Soft max function
        problem.substitutions['zero_max(x)'] = "0.5 * (abs(x) + x)"

        # AMD substitutions
        problem.substitutions['tr_uij'] = (
                "ux*ux + uy*uy + uz*uz + vx*vx + vy*vy + vz*vz + wx*wx + wy*wy + wz*wz")

        problem.substitutions['uik_ujk_Sij'] = (
                "   Δx**2 * (ux*ux*Sxx + vx*vx*Syy + wx*wx*Szz + 2*ux*vx*Sxy + 2*ux*wx*Sxz + 2*vx*wx*Syz)" + 
                " + Δy**2 * (uy*uy*Sxx + vy*vy*Syy + wy*wy*Szz + 2*uy*vy*Sxy + 2*uy*wy*Sxz + 2*vy*wy*Syz)" + 
                " + Δz**2 * (uz*uz*Sxx + vz*vz*Syy + wz*wz*Szz + 2*uz*vz*Sxy + 2*uz*wz*Sxz + 2*vz*wz*Syz)")

        if self.stratified:
            problem.substitutions['wk_bk'] = "Δx**2*wx*bx + Δy**2*wy*by + Δz**2*wz*bz"
        else:
            problem.substitutions['wk_bk'] = "0"

        problem.substitutions['ν_sgs'] = "zero_max(-C_poincare**2 * (uik_ujk_Sij - wk_bk) / tr_uij)"

        self.add_substitutions_subgrid_stress(problem)

        for c in tracers:
            # mod_Dc = |∇c|²
            mod_Dc = f"mod_D{c}"
            problem.substitutions[mod_Dc] = f"{c}x**2 + {c}y**2 + {c}z**2"

            # Dc_dot_ui = ∇c • ∂ᵢu
            Dc_dot_ux = f"D{c}_dot_ux"
            Dc_dot_uy = f"D{c}_dot_uy"
            Dc_dot_uz = f"D{c}_dot_uz"

            problem.substitutions[Dc_dot_ux] = f"ux*{c}x + vx*{c}y + wx*{c}z"
            problem.substitutions[Dc_dot_uy] = f"uy*{c}x + vy*{c}y + wy*{c}z"
            problem.substitutions[Dc_dot_uz] = f"uz*{c}x + vz*{c}y + wz*{c}z"

            # uik_ck_ci = Δₖ² ∂ₖuᵢ ∂ₖc ∂ᵢc = Δₖ² ∂ₖc (∇c • ∂ₖu) 
            uik_ck_ci = f"uik_{c}k_{c}i"
            problem.substitutions[uik_ck_ci] = (
                   f"Δx**2 * {c}x * {Dc_dot_ux}" + 
                f" + Δy**2 * {c}y * {Dc_dot_uy}" + 
                f" + Δz**2 * {c}z * {Dc_dot_uz}")

            # κ_sgs = -C^2 Δₖ² ∂ₖuᵢ ∂ₖc ∂ᵢc / |∇c|²
            κ_sgs = f"κ{c}_sgs"
            problem.substitutions[κ_sgs] = f"zero_max(-C_poincare**2 * {uik_ck_ci} / {mod_Dc})"
            self.add_substitutions_subgrid_flux(problem, c)
