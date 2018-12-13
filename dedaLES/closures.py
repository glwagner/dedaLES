class EddyViscosityClosure():
    """
    Generic LES closure.
    """
    def __init__(self):
        pass

    def substitute_strainratetensor(self, problem):
        # Strain-rate tensor of resolved flow
        problem.substitutions['Sxx'] = "ux"
        problem.substitutions['Syy'] = "vy"
        problem.substitutions['Szz'] = "wz"
        problem.substitutions['Sxy'] = "(uy + vx) / 2"
        problem.substitutions['Syz'] = "(vz + wy) / 2"
        problem.substitutions['Szx'] = "(wx + uz) / 2"
        problem.substitutions['Syx'] = "Sxy"
        problem.substitutions['Szy'] = "Syz"
        problem.substitutions['Sxz'] = "Szx"
        problem.substitutions['tr_S2'] = (
            "Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz")

    def substitute_subgridstress(self, problem):
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

    def substitute_subgridflux(self, problem, tracer):
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
    δ_const : float
        Constant-filter size
    Cs : float
        Poincare constant for grid-relative filter
    Sc : float
        Turbulent Schmidt number (Sc = ν_sgs / κ_sgs)

    Notes
    -----
    The subgrid viscosity is calculated with constant-filter and grid-relative
    filter sizes as

        ν_sgs = [δ_const**2 + (Cs*δ_grid)**2] * sqrt(2*tr(S**2))
        δ_grid = (δx * δy * δz)**(1/3)

    where the grid-based fitler lengths are determined from the grid spacing
    and dealiasing factors D as

        δx[i] = Dx * 2 * Δx[i] = Dx * (x[i+1] - x[i-1])

    """
    def __init__(self, δ_const=0, Cs=0.17, Sc=1):
        self.δ_const = δ_const
        self.Cs = Cs
        self.Sc = Sc

    def substitutions(self, problem, tracers=[]):
        # Construct grid-based filter field
        δx = problem.domain.bases[0].dealias * 2 * problem.domain.grid_spacing(0)
        δy = problem.domain.bases[1].dealias * 2 * problem.domain.grid_spacing(1)
        δz = problem.domain.bases[2].dealias * 2 * problem.domain.grid_spacing(2)

        self.δ = problem.domain.new_field()
        self.δ['g'] = (δx*δy*δz)**(1/3)

        # Add subgrid parameters to problem
        problem.parameters['δ0'] = self.δ_const
        problem.parameters['δ'] = self.δ
        problem.parameters['Cs'] = self.Cs
        problem.parameters['Sc_sgs'] = self.Sc

        # Add subgrid substitutions to problem
        self.substitute_strainratetensor(problem)
        problem.substitutions['ν_sgs'] = "(δ0**2 + (Cs*δ)**2) * sqrt(2*tr_S2)"
        self.substitute_subgridstress(problem)

        # Add tracer terms to problem
        for tracer in tracers:
            κ_sgs = f"κ{tracer}_sgs"
            problem.substitutions[κ_sgs] = "ν_sgs / Sc_sgs"
            self.substitute_subgridflux(problem, tracer)



class AnisotropicMinimumDissipation(EddyViscosityClosure):
    """
    Anisotropic minimum dissipation turbulence closure for Large Eddy Simulation.

    Parameters
    ----------
    δ_const : float
        Constant-filter size
    C : float
        Poincare constant for grid-relative filter
    """
    def __init__(self, C=0.2887, stratified=False):
        self.C = C
        self.stratified = stratified

    def substitutions(self, problem, tracers=[]):
        # Construct grid-based filter field
        δx = problem.domain.bases[0].dealias * 2 * problem.domain.grid_spacing(0)
        δy = problem.domain.bases[1].dealias * 2 * problem.domain.grid_spacing(1)
        δz = problem.domain.bases[2].dealias * 2 * problem.domain.grid_spacing(2)

        self.δx = problem.domain.new_field()
        self.δy = problem.domain.new_field()
        self.δz = problem.domain.new_field()
        self.δx['g'] = δx
        self.δy['g'] = δy
        self.δz['g'] = δz

        # Add subgrid parameters to problem
        problem.parameters['δx'] = self.δx
        problem.parameters['δy'] = self.δy
        problem.parameters['δz'] = self.δz
        problem.parameters['C'] = self.C

        # Add subgrid substitutions to problem
        self.substitute_strainratetensor(problem)

        # AMD substitutions
        problem.substitutions['tr_uij'] = (
                "ux*ux + uy*uy + uz*uz + vx*vx + vy*vy + vz*vz + wx*wx + wy*wy + wz*wz")

        problem.substitutions['uik_ujk_Sij'] = (
                "   δx**2 * (ux*ux*Sxx + vx*vx*Syy + wx*wx*Szz + 2*ux*vx*Sxy + 2*ux*wx*Sxz + 2*vx*wx*Syz)" + 
                " + δy**2 * (uy*uy*Sxx + vy*vy*Syy + wy*wy*Szz + 2*uy*vy*Sxy + 2*uy*wy*Sxz + 2*vy*wy*Syz)" + 
                " + δz**2 * (uz*uz*Sxx + vz*vz*Syy + wz*wz*Szz + 2*uz*vz*Sxy + 2*uz*wz*Sxz + 2*vz*wz*Syz)")

        if self.stratified:
            problem.substitutions['wk_bk'] = "δx**2*wx*bx + δy**2*wy*by + δz**2*wz*bz"
        else:
            problem.substitutions['wk_bk'] = "0"

        problem.substitutions['ν_sgs'] = "-C**2 * (uik_ujk_Sij - wk_bk) / tr_uij"

        self.substitute_subgridstress(problem)

        for c in tracers:
            # mod_Dc = |∇c|²
            mod_Dc = f"mod_D{c}"
            problem.substitutions[mod_Dc] = f"{c}x**2 + {c}y**2 + {c}**2"

            # Dc_dot_ui = ∇c • ∂ᵢu
            Dc_dot_ux = f"n{c}_dot_ux"
            Dc_dot_uy = f"n{c}_dot_uy"
            Dc_dot_uz = f"n{c}_dot_uz"

            problem.substitutions[Dc_dot_ux] = f"ux*{c}x + vx*{c}y + wx*{c}z"
            problem.substitutions[Dc_dot_uy] = f"uy*{c}x + vy*{c}y + wy*{c}z"
            problem.substitutions[Dc_dot_uz] = f"uz*{c}x + vz*{c}y + wz*{c}z"

            # uik_ck_ci = Δₖ² ∂ₖuᵢ ∂ₖc ∂ᵢc = Δₖ² ∂ₖc (∇c • ∂ₖu) 
            uik_ck_ci = f"uik_{c}k_{c}i"
            problem.substitutions[uik_ck_ci] = (
                   f"δx**2 * {c}x * {Dc_dot_ux}" + 
                f" + δy**2 * {c}y * {Dc_dot_uy}" + 
                f" + δz**2 * {c}z * {Dc_dot_uz}")

            # κ_sgs = -C^2 Δₖ² ∂ₖuᵢ ∂ₖc ∂ᵢc / |∇c|²
            κ_sgs = f"κ{c}_sgs"
            problem.substitutions[κ_sgs] = f"-C**2 * {uik_ck_ci} / {mod_Dc}"
            self.substitute_subgridflux(problem, c)

























