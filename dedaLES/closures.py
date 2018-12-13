class LESClosure():
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
        qx = "q%sx" %tracer
        qy = "q%sy" %tracer
        qz = "q%sz" %tracer

        problem.substitutions[qx] = "- κ_sgs * dx(%s)" %tracer
        problem.substitutions[qy] = "- κ_sgs * dy(%s)" %tracer
        problem.substitutions[qz] = "- κ_sgs * dz(%s)" %tracer
        problem.substitutions['F%s_sgs' %tracer] = "- dx(%s) - dy(%s) - dz(%s)" %(qx, qy, qz)

    def add_closure_substitutions(self):
        pass


class ConstantSmagorinsky(LESClosure):
    """
    Constant Smagorinsky closure for Large Eddy Simulation.

    Parameters
    ----------
    δ_const : float
        Constant-filter size
    Cs : float
        Poincare constant for grid-relative filter

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
    def __init__(self, δ_const=0, Cs=0.17):
        self.δ_const = δ_const
        self.Cs = Cs

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

        # Add subgrid substitutions to problem
        self.substitute_strainratetensor(problem)
        problem.substitutions['ν_sgs'] = "(δ0**2 + (Cs*δ)**2) * sqrt(2*tr_S2)"
        self.substitute_subgridstress(problem)

        # Add tracer terms to problem
        problem.substitutions['κ_sgs'] = "ν_sgs"
        for tracer in tracers:
            self.substitute_subgridflux(problem, tracer)
