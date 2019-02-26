import numpy as np
from mpi4py import MPI

from dedalus.extras.flow_tools import CFL
from dedalus.core.future import FutureField

def grid_stats(model, axis):
    """Compute minimum and maximum grid spacing of the `model` along `axis`."""
    Δmin = np.min(model.problem.domain.grid_spacing(axis))
    Δmax = np.max(model.problem.domain.grid_spacing(axis))
    return Δmin, Δmax

def min_spacing(model):
    """Compute the minimum grid spacing of the `model`."""
    Δmin = []
    for axis in range(model.problem.domain.dim):
        Δmin.append(model.problem.domain.grid_spacing(axis))
    return np.min(np.array(Δmin))

def random_noise(domain, amplitude=1, seed=23):
    """Generate `domain`-spanning `amplitude` random noise with random seed `seed`."""
    rand = np.random.RandomState(seed=seed)
    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    return amplitude*rand.standard_normal(gshape)[slices]

def mpiprint(msg):
    if MPI.COMM_WORLD.Get_rank() is 0:
        print(msg)

def add_parameters(problem, **params):
    """
    Add parameters to a dedalus problem programmatically.
    """
    for name, value in params.items():
        problem.parameters[name] = value

def add_substitutions(problem, **substitutions):
    """
    Add substitutions to a dedalus problem programmatically.
    """
    for name, value in substitutions.items():
        problem.substitutions[name] = value

def add_first_derivative_substitutions(problem, variables, dims):
    """
    Add first-derivative substitutions for `variables` to `problem`.

    Args
    ----
        variables : str or list
            List of variables to make subsitutions for

        dims : str or list
            Dimensions along which derivatives are taken

    Example
    -------
    >> add_first_derivative_substitutions(problem, 'u', 'x')

    is equivalent to

    >> problem.substitutions['ux'] = "dx(u)"
    """
    if type(dims) is str: dims = [dims]
    if type(variables) is str: variables = [variables]

    for x in dims:
        for u in variables:
            problem.substitutions[u+x] = f"d{x}({u})"

def bind_parameters(obj, **params):
    """
    Bind the name, value pairs in `params` as members of the class `obj`.
    """
    for name, value in params.items():
        setattr(obj, name, value)


class TimeStepGizmo(CFL):
    """ 
    Computes frequency-limited timestep from a set of frequencies, velocities, and diffusivities.

    Parameters
    ----------
    solver : solver object
        Problem solver

    initial_dt : float
        Initial timestep

    cadence : int, optional
        Iteration cadence for computing new timestep (default: 1)

    safety : float, optional
        Safety factor for scaling computed timestep (default: 1.)

    max_dt : float, optional
        Maximum allowable timestep (default: inf)

    min_dt : float, optional
        Minimum allowable timestep (default: 0.)

    max_change : float, optional
        Maximum fractional change between timesteps (default: inf)

    min_change : float, optional
        Minimum fractional change between timesteps (default: 0.)

    threshold : float, optional
        Fractional change threshold for changing timestep (default: 0.)

    Notes
    -----
    The new timestep is computed by summing across the provided frequencies
    for each grid point, and then reciprocating the maximum "total" frequency
    from the entire grid.

    """
    def __init__(self, *args, **kwargs):
        CFL.__init__(self, *args, **kwargs)

    def add_diffusivity(self, diffusivity):
        """Add an on-grid isotropic diffusivity."""
        diff = FutureField.parse(diffusivity, self.solver.evaluator.vars, self.solver.domain)    
        for axis in range(self.solver.domain.dim):
            freq = diff / self.grid_spacings[axis]**2 
            self.add_frequency(freq)
