import numpy as np
from mpi4py import MPI

def random_noise(domain, amplitude=1, seed=23):
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

def add_first_derivative_substitutions(problem, coordinate='x', variables=[]):
    """
    Add first-derivative substitutions for `variables` to problem.

    Args
    ----
        coordinate : (str)
            Coordinate along which the derivative is taken

        variables : (list)
            List of variables to make subsitutions for

    Example
    -------
    >> add_first_derivative_substitutions(problem, coordinate='x', variables=['u'])

    is equivalent to

    >> problem.substitutions['ux'] = "dx(u)"  
    """
    for var in variables:
        problem.substitutions[var+coordinate] = f"d{coordinate}({var})"

def bind_parameters(obj, **params):
    """
    Bind the name, value pairs in `params` as members of the class `obj`.
    """
    for name, value in params.items():
        setattr(obj, name, value)
