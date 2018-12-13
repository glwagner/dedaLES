import sys; sys.path.append("..")

import time, logging

import dedaLES


def benchmark(func):
    """
    A decorator for benchmarking.
    """
    def function_timer(*args, **kwargs):
        """
        A nested function for timing other functions.
        """
        start = time.time()
        value = func(*args, **kwargs)
        end = time.time()

        return end-start

    return function_timer

@benchmark
def benchmark_run(model, iterations=100, dt=1e-4):
    for i in range(iterations):
        model.solver.step(dt)

@benchmark
def benchmark_build(model):
    model.build_solver()

def init_rayleigh_benard_benchmark(nx=64, ny=64, nz=16, closure=None):
    Lx, Ly, Lz = 25.0, 25.0, 1.0 # Domain

    # Parameters
    Pr = 1.0       # Prandtl number
    f = 0.0        # Coriolis parameter
    κ = 1.0        # Thermal diffusivity 
    ε = 0.8        # Perturbation above criticality
    a = 1e-3       # Noise amplitude for initial condition

    # Constants
    Ra_critical = 1707.762
    Ra = Ra_critical + ε
    ν = Pr*κ                   # Viscosity 
    Bz = -Ra*Pr                 # Unstable buoyancy gradient

    # Construct model
    model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, 
                                          f=f, ν=ν, κ=κ, B0=Bz*Lz, a=a, closure=closure)

    model.set_bc("nopenetration", "top", "bottom")
    model.set_bc("freeslip", "top", "bottom")
    model.problem.add_bc("right(b) = B0")
    model.problem.add_bc("left(b) = 0")

    return model


def set_ic_rayleigh_benard_benchmark(model):
    noise = dedaLES.random_noise(model.domain)

    # Linear background + perturbations damped at walls
    zbottom, ztop = model.zbasis.interval
    a = model.problem.parameters['a']
    Bz = model.problem.parameters['B0']/(ztop-zbottom)
    z = model.z

    pert = a * noise * (z - ztop) * (z - zbottom)
    model.set_fields(b=Bz*(z - pert))
