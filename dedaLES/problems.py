import numpy as np
from numpy import pi
from mpi4py import MPI

from .closures import ConstantSmagorinsky
from .flows import BoussinesqChannelFlow


def tanhstep(z, z0, d):
    return 0.5*(np.tanh((z-z0)/d) + 1)


class StratifiedCouetteFlow(BoussinesqChannelFlow):
    def __init__(self,
        nx = 64,
        ny = 64,
        nz = 48,
        Lx = 4*pi,      # Domain length [m]
        Ly = 4*pi,      # Domain length [m]
        Lz = 2.0,       # Domain height [m]
        f = 0.0,        # Coriolis parameters [s⁻¹]
        κ = 7.0,        # Diffusivity [m²/s]
        ν = 1.0,        # Kinematic viscosity [m²/s]
        utop = 0.0,     # Buoyancy gradient [s⁻²]
        closure = None, # LES closure
        **params):

        BoussinesqChannelFlow.__init__(self, nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz, f=f, κ=κ, ν=ν, 
                                       closure=closure, **params)

        self.set_nopenetration_topandbottom()
        self.set_wall_velocity_top(utop=utop)
        self.set_noslip_bottom()
        self.set_noflux_topandbottom()

        self.build_solver()


class RayleighBernardConvection(BoussinesqChannelFlow):
    def __init__(self,
        nx = 128,
        ny = 128,
        nz = 16,
        Lx = 800.0,     # Domain length [m]
        Ly = 800.0,     # Domain length [m]
        Lz = 100.0,     # Domain height [m]
        f = 1e-4,       # Coriolis parameters [s⁻¹]
        κ = 1.43e-7,    # Diffusivity [m²/s]
        ν = 1e-6,       # Kinematic viscosity [m²/s]
        Ra = 2500,      # Rayleigh number
        Bz = None,      # Buoyancy gradient [s⁻²]
        closure = None, # LES closure
        **params):

        if Bz is None:
            Bz = Ra * ν * κ / (Lz**4)

        BoussinesqChannelFlow.__init__(self, nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz, f=f, κ=κ, ν=ν, 
                                       closure=closure, Bz=Bz, **params)
                                       

        self.set_nopenetration_topandbottom()
        self.set_noslip_topandbottom()
        self.problem.add_bc("left(b) = -left(Bz*z)")
        self.problem.add_bc("right(b) = -right(Bz*z)")

        self.build_solver()


    def set_unstable_ic(self, magnitude=1e-3):

        # Random perturbations, initialized globally for same results in parallel
        gshape = self.domain.dist.grid_layout.global_shape(scales=1)
        slices = self.domain.dist.grid_layout.slices(scales=1)
        rand = np.random.RandomState(seed=23)
        noise = rand.standard_normal(gshape)[slices]

        # Linear background + perturbations damped at walls
        zb, zt = self.zbasis.interval
        pert =  magnitude * noise * (zt - self.z) * (self.z - zb) / self.Lz
        self.b['g'] = -self.Bz*(self.z - pert)
        self.b.differentiate('z', out=self.bz)
