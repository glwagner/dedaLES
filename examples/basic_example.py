import sys; sys.path.append("..")

import dedaLES

from dedaLES import mpiprint

# Basics
mpiprint('')

# 1. DNS of Navier-Stokes
dns_model = dedaLES.NavierStokesTriplyPeriodicFlow(nx=4, ny=4, nz=4, ν=1)
dns_model.build_solver()

mpiprint("\nDNS Navier-Stokes model built!\n")

# 2. LES of Navier-Stokes
les_model = dedaLES.NavierStokesTriplyPeriodicFlow(nx=4, ny=4, nz=4, ν=1, closure=dedaLES.ConstantSmagorinsky())
les_model.build_solver()

mpiprint("\nLES Navier-Stokes model built!\n")

# 2. Boussinesq in a channel
bouss_model = dedaLES.BoussinesqChannelFlow(nx=4, ny=4, nz=4, ν=1, κ=0.1)

for bc in ["no penetration", "no slip", "no flux"]:
    bouss_model.set_bc(bc, "top", "bottom")

bouss_model.build_solver()

mpiprint("\nBoussinesq model built!\n")
