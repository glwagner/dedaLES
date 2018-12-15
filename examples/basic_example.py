import sys; sys.path.append("..")
import dedaLES

# Basics

# 1. Navier-Stokes
ns_model = dedaLES.NavierStokesTriplyPeriodicFlow(nx=4, ny=4, nz=4, ν=1)
ns_model.build_solver()

print("Navier-Stokes model built!")

# 2. Boussinesq in a channel
bouss_model = dedaLES.BoussinesqChannelFlow(nx=4, ny=4, nz=4, ν=1, κ=0.1)

for bc in ["nopenetration", "noslip", "noflux"]:
    bouss_model.set_bc(bc, "top", "bottom")

bouss_model.build_solver()


print("Boussinesq model built!")
