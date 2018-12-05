import sys; sys.path.append("..")

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

import dedaLES

# Test basic stuff
closure = dedaLES.ConstantSmagorinsky()
model = dedaLES.BoussinesqChannelFlow(Lx=2*pi, Ly=2*pi, Lz=2*pi, nx=32, ny=32, nz=32, ν=1, κ=1, closure=closure)
model.set_default_bcs()
model.build_solver()
model.run(iterations=10, logcadence=1)
