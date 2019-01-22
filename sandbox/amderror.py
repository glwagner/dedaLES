import sys; sys.path.append("..")
import logging
import numpy as np
from dedalus.extras import flow_tools
import dedaLES

logger = logging.getLogger(__name__)

a = 10.0
dt = 1e-4
nx = nz = 32
ny = 4

closure = dedaLES.AnisotropicMinimumDissipation()
#closure = dedaLES.ConstantSmagorinsky()
#closure = None
model = dedaLES.BoussinesqChannelFlow(nx=nx, ny=ny, nz=nz, ν=1, κ=1, closure=closure)

model.set_bc("no penetration", "top", "bottom")
model.set_bc("no slip", "top", "bottom")
model.set_bc("no flux", "top", "bottom")
model.build_solver()

# Initial condition
fields = {}
for field in ['u', 'v', 'w', 'b']:
    noise = dedaLES.random_noise(model.domain)
    pert = a * noise * model.z * (1 - model.z)
    fields[field] = pert

model.set_fields(**fields)
model.stop_at(iteration=3)

# Flow properties
flow = flow_tools.GlobalFlowProperty(model.solver)
flow.add_property("sqrt(u*u + v*v + w*w) / ν", name='Re_domain')

model.solver.step(dt)
logger.info('Max domain Re = %e' %flow.max("Re_domain"))

model.solver.step(dt)
logger.info('Max domain Re = %e' %flow.max("Re_domain"))
