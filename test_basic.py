import sys
import matplotlib.pyplot as plt
import numpy as np
import dedaLES

# Test basic stuff
closure = dedaLES.ConstantSmagorinsky()
model = dedaLES.BoussinesqModel(closure=closure)
model.set_default_bcs()
model.build_solver()
model.run()
