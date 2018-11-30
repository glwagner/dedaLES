import sys
import matplotlib.pyplot as plt
import numpy as np
import dedaLES

model = dedaLES.OceanModel()
model.set_default_bcs()
model.build_solver()
model.run()
