import sys; sys.path.append("..")

import matplotlib.pyplot as plt
import numpy as np

import dedaLES

model = dedaLES.RayleighBernardConvection()
model.set_unstable_ic()
model.run()
