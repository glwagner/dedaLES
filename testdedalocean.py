import sys
import matplotlib.pyplot as plt
import numpy as np
import dedalocean

model = dedalocean.DeepConvectionModel()
model.set_random_ic()
model.run()
