import sys
import matplotlib.pyplot as plt
import numpy as np
import fishalus

model = fishalus.DeepConvectionModel()
model.set_random_ic()
model.run()
