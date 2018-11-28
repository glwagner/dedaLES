import matplotlib.pyplot as plt
import numpy as np
import fishalus

model = fishalus.RayleighBernard()
model.set_random_ic()
model.run()
