import h5py
import numpy as np
import matplotlib.pyplot as plt

filename = "snapshots/snapshots_s10/snapshots_s10_p0.h5" #
f = h5py.File(filename, 'r')

# List all groups
n1 = len(list(f.keys()))
print(str(n1) +' clubs')
print()
clubs = list(f.keys())
print(clubs)
# Get the data
print('members of club 2')
members2 = list(f[clubs[1]])
print(list(members2))


print('members of club 1')
members1 = list(f[clubs[0]])
print(members1)

#output the values of all members
for j in members1:
    print(j)
    print(list(f[clubs[0]][j]))

#plot the horizontal component of velocity at the very last time-step
uv = np.array(f['tasks']['u'])
print(np.shape(uv))
plt.pcolormesh(uv[-1], cmap='RdBu_r')
plt.show()

