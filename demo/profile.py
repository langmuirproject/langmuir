from langmuir import *
import numpy as np
import matplotlib.pyplot as plt

elec = Electron(n=1e11, T=1000)
l = 40e-3
z = np.linspace(0,l,100)

geo = Cylinder(r=0.25e-3, l=l)
i = finite_length_current_density(geo, elec, V=5, z=z)

plt.plot(z*1e3, -i*1e6)
plt.xlabel('z [mm]')
plt.ylabel('i [nA/mm]')
plt.show()
