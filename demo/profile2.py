from langmuir import *
import numpy as np
import matplotlib.pyplot as plt

elec = Electron(n=1e11, T=1000)

l_probe = 40e-3
l_guard = 25e-3
z_probe = np.linspace(0, l_probe, 100)
z_guard = np.linspace(-l_guard, 0, 100)

geo = Cylinder(r=0.25e-3, l=l_probe, lguard=True)
i_probe = finite_length_current_density(geo, elec, V=5, z=z_probe)
i_guard = finite_length_current_density(geo, elec, V=5, z=z_guard)

plt.plot(z_probe*1e3, -i_probe*1e6)
plt.plot(z_guard*1e3, -i_guard*1e6, ':C0')
plt.xlabel('z [mm]')
plt.ylabel('i [nA/mm]')
plt.show()
