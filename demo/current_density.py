from langmuir import *
import matplotlib.pyplot as plt

elec = Species(n=35e10, eV=0.08)
geo = Cylinder(r=1e-4, l=25e-3)
z = np.linspace(0, geo.l, 1000)

I = finite_length_current_density(geo, elec, z, eta=-25, normalize=True)

plt.plot(z, I)
plt.grid()
plt.show()
