from langmuir import *
import numpy as np
import matplotlib.pyplot as plt

# c = Cylinder(1e-3, 50e-3)
# p = []
# p.append(Species(n=1e11, T=1000))
# p.append(Species('proton', 'kappa', n=1e11, T=1000, kappa=4))
# # Vs = [0.5,1,1.5,2,2.5]
# eta = [1, 10, 20]
# Is = tabulated_current(c, p, eta=eta, table='laframboise')
# print(Is)
# Is = OML_current(c, p, eta=eta)
# print(Is)

p = []
p.append(Species(n=1e11, T=1000))
p.append(Species('proton', 'kappa', n=1e11, T=1000, kappa=4))
c = Cylinder(2*p[0].debye, 1)
eta = -np.linspace(0,25,100)

I_oml = OML_current(c, p, eta=eta, normalize=True)
I_laf = tabulated_current(c, p, eta=eta, normalize=True, table='laframboise')

plt.plot(np.abs(eta), I_oml, label='OML')
plt.plot(np.abs(eta), I_laf, label='Laframboise')
plt.legend()
plt.show()
