from langmuir import *
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt

eta = np.linspace(-25, 25, 100)
geo = Sphere(1.5e-3)
eta = -20

kappas = np.linspace(4, 10, 200)
I_OML = np.zeros_like(kappas)
I1 = np.zeros_like(kappas)
I2 = np.zeros_like(kappas)

for i, kappa in enumerate(kappas):
    sp = Species('kappa', n=1e11, T=1000, kappa=kappa)
    I_OML[i] = OML_current(geo, sp, eta=eta)
    I1[i] = finite_radius_current(geo, sp, eta=eta, table='darian-marholm')
    I2[i] = finite_radius_current(geo, sp, eta=eta, table='darian-marholm', flip_kappa=True)

plt.plot(kappas, -I_OML, label='OML')
plt.plot(kappas, -I1, label='kappa')
plt.plot(kappas, -I2, label='1/kappa')
plt.grid()
plt.legend()
plt.show()

# I_OML = OML_current(geo, sp, eta=eta)
# I1 = finite_radius_current(geo, sp, eta=eta, table='darian-marholm', flip_kappa=False)
# # I2 = finite_radius_current(geo, sp, eta=eta, table='darian-marholm', flip_kappa=True)

# plt.plot(-eta, -I_OML, label='OML')
# plt.plot(-eta, -I1, label='kappa')
# # plt.plot(-eta, -I2, '.', label='1/kappa')
# plt.legend()
# plt.grid()
# plt.show()
