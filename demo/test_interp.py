from langmuir import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def power_law(eta, c, beta):
    return -c*eta**beta

elec = Species(n=1e11, T=1000)
geo = Cylinder(l=40e-3, r=1e-3, lguard=0.77)
geo = Cylinder(l=20e-3, r=1e-3, lguard=0.077)
geo = Cylinder(l=40e-3, r=1e-3, lguard=0.077)
# geo = Cylinder(l=5e-3, r=1e-3, lguard=100e-3)
geo = Cylinder(l=10*elec.debye, r=1e-3, lguard=200*elec.debye)
etas = np.linspace(10,100,100)

# Is = finite_length_current(geo, elec, eta=etas)
# plt.plot(etas, Is)
# plt.show()

Vs = np.linspace(0, 8, 100)
Is_I = -finite_length_current(geo, elec, eta=-etas)
Is_g = -finite_length_current(geo, elec, eta=-etas, interpolate='g')

popt_g, pcov_g = curve_fit(power_law, etas, Is_g)
Is_g_fit = power_law(etas, *popt_g)
print("g beta:", popt_g[1])

popt_I, pcov_I = curve_fit(power_law, etas, Is_I)
Is_I_fit = power_law(etas, *popt_I)
print("I beta:", popt_I[1])

plt.plot(etas, Is_I, label='I')
plt.plot(etas, Is_g, label='g')
plt.plot(etas, Is_g_fit, label='g_fit')
plt.plot(etas, Is_I_fit, label='I_fit')
plt.legend()
plt.show()

# TODO:
# Potentials <0 and >100
