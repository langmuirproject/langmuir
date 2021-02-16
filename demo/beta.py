from langmuir import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

elec = Electron()
eta = np.linspace(10,100,100)
geo = Cylinder(r=0.1*elec.debye, l=10*elec.debye)
I = finite_length_current(geo, elec, eta=eta, normalization='th')

def power_law(eta, c, beta):
    return c*eta**beta

popt, pcov = curve_fit(power_law, eta, I)

plt.plot(eta, I, label='Finite-length model')
plt.plot(eta, power_law(eta, *popt), ':k',
         label=r'Power law ($c={:.2f}, \beta={:.2f})$'.format(*popt))

plt.xlabel(r'$\eta$')
plt.ylabel(r'$I/I_\mathrm{th}$')
plt.legend()
plt.show()
