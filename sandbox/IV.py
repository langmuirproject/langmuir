#!/usr/bin/env python
from langmuir import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def power_law(eta, c, beta):
    return -c*eta**beta

elec = Species(n=35e10, eV=0.08)
for l in np.array([6, 10, 25, 50, 100, 500, 1000])*1e-3:
    geo = Cylinder(r=0.1*elec.debye, l=l)
    eta = np.linspace(10,100,100)
    I = -finite_length_current(geo, elec, eta=-eta)*1e6
    popt, pcov = curve_fit(power_law, eta, I)
    print(popt)
    plt.plot(eta,I)
plt.grid()
plt.show()
