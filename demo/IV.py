#!/usr/bin/env python
from langmuir import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def power_law(eta, c, beta):
    return -c*eta**beta

elec = Species(n=35e10, eV=0.08)
# lambds = np.array([2, 3.5, 5, 7.5, 10, 14, 20, 25, 30])*1e-3/elec.debye
lambds = np.array([2, 3.5, 5, 7.5, 10, 15, 20])*1e-3/elec.debye

for lambd in lambds:
    geo = Cylinder(r=0.1*elec.debye, l=lambd*elec.debye)
    eta = np.linspace(2,100,100)
    I = finite_length_current(geo, elec, eta=-eta, normalize='th')#*1e6
    popt, pcov = curve_fit(power_law, eta, I)
    print(popt)
    line, = plt.plot(eta, I, label='$\lambda={:.2f}, \\beta={:.2f}$'.format(lambd*elec.debye*1e3, popt[1]))
    plt.plot(eta, power_law(eta, *popt), '--', color=line.get_color())

plt.ylabel('$I/I_\mathrm{th}$')
plt.xlabel('$\eta$')
plt.title('$IV$-curves of finite-length probes and fits (dashed) to $I=cV^\\beta$')
plt.legend()
plt.grid()
plt.show()

elec = Species(n=35e10, eV=0.08)
l = 50e-3
eta = np.linspace(10,100,100)

geo = Cylinder(r=0.1*elec.debye, l=l)
I = -finite_length_current(geo, elec, eta=-eta)*1e6
popt, pcov = curve_fit(power_law, eta, I)
print(popt)
line, = plt.plot(eta, I, label='No guard, $\\beta={:.2f}$'.format(popt[1]))
plt.plot(eta, power_law(eta, *popt), '--', color=line.get_color())

geo = Cylinder(r=0.1*elec.debye, l=l, lguard=25e-3)
I = -finite_length_current(geo, elec, eta=-eta)*1e6
popt, pcov = curve_fit(power_law, eta, I)
print(popt)
line, = plt.plot(eta, I, label='25mm guard, $\\beta={:.2f}$'.format(popt[1]))
plt.plot(eta, power_law(eta, *popt), '--', color=line.get_color())

geo = Cylinder(r=0.1*elec.debye, l=l, lguard=500e-3)
I = -finite_length_current(geo, elec, eta=-eta)*1e6
popt, pcov = curve_fit(power_law, eta, I)
print(popt)
line, = plt.plot(eta, I, label='Ideal guard, $\\beta={:.2f}$'.format(popt[1]))
plt.plot(eta, power_law(eta, *popt), '--', color=line.get_color())

geo = Cylinder(r=0.1*elec.debye, l=l, lguard=500e-3, rguard=500e-3)
I = -finite_length_current(geo, elec, eta=-eta)*1e6
popt, pcov = curve_fit(power_law, eta, I)
print(popt)
line, = plt.plot(eta, I, label='Two ideal guards, $\\beta={:.2f}$'.format(popt[1]))
plt.plot(eta, power_law(eta, *popt), '--', color=line.get_color())

plt.ylabel('$I\,[\mathrm{\mu A}]$')
plt.xlabel('$\eta$')
plt.title('$IV$-curves of guarded probes and fits (dashed) to $I=cV^\\beta$.')
plt.legend()
plt.grid()
plt.show()

