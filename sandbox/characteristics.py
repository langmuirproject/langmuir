#!/usr/bin/env python
import matplotlib as mp

mp.rc('text', usetex=True)
mp.rc('font', family='serif', size=8)
mp.rc('axes', titlesize=8)
mp.rc('axes', labelsize=8)
mp.rc('xtick', labelsize=8)
mp.rc('ytick', labelsize=8)
mp.rc('legend', fontsize=8)
mp.rc('figure', titlesize=8)

from langmuir import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def power_law(eta, c, beta):
    return -c*eta**beta

elec = Species(n=1, eV=1)
lambds = np.array([1, 3, 10, 30, 600])

figwidth = 85
figheight = 70
dpi = 300

fig, ax = plt.subplots(figsize=(figwidth/25.4, figheight/25.4),
                       constrained_layout=True, dpi=dpi)

markers = list('o^spd')

eta = np.linspace(10,100,100)
eta_m = np.array([30, 50, 70, 90])

lw = 1.2
ms = 5

for i, lambd in enumerate(lambds):

    r=0.1*elec.debye
    l=lambd*elec.debye
    geo = Cylinder(r=r, l=l)

    I = finite_length_current(geo, elec, eta=-eta, normalization='th')
    I_m = finite_length_current(geo, elec, eta=-eta_m, normalization='th')

    popt, pcov = curve_fit(power_law, eta, I)

    line, = ax.plot(eta, I, '-', linewidth=lw, markersize=3, marker=markers[i],
                    markevery=(1000,1000),
                    label='$\lambda={:.0f}, \\beta={:.2f}$'.format(lambd, popt[1]))
    ax.plot(eta, power_law(eta, *popt), ':', color='black', linewidth=lw)
    ax.plot(eta_m, I_m, markers[i], markersize=ms, color=line.get_color())

ax.set_ylabel('$I/I_\mathrm{th}$')
ax.set_xlabel('$\eta$')
ax.legend()

fig.savefig('characteristics.png')
plt.show()
