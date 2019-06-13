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
lambds = np.array([1, 3, 10, 30, 500])

figwidth = 85
figheight = 85/1.618
dpi = 300

lambd_p = 30
lambd_l = 5
r=0.1*elec.debye
l=lambd_p*elec.debye
lguard=lambd_l*elec.debye

fig, ax = plt.subplots(figsize=(figwidth/25.4, figheight/25.4),
                       constrained_layout=True, dpi=dpi)

lw = 1.2
ms = 3
markers = list('o^spd')
alpha = 0.08
color = '#000000'

zeta = np.linspace(-lambd_l, lambd_p, 1000)
zeta_n = zeta[np.where(zeta<0)[0]]
zeta_p = zeta[np.where(zeta>=0)[0]]
zeta_m = np.arange(-2.5, lambd_p, 5)
eta = np.linspace(10, 100, 100)

# NO GUARD

geo = Cylinder(r=r, l=l)

I = finite_length_current(geo, elec, eta=-eta, normalize='th')
popt, pcov = curve_fit(power_law, eta, I)
beta = popt[1]

i = finite_length_current_density(geo, elec, zeta=zeta_p, eta=-10, normalize='OML')
ax.fill_between(zeta_p, i, color=color, alpha=alpha, lw=0)
line, = ax.plot(zeta_p, i, lw=lw,
                label='$\\lambda_g={:.0f}, \\beta={:.2f}$'.format(0, beta),
                marker=markers[0], markersize=ms, markevery=(1000,1000))
i = finite_length_current_density(geo, elec, zeta=zeta_m, eta=-10, normalize='OML')
ax.plot(zeta_m[1:], i[1:], markers[0], markersize=ms, color=line.get_color())

# FINITE GUARD

geo = Cylinder(r=r, l=l, lguard=lguard)

I = finite_length_current(geo, elec, eta=-eta, normalize='th')
popt, pcov = curve_fit(power_law, eta, I)
beta = popt[1]

i = finite_length_current_density(geo, elec, zeta=zeta_p, eta=-10, normalize='OML')
ax.fill_between(zeta_p, i, color=color, alpha=alpha, lw=0)
line, = ax.plot(zeta_p, i, lw=lw,
                label='$\\lambda_g={:.0f}, \\beta={:.2f}$'.format(lambd_l, beta),
                marker=markers[1], markersize=ms, markevery=(1000,1000))
i = finite_length_current_density(geo, elec, zeta=zeta_n, eta=-10, normalize='OML')
ax.plot(zeta_n, i, ':', color=line.get_color(), lw=lw)
i = finite_length_current_density(geo, elec, zeta=zeta_m, eta=-10, normalize='OML')
ax.plot(zeta_m, i, markers[1], markersize=ms, color=line.get_color())

# INFINITE GUARD

geo = Cylinder(r=r, l=l, lguard=float('inf'))

I = finite_length_current(geo, elec, eta=-eta, normalize='th')
popt, pcov = curve_fit(power_law, eta, I)
beta = popt[1]

i = finite_length_current_density(geo, elec, zeta=zeta_p, eta=-10, normalize='OML')
ax.fill_between(zeta_p, i, color=color, alpha=alpha, lw=0)
line, = ax.plot(zeta_p, i, lw=lw,
                label='$\\lambda_g=\\infty, \\beta={:.2f}$'.format(beta),
                marker=markers[2], markersize=ms, markevery=(1000,1000))
i = finite_length_current_density(geo, elec, zeta=zeta_n, eta=-10, normalize='OML')
ax.plot(zeta_n, i, ':', color=line.get_color(), lw=lw)
i = finite_length_current_density(geo, elec, zeta=zeta_m, eta=-10, normalize='OML')
ax.plot(zeta_m, i, markers[2], markersize=ms, color=line.get_color())

ax.set_ylim([0.95, 1.95])
ax.set_ylabel('$g$')
ax.set_xlabel('$\zeta-\lambda_g$')
ax.legend(bbox_to_anchor=(0.565, 0.6))

fig.savefig('guard.png')
plt.show()
