import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from langmuir import *
from scipy.special import erf

normalization='OML'

elec = Species(n=35e10, eV=0.08)
geo = Cylinder(r=0.1*elec.debye, l=20*elec.debye)
z = np.linspace(0, 20, 1000)
I = finite_length_current_density(geo, elec, zeta=z, eta=-25, normalization=normalization)

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)
line, = plt.plot(z, I)
if normalization=='th':
    ax.axis([-5, 105, 0.5, 20])
else:
    ax.axis([-5, 105, 0.5, 4.5])
ax.grid()
plt.xlabel(r'$z/\lambda_D$')
plt.ylabel(r'$I/I_\mathrm{OML}$')

ax_len = plt.axes([0.1, 0.10, 0.8, 0.03])
ax_eta = plt.axes([0.1, 0.05, 0.8, 0.03])
sl_len = Slider(ax_len, r'$l/\lambda_D$', 1, 100, valinit=20)
sl_eta = Slider(ax_eta, r'$\eta$', -10, 100, valinit=25)

def update(val):
    eta = sl_eta.val
    l = sl_len.val
    z = np.linspace(0, l, 1000)
    geo = Cylinder(r=0.1*elec.debye, l=l*elec.debye)
    I = finite_length_current_density(geo, elec, zeta=z, eta=-eta, normalization=normalization)
    line.set_ydata(I)
    line.set_xdata(z)

sl_len.on_changed(update)
sl_eta.on_changed(update)

plt.show()
