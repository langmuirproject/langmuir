from langmuir import *
import numpy as np
import matplotlib.pyplot as plt

elec = Electron()

lambds = np.array([10, 30, 100])
eta = np.linspace(0, 100, 100)

for lambd in lambds:
    geometry = Cylinder(r=0.2*elec.debye, l=lambd*elec.debye)
    current = finite_length_current(geometry, elec, eta=eta, normalization='th')
    plt.plot(eta, current, label='$\lambda={}$'.format(lambd))

current = OML_current(geometry, elec, eta=eta, normalization='th')
plt.plot(eta, current, label='OML')

plt.xlabel('$\eta$')
plt.ylabel('$I/I_\mathrm{th}$')
plt.legend()
plt.show()
