from langmuir import *
import numpy as np
import matplotlib.pyplot as plt

elec = Electron()

eta = np.linspace(0, 100, 100)

geometry = Cylinder(r=0.2*elec.debye, l=10*elec.debye)
current = OML_current(geometry, elec, eta=eta, normalization='th')

plt.plot(eta, current)
plt.xlabel('$\eta$')
plt.ylabel('$I/I_\mathrm{th}$')
plt.show()
