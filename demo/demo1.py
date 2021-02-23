from langmuir import *
import numpy as np
import matplotlib.pyplot as plt

V = np.linspace(-2, 2, 100)
I = OML_current(Cylinder(r=1e-3, l=25e-3), Electron(n=1e11, T=1000), V)

plt.plot(V, -I*1e6)
plt.xlabel('V [V]')
plt.ylabel('I [uA]')
plt.show()
