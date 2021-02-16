from langmuir import *
import numpy as np
import matplotlib.pyplot as plt

plasma = Electron(n=1e11, T=1000)
geometry = Sphere(r=5*debye(plasma))

V = np.linspace(0, 2, 100)

I_OML = OML_current(geometry, plasma, V)
I_FR = finite_radius_current(geometry, plasma, V)

plt.plot(V, -I_OML*1e6, label='OML')
plt.plot(V, -I_FR*1e6, label='FR')
plt.xlabel('V [V]')
plt.ylabel('I [uA]')
plt.legend()
plt.show()
