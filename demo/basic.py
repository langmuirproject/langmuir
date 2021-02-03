from langmuir import *
import numpy as np
import matplotlib.pyplot as plt

plasma = [
    Electron(n=1e11, T=1000),
    Hydrogen(n=1e11, T=1000)
    ]

geometry = Cylinder(r=1e-3, l=60e-3, lguard=True)

V = np.linspace(-2, 2, 100)

I_OML = OML_current(geometry, plasma, V)
I_FL = finite_length_current(geometry, plasma, V)

plt.plot(V, -I_OML*1e6, label='OML')
plt.plot(V, -I_FL*1e6, label='FL')
plt.xlabel('V [V]')
plt.ylabel('I [uA]')
plt.legend()
plt.show()
