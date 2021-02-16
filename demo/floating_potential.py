from langmuir import *
from scipy.constants import value as constants
from scipy.optimize import root_scalar

n=1e11
T=1000
e = constants('elementary charge')
k = constants('Boltzmann constant')

plasma = [Electron(n=n, T=T),
          Hydrogen(n=n, T=T)]

geometry = Sphere(r=0.2*debye(plasma))

def res(V):
    return OML_current(geometry, plasma, V)

sol = root_scalar(res, x0=-3, x1=0)

print(sol.root)
print(-2.5*k*T/e)
