from langmuir import *
from scipy.optimize import leastsq

geo = Cylinder(1e-3, 25e-3)
V0 = -0.5
V = np.array([2.5, 4.0, 5.5, 7.0])
I = OML_current(geo, Species(n=120e10, T=1000), V+V0)

def residual(x):
    n, V0 = x
    return OML_current(geo, Species(n=n, T=1500), V+V0) - I

x0 = [10e10, -0.3]
x, c = leastsq(residual, x0)
n, V0 = x
print(n)
print(V0)
