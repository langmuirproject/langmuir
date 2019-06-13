from langmuir import *
from scipy.optimize import leastsq

sp = Species(n=1e11, T=1000)
geo = Cylinder(1e-3, 25e-3)
I = -0.4e-6

def residual(V):
    return finite_radius_current(geo, sp, V) - I

x, c = leastsq(residual, 0)
print(x[0])
