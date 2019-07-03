from langmuir import *
from scipy.optimize import fsolve, leastsq, curve_fit

geo = Cylinder(1e-3, 25e-3)
V0 = -0.5
V = np.array([2.5, 4.0, 5.5, 7.0])
I = OML_current(geo, Species(n=120e10, T=1000, kappa=5), V+V0)

def func(x):
    n, V0, kappa = x
    res =  I - OML_current(geo, Species(n=n, T=1000, kappa=kappa), V+V0)
    return res

x0 = [10e10, -0.3, 5]
x, pcov = leastsq(func, x0)
print(x[0], x[2])
