from langmuir import *
from scipy.optimize import fsolve, leastsq, curve_fit

geo = Cylinder(1e-3, 25e-3)
V = np.array([2.5, 4.0, 5.5])-0.5
I = OML_current(geo, Species(n=30e10, T=1500), V)

def func(x):
    n, T, V0 = x
    return I - OML_current(geo, Species(n=n, T=T), V0+V)

x = fsolve(func, (10e10, 1000, -1), xtol=1e-10)
x = fsolve(func, (10e10, 1000, -1), xtol=1e-10)

# x, pcov = leastsq(func, (10e10, 1000, -1))
print(x)
print(func(x))
n, T, V0 = x
I2 = OML_current(geo, Species(n=n, T=T), V0+V)
print(I)
print(I2)
