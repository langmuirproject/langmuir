from langmuir import *
from scipy.optimize import least_squares

elec = Species(n=120e10, T=1000)
geo = Cylinder(r=5*elec.debye, l=25e-3)
V0 = -0.2
V = np.array([0.9, 1.1])

I = finite_radius_current(geo, elec, V+V0)
print("Probe currents:", I)

def residual(x):
    n, V0 = x
    return (finite_radius_current(geo, Species(n=n*1e10, T=1000), V+V0) - I)*1e6

x0 = [10, -0.3]

res = least_squares(residual, x0, bounds=([10, -1],[200, 0]))
n, V0 = res.x

print("n:", n)
print("V0:", V0)
