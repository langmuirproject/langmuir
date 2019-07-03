from langmuir import *
from scipy.optimize import least_squares

elec = Species(n=120e10, T=1000)        # Actual plasma electron parameters
geo = Cylinder(r=5*elec.debye, l=25e-3) # Probe geometry
V0 = -0.2                               # Actual floating potential
V = np.array([0.9, 1.1])                # Probe biases wrt. floating potential

I = finite_radius_current(geo, elec, V+V0)
print("Probe currents:", I)

# Residual function for a given (n, V0)-tuple. 
# Parameters are scaled to make all quantities sufficiently close to unity
# Temperature is near-singular for infinite cylindrical probes,
# so we make no attempt at solving the (n, T, V0)-tuple, but assume T is known.
def residual(x):
    n, V0 = x
    return (finite_radius_current(geo, Species(n=n*1e10, T=1000), V+V0) - I)*1e6

x0 = [10, -0.3] # Initial guess of (n, V0)

# Least-squares solution using Scipy
res = least_squares(residual, x0, bounds=([10, -1],[200, 0]))
n, V0 = res.x

print("n:", n)
print("V0:", V0)
