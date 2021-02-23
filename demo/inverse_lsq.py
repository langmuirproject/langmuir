from langmuir import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from tqdm import tqdm

geometry = Cylinder(r=0.255e-3, l=25e-3, lguard=True)
V = np.array([2,3,4,5])

model_truth = OML_current
model_pred = OML_current

data = generate_synthetic_data(geometry, V, model=model_truth)

alt = data['alt']
ne_jacobsen = jacobsen_density(geometry, V, data['I'])
ne_fit = np.zeros_like(alt)
V0_fit = np.zeros_like(alt)

def residual(x, I, T):
    n, V0 = x
    return (model_pred(geometry, Electron(n=n, T=T), V+V0) - I)*1e6

x0 = [10e10, -0.3]
x_scale = [1e10, 1]
for i in tqdm(range(len(alt))):
    I = data['I'][i]
    T = data['Te'][i]
    T = 2000
    res = least_squares(residual, x0, args=(I,T), x_scale=x_scale)
    ne_fit[i], V0_fit[i] = res.x
    x0 = res.x

plot = plt.plot

plt.figure()
plot(data['ne'], alt, label='Ground truth')
plot(ne_jacobsen, alt, label='Jacobsen')
plot(ne_fit, alt, label='Fit')
plt.xlabel('Density $[\mathrm{m}^{-3}]$')
plt.ylabel('Altitude $[\mathrm{km}]$')
plt.legend()

plt.figure()
plot(data['V0'], alt, label='Ground truth')
plot(V0_fit, alt, label='Fit')
plt.xlabel('Floating potential $[\mathrm{V}]$')
plt.ylabel('Altitude $[\mathrm{km}]$')
plt.legend()

plt.show()
