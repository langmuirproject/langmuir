#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import langmuir as l
from tqdm import tqdm
from itertools import count
from localreg import RBFnet, plot_corr
from localreg.metrics import rms_rel_error

geo = l.Cylinder(r=0.255e-3, l=25e-3, lguard=float('inf'))
model = l.finite_length_current
Vs = np.array([2, 3, 4, 5]) # bias voltages

"""
PART 1: GENERATE SYNTHETIC DATA USING LANGMUIR
"""

def rand_uniform(N, range):
    """Generate N uniformly distributed numbers in range"""
    return range[0]+(range[1]-range[0])*np.random.rand(N)

def rand_log(N, range):
    """Generate N logarithmically distributed numbers in range"""
    x = rand_uniform(N, np.log(range))
    return np.exp(x)

N = 1000
ns  = rand_log(N, [1e11, 12e11])  # densities
Ts = rand_uniform(N, [800, 2500]) # temperatures
V0s = rand_uniform(N, [-3,  0])   # floating potentials

# Generate probe currents corresponding to plasma parameters
Is = np.zeros((N,len(Vs)))
for i, n, T, V0 in zip(count(), ns, Ts, tqdm(V0s)):
    Is[i] = model(geo, l.Electron(n=n, T=T), V=V0+Vs)

"""
PART 2: TRAIN AND TEST THE REGRESSION NETWORK
"""

# Use M first data points for training, the rest for testing.
M = int(0.8*N)

# Train by minimizing relative error in density
net = RBFnet()
net.train(Is[:M], ns[:M], num=20, relative=True, measure=rms_rel_error)

# Plot and print error metrics on test data
pred = net.predict(Is[M:])

fig, ax = plt.subplots()
plot_corr(ax, ns[M:], pred, log=True)
print("RMS of relative density error: {:.1f}%".format(100*rms_rel_error(ns[M:], pred)))

"""
PART 3: PREDICT PLASMA PARAMETERS FROM ACTUAL DATA
"""

data = l.generate_synthetic_data(geo, Vs, model=model)
pred = net.predict(data['I'])

plt.figure()
plt.plot(data['ne'], data['alt'], label='Ground truth')
plt.plot(pred, data['alt'], label='Predicted')
plt.xlabel('Density $[\mathrm{m}^{-3}]$')
plt.ylabel('Altitude $[\mathrm{km}]$')
plt.legend()

plt.show()
