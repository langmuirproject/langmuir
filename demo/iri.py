from scipy.constants import value as constants
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import pyglow
from tqdm import tqdm

lat = 45.
lon = 0.
alts = np.linspace(60., 1300., 1001)
dn = datetime(2008, 3, 1, 12, 00)

ns = {}
Ts = {}

for alt in tqdm(alts):
    pt = pyglow.Point(dn, lat, lon, alt)
    pt.run_iri()
    pt.run_msis()

    n = pt.ni
    n['e'] = pt.ne
    for key in n:
        value = n[key]*1e6
        if value < 1e8: value = float('nan')
        if not key in ns:
            ns[key] = []
        ns[key].append(value)

    T = {'e': pt.Te, 'i': pt.Ti}
    for key in T:
        value = T[key]
        if value < 1: value = float('nan')
        if not key in Ts:
            Ts[key] = []
        Ts[key].append(value)

np.savez_compressed('iri.npz',
                    alt=alts,
                    Te=np.array(Ts['e']),
                    Ti=np.array(Ts['i']),
                    ne=np.array(ns['e']),
                    nO=np.array(ns['O+']),
                    nO2=np.array(ns['O2+']),
                    nNO=np.array(ns['NO+']),
                    nH=np.array(ns['H+'])
                   )


plt.figure()
plot = plt.loglog
plot(ns['e'], alts, label='e')
plot(ns['O+'], alts, label='O+')
plot(ns['O2+'], alts, label='O2+')
plot(ns['NO+'], alts, label='NO+')
plot(ns['H+'], alts, label='H+')
plt.legend()

plt.figure()
plot = plt.semilogy
plot(Ts['e'], alts, label='e')
plot(Ts['i'], alts, label='i')
plt.legend()

plt.show()
