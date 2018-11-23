"""
Copyright 2018
    Sigvald Marholm <marholm@marebakken.com>
    Diako Darian <diakod@math.uio.no>

This file is part of langmuir.

langmuir is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

langmuir is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with langmuir.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import division
from langmuir import *
from scipy.constants import value as constants
from pytest import approx
import pytest
import numpy as np

def test_species_defaults():
    T = 1000
    n = 1e11
    e = constants('elementary charge')
    m = constants('electron mass')
    k = constants('Boltzmann constant')
    s = Species(n=n, T=T)
    assert(s.m == approx(m))
    assert(s.q == approx(-e))
    assert(s.n == approx(n))
    assert(s.T == approx(T))
    assert(s.vth == approx(np.sqrt(k*1000/m)))
    assert(s.alpha == 0)
    assert(s.kappa == float('inf'))

def test_species_electron():
    s = Species('electron', n=1e11, T=1000)
    assert(s.m == approx(constants('electron mass')))
    assert(s.q == approx(-constants('elementary charge')))

def test_species_proton():
    s = Species('proton', n=1e11, T=1000)
    assert(s.m == approx(constants('proton mass')))
    assert(s.q == approx(constants('elementary charge')))

def test_species_positron():
    s = Species('positron', n=1e11, T=1000)
    assert(s.m == approx(constants('electron mass')))
    assert(s.q == approx(constants('elementary charge')))

def test_species_maxwellian():
    s = Species(n=1e11, T=1000)
    assert(s.kappa == float('inf'))
    assert(s.alpha == approx(0))

def test_species_kappa():
    s = Species(n=1e11, T=1000, kappa=4)
    assert(s.kappa == approx(4))
    assert(s.alpha == approx(0))

def test_species_cairns():
    s = Species(n=1e11, T=1000, alpha=0.2)
    assert(s.kappa == float('inf'))
    assert(s.alpha == approx(0.2))

def test_species_kappa_cairns():
    s = Species(n=1e11, T=1000, kappa=4, alpha=0.2)
    assert(s.kappa == approx(4))
    assert(s.alpha == approx(0.2))

def test_species_q():
    s = Species(n=1e11, T=1000, q=123)
    assert(s.q == approx(123))

def test_species_Z():
    s = Species(n=1e11, T=1000, Z=2)
    assert(s.q == approx(2*constants('elementary charge')))

def test_species_m():
    s = Species(n=1e11, T=1000, m=234)
    assert(s.m == approx(234))

def test_species_amu():
    s = Species(n=1e11, T=1000, amu=16)
    assert(s.m == approx(16*constants('atomic mass constant')))

def test_species_vth():
    s = Species(n=1e11, vth=1e6)
    assert(s.T == approx(1e6**2*s.m/constants('Boltzmann constant')))
    assert(s.vth == approx(1e6))

def test_species_eV():
    s = Species(n=1e11, eV=0.2)
    e = constants('elementary charge')
    k = constants('Boltzmann constant')
    assert(s.T == approx(0.2*e/k))

def test_species_repr():
    # A good representation of a (small) object is one you can use to recreate
    # the object
    s  = Species('kappa-cairns', q=1, m=3, n=1000, T=1000, kappa=6, alpha=0.2)
    s2 = eval(str(s))
    assert(s2.q == approx(s.q))
    assert(s2.m == approx(s.m))
    assert(s2.n == approx(s.n))
    assert(s2.T == approx(s.T))
    assert(s2.kappa == approx(s.kappa))
    assert(s2.alpha == approx(s.alpha))

def test_species_convenience_values():
    s = Species(n=1e11, T=1000)
    omega_p = 17839863.646512896
    eV = 0.08617330337217212
    vth = 123111.06021431087
    assert(s.omega_p == approx(omega_p))
    assert(s.freq_p == approx(omega_p/(2*np.pi)))
    assert(s.period_p == approx(2*np.pi/omega_p))
    assert(s.eV == approx(eV))
    assert(s.vth == approx(vth))

def test_species_convenience_methods():
    s = Species(n=1e11, T=1000)
    B = 50e-6
    omega_c = 8794100.11801062
    assert(s.omega_c(B) == approx(omega_c))
    assert(s.freq_c(B) == approx(omega_c/(2*np.pi)))
    assert(s.period_c(B) == approx(2*np.pi/omega_c))
    assert(s.larmor(B) == approx(s.vth/omega_c))

def test_species_debye_length():
    sp_m = Species(n=1e11, T=1000)
    sp_k = Species('kappa', n=1e11, T=1000, kappa=4)
    sp_a = Species('cairns', n=1e11, T=1000, alpha=0.2)
    sp_ka = Species('kappa-cairns', n=1e11, T=1000, kappa=4, alpha=0.2)
    assert(sp_m.debye == approx(6.900896931e-03))
    assert(sp_k.debye == approx(0.0058323223992018209))
    assert(sp_a.debye == approx(0.010911276093510305))
    assert(sp_ka.debye == approx(0.011952704095164884))

def test_debye():
    plasma = []
    plasma.append(Species('electron', n=1e11, T=1000))
    plasma.append(Species('proton',   n=1e11, T=1000))
    assert(debye(plasma[0])==plasma[0].debye)
    assert(debye(plasma) == approx(0.004879671013271479))
