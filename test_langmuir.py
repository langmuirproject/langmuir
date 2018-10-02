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

"""
A NOTE ON FLOATING POINT COMPARISON IN THESE TESTS:

It is often said that you should _never_ compare two floating point numbers for
equality using the == operator. I think this is a good rule, but you _can_
violate it if you understand what you're doing. Floating point errors occur in
computations, and in conversion between bases. Consider this example::

    >>> a = 0.1
    >>> b = 0.1
    >>> c = b
    >>> "{:.50f}".format(a)
    '0.10000000000000000555111512312578270211815834045410'

Here you have a conversion from base 10 to base 2 with as many bits as dictated by the programming language. Since 0.1 cannot be represented exactly in base 2, you get round-off errors. However, you should have the exact same error for b. And bits do not magically flip when moving the contents of one variable to another. Hence, comparing the three variables a, b, and c with == should be safe.
"""

from langmuir import *
from scipy.constants import value as constants
from pytest import approx
import pytest
import numpy as np

def test_lafr_attr_current_geometry():
    with pytest.raises(ValueError) as e_info:
        f = lafr_attr_current('Pentagon')

def test_lafr_norm_current_geometry():
    with pytest.raises(ValueError) as e_info:
        f = lafr_norm_current('Pentagon', 1, 1e11, 1e3)

def test_species_defaults():
    T = 1000
    n = 1e11
    e = constants('elementary charge')
    m = constants('electron mass')
    k = constants('Boltzmann constant')
    s = Species(n=n, T=T)
    assert(s.dist == 'maxwellian')
    assert(s.m == approx(m))
    assert(s.q == approx(-e))
    assert(s.n == approx(n))
    assert(s.T == approx(T))
    assert(s.vth == approx(np.sqrt(k*1000/m)))

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
    s = Species('maxwellian', n=1e11, T=1000)
    assert(s.dist == 'maxwellian')
    assert(s.kappa == float('inf'))
    assert(s.alpha == approx(0))

def test_species_kappa():
    s = Species('kappa', n=1e11, T=1000, kappa=4)
    assert(s.dist == 'kappa')
    assert(s.kappa == approx(4))
    assert(s.alpha == approx(0))

def test_species_cairns():
    s = Species('cairns', n=1e11, T=1000, alpha=0.2)
    assert(s.dist == 'cairns')
    assert(s.kappa == float('inf'))
    assert(s.alpha == approx(0.2))

def test_species_kappa_cairns():
    s = Species('kappa-cairns', n=1e11, T=1000, kappa=4, alpha=0.2)
    assert(s.dist == 'kappa-cairns')
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
