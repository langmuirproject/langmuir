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
    assert(s2.dist == s.dist)

def test_plane_repr():
    # A good representation of a (small) object is one you can use to recreate
    # the object
    g  = Plane(3)
    g2 = eval(str(g))
    assert(g2.A == approx(g.A))

def test_cylinder_repr():
    # A good representation of a (small) object is one you can use to recreate
    # the object
    g  = Cylinder(3,2)
    g2 = eval(str(g))
    assert(g2.r == approx(g.r))
    assert(g2.l == approx(g.l))

def test_sphere_repr():
    # A good representation of a (small) object is one you can use to recreate
    # the object
    g  = Sphere(3)
    g2 = eval(str(g))
    assert(g2.r == approx(g.r))

def test_normalization_current():
    sph = Sphere(1e-3)
    cyl = Cylinder(1e-3,10e-3)
    sp  = Species(n=1e11, T=1000)
    assert(normalization_current(sph, sp) == approx(-9.888431090271652e-09))
    assert(normalization_current(cyl, sp) == approx(-4.944215545135826e-08))

    with pytest.raises(ValueError):
        normalization_current(None, sp)

def test_normalization_current_multiple_species():
    cyl = Cylinder(1e-3,10e-3)
    plasma = []
    plasma.append(Species('electron', n=1e11, T=1000))
    plasma.append(Species('proton',   n=1e11, T=1000))
    I0_n = normalization_current(cyl, plasma[0])
    I0_p = normalization_current(cyl, plasma[1])
    assert(normalization_current(cyl, plasma) == approx(I0_n+I0_p))

def test_OML_current_not_normalized_by_default():
    geometry = Cylinder(0.255e-3, 25e-3)
    species  = Species(n=10e10, eV=0.26)
    I1 = OML_current(geometry, species, 5)
    I2 = OML_current(geometry, species, 5, normalize=False)
    assert(I1 == approx(I2))


def test_Debye_length():
    sp_m = Species(n=1e11, T=1000)
    sp_k = Species('kappa', n=1e11, T=1000, kappa=4)
    sp_a = Species('cairns', n=1e11, T=1000, alpha=0.2)
    sp_ka = Species('kappa-cairns', n=1e11, T=1000, kappa=4, alpha=0.2)
    assert(sp_m.debye == approx(6.900896931e-03))
    assert(sp_k.debye == approx(0.0058323223992018209))
    assert(sp_a.debye == approx(0.010911276093510305))
    assert(sp_ka.debye == approx(0.011952704095164884))
