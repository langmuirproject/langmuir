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
