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
from pytest import approx

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
