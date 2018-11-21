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

from langmuir import *
from pytest import approx

def test_thermal_current_maxwellian_normalized():
    sp = Species(n=1e11, T=1000)

    geo = Cylinder(3*sp.debye, 1)
    I = thermal_current(geo, sp, normalize=True)
    assert(I == approx(1))

    geo = Sphere(3*sp.debye)
    I = thermal_current(geo, sp, normalize=True)
    assert(I == approx(1))

def test_thermal_OML_current_maxwellian_normalized():
    sp = Species(n=1e11, T=1000)

    geo = Cylinder(3*sp.debye, 1)
    I = OML_current(geo, sp, eta=0, normalize=True)
    assert(I == approx(1))

    geo = Sphere(3*sp.debye)
    I = OML_current(geo, sp, eta=0, normalize=True)
    assert(I == approx(1))

def test_OML_current_maxwellian():
    I = OML_current(Cylinder(1e-3, 25e-3), Species(n=1e11, T=1000), 5)
    assert(I == approx(-1.0714855994312037e-06))
