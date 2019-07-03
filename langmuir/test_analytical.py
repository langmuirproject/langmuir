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
import pytest

def test_thermal_current_maxwellian_normalized():
    sp = Species(n=1e11, T=1000)

    geo = Cylinder(3*sp.debye, 1)
    I = thermal_current(geo, sp, normalization='thmax')
    assert(I == approx(1))

    geo = Sphere(3*sp.debye)
    I = thermal_current(geo, sp, normalization='thmax')
    assert(I == approx(1))

def test_thermal_current_cairns_normalized():
    sp = Species(n=1e11, T=1000, alpha=0.2)
    sp_m = Species(n=1e11, T=1000, alpha=0.0)

    geo = Cylinder(3 * sp.debye, 1)
    I = thermal_current(geo, sp, normalization='thmax')
    I_m = thermal_current(geo, sp_m, normalization='thmax')
    assert(I == approx(0.4))
    assert(I_m == approx(1))

    geo = Sphere(3 * sp.debye)
    I = thermal_current(geo, sp, normalization='thmax')
    I_m = thermal_current(geo, sp_m, normalization='thmax')
    assert(I == approx(1.45))
    assert(I_m == approx(1))

def test_thermal_current_kappa_normalized():
    sp = Species(n=1e11, T=1000, kappa=4)

    geo = Cylinder(3 * sp.debye, 1)
    I = thermal_current(geo, sp, normalization='thmax')
    assert(I == approx(1.8446619684315548))

    geo = Sphere(3 * sp.debye)
    I = thermal_current(geo, sp, normalization='thmax')
    assert(I == approx(0.95153286194814457))

def test_thermal_current_kappa_cairns_normalized():
    sp = Species(n=1e11, T=1000, alpha=0.2, kappa=4)
    sp_k = Species(n=1e11, T=1000, alpha=0, kappa=4)

    geo = Cylinder(3 * sp.debye, 1)
    I = thermal_current(geo, sp, normalization='thmax')
    I_k = thermal_current(geo, sp_k, normalization='thmax')
    assert(I == approx(0.43920523057894167))
    assert(I_k == approx(1.8446619684315548))

    geo = Sphere(3 * sp.debye)
    I = thermal_current(geo, sp, normalization='thmax')
    I_k = thermal_current(geo, sp_k, normalization='thmax')
    assert(I == approx(2.5374209651950519))
    assert(I_k == approx(0.95153286194814457))

def test_thermal_OML_current_maxwellian_normalized():
    sp = Species(n=1e11, T=1000)

    geo = Cylinder(3*sp.debye, 1)
    I = OML_current(geo, sp, eta=0, normalization='thmax')
    assert(I == approx(1))

    geo = Sphere(3*sp.debye)
    I = OML_current(geo, sp, eta=0, normalization='thmax')
    assert(I == approx(1))

def test_OML_current_maxwellian():
    I = OML_current(Cylinder(1e-3, 25e-3), Species(n=1e11, T=1000), 5)
    assert(I == approx(-1.0714855994312037e-06))

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

def test_OML_current_not_normalized():
    geometry = Cylinder(0.255e-3, 25e-3)
    species  = Species(n=10e10, eV=0.26)
    I1 = OML_current(geometry, species, 5)
    assert(I1 == approx(-2.777872191401695e-07))

def test_OML_current_cairns():
    sp = Species(n=1e11, T=1000, alpha=0.2)
    sp_m = Species(n=1e11, T=1000, alpha=0.0)

    # Cylindrical geometry
    geo = Cylinder(0.255e-3, 25e-3)
    # Attracted particles - Cairns VDF
    I = OML_current(geo, sp, eta=-5)
    assert(I == approx(-9.3492e-8))
    #  Attracted particles - Maxwellian VDF
    I = OML_current(geo, sp_m, eta=-5)
    assert(I == approx(-8.68504e-8))
    # Repelled particles - Cairns VDF
    I = OML_current(geo, sp, eta=5)
    assert(I == approx(-2.21933e-9))
    # Repelled particles - Maxwellian VDF
    I = OML_current(geo, sp_m, eta=5)
    assert(I == approx(-2.12376e-10))

    # Spherical geometry
    geo = Sphere(0.255e-3)
    # Attracted particles - Cairns VDF
    I = OML_current(geo, sp, eta=-5)
    assert(I == approx(-3.02208e-9))
    # Attracted particles - Maxwellian VDF
    I = OML_current(geo, sp_m, eta=-5)
    assert(I == approx(-3.85797e-9))
    # Repelled particles - Cairns VDF
    I = OML_current(geo, sp, eta=5)
    assert(I == approx(-4.52743e-11))
    # Repelled particles - Maxwellian VDF
    I = OML_current(geo, sp_m, eta=5)
    assert(I == approx(-4.33247e-12))

def test_OML_current_kappa_cairns():
    sp = Species(n=1e11, T=1000, kappa=4.0, alpha=0.2)
    sp_k = Species(n=1e11, T=1000, kappa=4.0, alpha=0.0)

    # Cylindrical geometry
    geo = Cylinder(0.255e-3, 25e-3)
    # Attracted particles - Kappa-Cairns VDF
    I = OML_current(geo, sp, eta=-5)
    assert(I == approx(-1.19482e-07))
    # Attracted particles - Kappa VDF
    I = OML_current(geo, sp_k, eta=-5)
    assert(I == approx(-8.66487e-08))
    # Repelled particles - Kappa-Cairns VDF
    I = OML_current(geo, sp, eta=5)
    assert(I == approx(-3.9989e-8))
    # Repelled particles - Kappa VDF
    I = OML_current(geo, sp_k, eta=5)
    assert(I == approx(-1.1108e-9))

    # Spherical geometry
    geo = Sphere(0.255e-3)
    # Attracted particles - Kappa-Cairns VDF
    I = OML_current(geo, sp, eta=-5)
    assert(I == approx(-3.2631e-9))
    # Attracted particles - Kappa VDF
    I = OML_current(geo, sp_k, eta=-5)
    assert(I == approx(-4.28282e-9))
    # Repelled particles - Kappa-Cairns VDF
    I = OML_current(geo, sp, eta=5)
    assert(I == approx(-8.15775e-10))
    # Repelled particles - Kappa VDF
    I = OML_current(geo, sp_k, eta=5)
    assert(I == approx(-2.26604e-11))
