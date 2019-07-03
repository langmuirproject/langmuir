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

@pytest.fixture
def electron():
    return Species(n=1e11, T=1000)

@pytest.fixture
def eta_fl():
    return np.array([2, 10, 17, 25, 32, 50, 75, 100], dtype=np.float)

def test_current_normalize_OML(electron, eta_fl):

    geo = Cylinder(r=0.1*electron.debye, l=10*electron.debye)
    g = finite_length_current(geo, electron, eta=-eta_fl, normalization='OML')
    I = finite_length_current(geo, electron, eta=-eta_fl)
    I0 = OML_current(geo, electron, eta=-eta_fl)

    assert np.allclose(g, I/I0)

def test_current_normalize_th(electron, eta_fl):

    geo = Cylinder(r=0.1*electron.debye, l=10*electron.debye)
    g = finite_length_current(geo, electron, eta=-eta_fl, normalization='th')
    I = finite_length_current(geo, electron, eta=-eta_fl)
    I0 = thermal_current(geo, electron)

    assert np.allclose(g, I/I0)

def test_current_density_normalize_OML(electron, eta_fl):

    geo = Cylinder(r=0.1*electron.debye, l=10*electron.debye)
    geonorm = Cylinder(r=0.1*electron.debye, l=1)
    zeta = np.linspace(1,10,10)
    for eta in eta_fl:
        g = finite_length_current_density(geo, electron, eta=-eta, zeta=zeta, normalization='OML')
        i = finite_length_current_density(geo, electron, eta=-eta, zeta=zeta)
        i0 = OML_current(geonorm, electron, eta=-eta)
        gf = i/i0

        assert np.allclose(g, gf)

def test_current_density_normalize_th(electron, eta_fl):

    geo = Cylinder(r=0.1*electron.debye, l=10*electron.debye)
    geonorm = Cylinder(r=0.1*electron.debye, l=1)
    zeta = np.linspace(1,10,10)
    for eta in eta_fl:
        g = finite_length_current_density(geo, electron, eta=-eta, zeta=zeta, normalization='th')
        i = finite_length_current_density(geo, electron, eta=-eta, zeta=zeta)
        i0 = thermal_current(geonorm, electron)
        gf = i/i0

        assert np.allclose(g, gf)

def test_current_vs_OML(electron):

    # Agrees to within 5% of OML for infinite probes

    geo = Cylinder(r=electron.debye, l=10*electron.debye,
                  # lguard=float('inf'), rguard=500*electron.debye)
                  lguard=float('inf'), rguard=float('inf'))
    for eta in [2, 6, 10, 17, 25, 32, 50, 75, 100]:
        I = finite_length_current(geo, electron, eta=-eta, normalization='OML')
        assert np.allclose(I, 1, 0.05, 0), "eta={}".format(-eta)

def test_current_density_vs_OML(electron):

    # Agrees to within 5% of OML for infinite probes

    geo = Cylinder(r=electron.debye, l=1,
                   lguard=float('inf'), rguard=float('inf'))
    for eta in [2, 6, 10, 17, 25, 32, 50, 75, 100]:
        I_OML = OML_current(geo, electron, eta=-eta)
        I_fl = finite_length_current_density(geo, electron, eta=-eta, z=0.5)
        assert I_OML == approx(I_fl, 0.05), "eta={}".format(-eta)

def test_current_guard(electron):

    r = electron.debye
    l = 10*electron.debye
    long = 1000*electron.debye
    eta = 17

    Ino   = finite_length_current(Cylinder(r, l), electron, eta=-eta)

    Iinf  = finite_length_current(Cylinder(r, l, lguard=np.inf), electron, eta=-eta)
    Ilong = finite_length_current(Cylinder(r, l, lguard=long), electron, eta=-eta)
    assert Iinf == approx(Ilong)
    assert not Iinf == approx(Ino)

    Iinf  = finite_length_current(Cylinder(r, l, rguard=np.inf), electron, eta=-eta)
    Ilong = finite_length_current(Cylinder(r, l, rguard=long), electron, eta=-eta)
    assert Iinf == approx(Ilong)
    assert not Iinf == approx(Ino)

    Iinf  = finite_length_current(Cylinder(r, l, lguard=l, rguard=np.inf), electron, eta=-eta)
    Ilong = finite_length_current(Cylinder(r, l, lguard=l, rguard=long), electron, eta=-eta)
    assert Iinf == approx(Ilong)

def test_current_density_guard(electron):

    r = electron.debye
    l = 10*electron.debye
    long = 1000*electron.debye
    eta = 17

    Ino   = finite_length_current_density(Cylinder(r, l), electron, eta=-eta)

    Iinf  = finite_length_current_density(Cylinder(r, l, lguard=np.inf), electron, eta=-eta)
    Ilong = finite_length_current_density(Cylinder(r, l, lguard=long), electron, eta=-eta)
    assert Iinf == approx(Ilong)
    assert not Iinf == approx(Ino)

    Iinf  = finite_length_current_density(Cylinder(r, l, rguard=np.inf), electron, eta=-eta)
    Ilong = finite_length_current_density(Cylinder(r, l, rguard=long), electron, eta=-eta)
    assert Iinf == approx(Ilong)
    assert not Iinf == approx(Ino)

    Iinf  = finite_length_current_density(Cylinder(r, l, lguard=l, rguard=np.inf), electron, eta=-eta)
    Ilong = finite_length_current_density(Cylinder(r, l, lguard=l, rguard=long), electron, eta=-eta)
    assert Iinf == approx(Ilong)

def test_current_samples():

    # Testing that a few random sample points don't change. All points are in
    # the original data, so the test should be independent of interpolation

    elec = Species(n=35e10, eV=0.08)

    r = 1e-3
    ls = np.array([5, 30, 80, 400, 2000], dtype=np.float)*1e-3
    etas = np.array([100, 2, 32, 17, 25], dtype=np.float)
    Ifs = np.array([-7.311843622376143e-06,
                    -1.2270514750624203e-06,
                    -1.332414356159239e-05,
                    -3.5970798600771154e-05,
                    -0.00020182061197445252])

    for l, eta, If in zip(ls, etas, Ifs):
        geo = Cylinder(r=r, l=l)
        I = finite_length_current(geo, elec, eta=-eta)
        assert I == approx(If)

def test_current_density_samples():

    # Testing that a few random sample points don't change. All points are in
    # the original data, so the test should be independent of interpolation

    elec = Species(n=35e10, eV=0.08)

    r = 1e-3
    ls = np.array([5, 30, 80, 400, 2000], dtype=np.float)*1e-3
    etas = np.array([100, 2, 32, 17, 25], dtype=np.float)
    zs = np.array([0, .22, .1, .75, 0.91, 1])
    Ifs = np.array([[-0.00126263, -0.00147027, -0.00137461, -0.00148828, -0.00136481, -0.00126263],
                    [-4.02337566e-05, -4.10898286e-05, -4.26822956e-05, -4.07257345e-05, -4.27430040e-05, -4.02337566e-05],
                    [-0.00021501, -0.00013982, -0.00021829, -0.00013476, -0.00023155, -0.00021501],
                    [-1.41492521e-04, -8.50936721e-05, -8.51884293e-05, -8.50936707e-05, -8.53236120e-05, -1.41492521e-04],
                    [-1.76726210e-04, -9.93476796e-05, -9.93476796e-05, -9.93476796e-05, -9.93476796e-05, -1.76726210e-04]])

    for l, eta, If in zip(ls, etas, Ifs):
        geo = Cylinder(r=r, l=l)
        I = finite_length_current_density(geo, elec, eta=-eta, z=zs*l)
        assert np.allclose(I, If)

def test_current_interpolated_samples():

    # Testing that a few random sample points don't change. These points depen
    # upon the interpolation scheme. 

    elec = Species(n=35e10, eV=0.08)

    r = 1e-3
    ls = np.array([7, 70, 700], dtype=np.float)*1e-3
    etas = np.array([65, 1, 90], dtype=np.float)
    Ifs = np.array([-6.2794005496528275e-06,
                    -1.9760221828553347e-06,
                    -0.00014456804182897328])

    for l, eta, If in zip(ls, etas, Ifs):
        geo = Cylinder(r=r, l=l)
        I = finite_length_current(geo, elec, eta=-eta)
        assert I == approx(If)

def test_current_density_interpolated_samples():

    # Testing that a few random sample points don't change. These points depen
    # upon the interpolation scheme. 

    elec = Species(n=35e10, eV=0.08)

    r = 1e-3
    ls = np.array([7, 70, 700], dtype=np.float)*1e-3
    etas = np.array([65, 1, 90], dtype=np.float)
    zs = np.array([0, .4, .29, 1])
    Ifs = np.array([[-0.00045615, -0.00106614, -0.00099727, -0.00045615],
                    [-2.92817683e-05, -2.75155255e-05, -2.76705163e-05, -2.92817683e-05],
                    [-0.00036732, -0.00019135, -0.00019135, -0.00036732]])

    for l, eta, If in zip(ls, etas, Ifs):
        geo = Cylinder(r=r, l=l)
        I = finite_length_current_density(geo, elec, eta=-eta, z=zs*l)
        assert np.allclose(I, If)

def test_current_non_maxwellian(caplog):

    elec = Species(n=35e10, eV=0.08)
    geo = Cylinder(0.1*elec.debye, 10*elec.debye)
    I = finite_length_current(geo, Species(n=35e10, T=1000, alpha=0.3), 3.0)
    I = finite_length_current(geo, Species(n=35e10, T=1000, kappa=5), 3.0)

    assert(len(caplog.records) == 2)

def test_current_density_non_maxwellian(caplog):

    elec = Species(n=35e10, eV=0.08)
    geo = Cylinder(0.1*elec.debye, 10*elec.debye)
    I = finite_length_current_density(geo, Species(n=35e10, T=1000, alpha=0.3), 3.0)
    I = finite_length_current_density(geo, Species(n=35e10, T=1000, kappa=5), 3.0)

    assert(len(caplog.records) == 2)

def test_current_density_matrix(electron):

    geo = Cylinder(0.1*electron.debye, 10*electron.debye)
    eta = np.array([-1, 0, 1, 10])
    zeta = np.linspace(0, 10, 8)
    I_matrix = finite_length_current_density(geo, electron, eta=-eta, zeta=zeta)
    for i, e in enumerate(eta):
        I_row = finite_length_current_density(geo, electron, eta=-e, zeta=zeta)
        assert np.allclose(I_row, I_matrix[i])

def test_current_repelled(electron):

    geo = Cylinder(0.1*electron.debye, electron.debye)
    V = -np.array([0.1, 1.2, 56])
    I = finite_length_current(geo, electron, V)
    I_OML = OML_current(geo, electron, V)
    assert np.allclose(I, I_OML)

def test_current_density_repelled(electron):

    geo = Cylinder(0.1*electron.debye, 1)
    V = -np.array([0.1, 1.2, 56])
    I = finite_length_current_density(geo, electron, V)
    I_OML = OML_current(geo, electron, V)
    assert np.allclose(I, I_OML)
