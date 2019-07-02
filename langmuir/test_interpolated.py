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

# I generally try to test one point in each table, with unequal row/column
# indices and where nearby elements differ from the one I test.

@pytest.fixture
def electron():
    return Species(n=1e11, T=1000)

@pytest.fixture
def eta_fl():
    return np.array([2, 10, 17, 25, 32, 50, 75, 100], dtype=np.float)

def test_laframboise_sphere(electron):
    I = finite_radius_current(Sphere(electron.debye), electron, eta=-0.6,
                              normalize=True, table='laframboise')
    assert(I == approx(1.595))

def test_laframboise_cylinder(electron):
    I = finite_radius_current(Cylinder(3*electron.debye, 1), electron, eta=-2.0,
                              normalize=True, table='laframboise')
    assert(I == approx(1.928))

def test_darian_marholm_uncomplete_sphere(electron):

    def curr(kappa, alpha):
        sp = Species(n=1e11, T=1000, kappa=kappa, alpha=alpha)
        return finite_radius_current(Sphere(2*sp.debye), sp, eta=-1,
                                     normalize=True,
                                     table='darian-marholm uncomplete')

    # Testing one element for each alpha/kappa table
    assert(curr(float('inf'), 0  ) == approx(1.958))
    assert(curr(float('inf'), 0.2) == approx(2.069))
    assert(curr(6           , 0  ) == approx(1.982))
    assert(curr(6           , 0.2) == approx(2.369))
    assert(curr(4           , 0  ) == approx(1.990))
    assert(curr(4           , 0.2) == approx(2.846))

    # Testing that V=0 is the inaccurate one and that R=0 doesn't exist
    I = finite_radius_current(Sphere(0.2*electron.debye), electron, 0,
                              normalize=True, table='darian-marholm uncomplete')
    assert(I == approx(0.965))

def test_darian_marholm_sphere(electron):

    # Testing one from the uncomplete subset
    sp = Species(n=1e11, T=1000, kappa=6)
    I = finite_radius_current(Sphere(2*sp.debye), sp, eta=-1,
                              normalize=True, table='darian-marholm')
    assert(I == approx(1.982))

    # Testing that V=0 is the accurate one (not the one from uncomplete)
    I = finite_radius_current(Sphere(0.2*electron.debye), electron, eta=0,
                              normalize=True)
    assert(I == approx(1.0))

    # Testing that R=0 exists
    I = finite_radius_current(Sphere(0), electron, eta=-3,
                              normalize=True, table='darian-marholm')
    assert(I == approx(4.0))

def test_darian_marholm_uncomplete_cylinder(electron):

    def curr(kappa, alpha):
        sp = Species(n=1e11, T=1000, kappa=kappa, alpha=alpha)
        return finite_radius_current(Cylinder(3*sp.debye, 1), sp, eta=-1,
                                     normalize=True,
                                     table='darian-marholm uncomplete')

    assert(curr(float('inf'), 0  ) == approx(1.538))
    assert(curr(float('inf'), 0.2) == approx(1.914))
    assert(curr(6           , 0  ) == approx(1.509))
    assert(curr(6           , 0.2) == approx(2.165))
    assert(curr(4           , 0  ) == approx(1.510))
    assert(curr(4           , 0.2) == approx(2.547))

    # Testing that V=0 is the inaccurate one and that R=0 doesn't exist
    I = finite_radius_current(Cylinder(1.0*electron.debye, 1), electron, 0,
                              normalize=True, table='darian-marholm uncomplete')
    assert(I == approx(0.974))

def test_darian_marholm_cylinder(electron):

    # Testing one from the uncomplete subset
    sp = Species(n=1e11, T=1000, kappa=6)
    I = finite_radius_current(Cylinder(3*sp.debye, 1), sp, eta=-1,
                              normalize=True, table='darian-marholm')
    assert(I == approx(1.509))

    # Testing that V=0 is the accurate one (not the one from uncomplete)
    I = finite_radius_current(Cylinder(electron.debye, 1), electron, eta=0,
                              normalize=True)
    assert(I == approx(1.0))

    # Testing that R=0 exists
    I = finite_radius_current(Cylinder(0,1), electron, eta=-3,
                              normalize=True, table='darian-marholm')
    assert(I == approx(2.2417, 1e-3))

def test_laframboise_darian_marholm_sphere(electron):

    # Testing one from each of Darian-Marholm and Laframboise

    sp = Species(n=1e11, T=1000, kappa=4, alpha=0.2)
    I = finite_radius_current(Sphere(3*sp.debye), sp, eta=-5, normalize=True,
                              table='laframboise-darian-marholm')
    assert(I == approx(4.535))

    I = finite_radius_current(Sphere(3*electron.debye), electron, eta=-5,
                              normalize=True,
                              table='laframboise-darian-marholm')
    assert(I == approx(4.640))

def test_laframboise_darian_marholm_cylinder(electron):

    # Testing one from each of Darian-Marholm and Laframboise
    sp = Species(n=1e11, T=1000, kappa=4, alpha=0.2)
    I = finite_radius_current(Cylinder(3*sp.debye, 1), sp, eta=-5,
                              normalize=True,
                              table='laframboise-darian-marholm')
    assert(I == approx(3.465))

    I = finite_radius_current(Cylinder(3*electron.debye, 1), electron, eta=-5,
                              normalize=True,
                              table='laframboise-darian-marholm')
    assert(I == approx(2.701))

def test_finite_radius_vs_OML():

    sp = Species(n=1e11, T=1000, alpha=0.2, kappa=6)
    geo = Sphere(0.2*sp.debye)
    I_OML = OML_current(geo, sp, eta=-15)
    I_fr = finite_radius_current(geo, sp, eta=-15)
    assert(I_OML == approx(I_fr, 0.03))

def test_discard_spectral_indices(caplog):

    sp = Species(n=1e11, T=1000, alpha=0.2, kappa=6)
    geo = Sphere(0.2*sp.debye)

    I = finite_radius_current(geo, sp, eta=-15, table='laframboise')
    assert(len(caplog.records) == 1)

def test_discard_repelled_particles(caplog):

    electron = Species(n=1e11, T=1000)
    I = finite_radius_current(Sphere(electron.debye), electron, eta=[-1,1])
    assert(len(caplog.records) == 1)
    assert('negative' in caplog.text)
    assert('positive' not in caplog.text)

    ion = Species('proton', n=1e11, T=1000)
    I = finite_radius_current(Sphere(ion.debye), ion, eta=[-1,1])
    assert(len(caplog.records) == 2)
    assert('positive' in caplog.text)

def test_domain(caplog, electron):

    # Too high voltage
    I = finite_radius_current(Sphere(electron.debye), electron, eta=-30)
    assert(np.isnan(I))

    # Too large radius
    I = finite_radius_current(Sphere(200*electron.debye), electron, eta=-20)
    assert(np.isnan(I))

    # Too high alpha
    sp = Species(n=1e11, T=1000, alpha=0.3)
    I = finite_radius_current(Sphere(sp.debye), sp, eta=-20)
    assert(np.isnan(I))

    # Too low kappa
    sp = Species(n=1e11, T=1000, kappa=0.2)
    I = finite_radius_current(Sphere(sp.debye), sp, eta=-20)
    assert(np.isnan(I))

    assert(len(caplog.records) == 4)

def test_finite_length_normalize_OML(electron, eta_fl):

    geo = Cylinder(r=0.1*electron.debye, l=10*electron.debye)
    g = finite_length_current(geo, electron, eta=-eta_fl, normalize='OML')
    I = finite_length_current(geo, electron, eta=-eta_fl)
    I0 = OML_current(geo, electron, eta=-eta_fl)

    assert np.allclose(g, I/I0)

def test_finite_length_normalize_th(electron, eta_fl):

    geo = Cylinder(r=0.1*electron.debye, l=10*electron.debye)
    g = finite_length_current(geo, electron, eta=-eta_fl, normalize='th')
    I = finite_length_current(geo, electron, eta=-eta_fl)
    I0 = thermal_current(geo, electron)

    assert np.allclose(g, I/I0)

def test_finite_length_density_normalize_OML(electron, eta_fl):

    geo = Cylinder(r=0.1*electron.debye, l=10*electron.debye)
    geonorm = Cylinder(r=0.1*electron.debye, l=1)
    zeta = np.linspace(1,10,10)
    for eta in eta_fl:
        g = finite_length_current_density(geo, electron, eta=-eta, zeta=zeta, normalize='OML')
        i = finite_length_current_density(geo, electron, eta=-eta, zeta=zeta)
        i0 = OML_current(geonorm, electron, eta=-eta)
        gf = i/i0

        assert np.allclose(g, gf)

def test_finite_length_density_normalize_th(electron, eta_fl):

    geo = Cylinder(r=0.1*electron.debye, l=10*electron.debye)
    geonorm = Cylinder(r=0.1*electron.debye, l=1)
    zeta = np.linspace(1,10,10)
    for eta in eta_fl:
        g = finite_length_current_density(geo, electron, eta=-eta, zeta=zeta, normalize='th')
        i = finite_length_current_density(geo, electron, eta=-eta, zeta=zeta)
        i0 = thermal_current(geonorm, electron)
        gf = i/i0

        assert np.allclose(g, gf)

def test_finite_length_vs_OML(electron):

    # Agrees to within 5% of OML for infinite probes

    geo = Cylinder(r=electron.debye, l=10*electron.debye,
                  lguard=float('inf'), rguard=500*electron.debye)
    for eta in [2, 6, 10, 17, 25, 32, 50, 75, 100]:
        I = finite_length_current(geo, electron, eta=-eta, normalize='OML')
        assert np.allclose(I, 1, 0.05, 0), "eta={}".format(-eta)

def test_finite_length_density_vs_OML(electron):

    # Agrees to within 5% of OML for infinite probes

    lambd = 500
    geo_fl  = Cylinder(r=electron.debye, l=lambd*electron.debye)
    geo_OML = Cylinder(r=electron.debye, l=1)
    for eta in [2, 6, 10, 17, 25, 32, 50, 75, 100]:
        I_OML = OML_current(geo_OML, electron, eta=-eta)
        I_fl = finite_length_current_density(geo_fl, electron, eta=-eta, zeta=lambd/2)
        assert I_OML == approx(I_fl, 0.05), "eta={}".format(-eta)

def test_finite_length_samples():

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

def test_finite_length_density_samples():

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

def test_finite_length_interpolated_samples():

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

def test_finite_length_density_interpolated_samples():

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
