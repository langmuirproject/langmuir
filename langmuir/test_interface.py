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
from scipy.constants import value as constants
try:
    from inspect import getfullargspec as getargspec # Python 3
except ImportError:
    from inspect import getargspec # Python 2

current_models = [
    thermal_current,
    normalization_current,
    OML_current,
    finite_radius_current,
    finite_length_current,
    finite_length_current_density
]

@pytest.fixture
def electron():
    return Species(n=1e11, T=1000)

@pytest.fixture
def proton():
    return Species('proton', n=1e11, T=1000)

@pytest.mark.parametrize("current", current_models)
def test_multiple_species(current, electron, proton):

    # Add a voltage for those models having a V argument.
    # This must be small to not make repelled current negligible.
    args = getargspec(current).args
    kwargs = {}
    if 'V' in args: kwargs['V'] = 0.1

    # geo = Sphere(electron.debye)
    geo = Cylinder(electron.debye, 10*electron.debye)
    Ie = current(geo, electron, **kwargs)
    Ii = current(geo, proton, **kwargs)
    I  = current(geo, [electron, proton], **kwargs)

    assert Ie+Ii == approx(I)

@pytest.mark.parametrize("current", current_models)
def test_normalization_thmax(current, electron):

    if current != normalization_current:

        args = getargspec(current).args
        kwargs = {}
        if 'V' in args: kwargs['V'] = 0.1

        geo = Cylinder(0.1*electron.debye, 1)
        In = current(geo, electron, normalization='thmax', **kwargs)
        I  = current(geo, electron, **kwargs)
        I0 = normalization_current(geo, electron)

        assert I0*In == approx(I)

@pytest.mark.parametrize("current", current_models)
def test_normalization_th(current, electron):

    if current not in [normalization_current, thermal_current]:

        args = getargspec(current).args
        kwargs = {}
        if 'V' in args: kwargs['V'] = 0.1

        geo = Cylinder(0.1*electron.debye, 1)
        In = current(geo, electron, normalization='th', **kwargs)
        I  = current(geo, electron, **kwargs)
        I0 = thermal_current(geo, electron)

        assert I0*In == approx(I)

@pytest.mark.parametrize("current", current_models)
def test_normalization_oml(current, electron):

    if current not in [normalization_current, thermal_current, OML_current]:

        args = getargspec(current).args
        kwargs = {}
        if 'V' in args: kwargs['V'] = 0.1

        geo = Cylinder(0.1*electron.debye, 1)
        In = current(geo, electron, normalization='oml', **kwargs)
        I  = current(geo, electron, **kwargs)
        I0 = OML_current(geo, electron, 0.1)

        assert I0*In == approx(I)

@pytest.mark.parametrize("current", current_models)
def test_eta(current, electron):

    args = getargspec(current).args
    kwargs = {}
    if 'V' in args or 'eta' in args:

        # Should have both or none
        assert 'V' in args
        assert 'eta' in args

        e = constants('elementary charge')
        k = constants('Boltzmann constant')
        T = 1000
        eta = -10

        geometry = Cylinder(r=electron.debye, l=10*electron.debye)
        I_eta = current(geometry, electron, eta=eta)
        I_V   = current(geometry, electron, V=-eta*k*T/e)

        assert I_eta == approx(I_V)

@pytest.mark.parametrize("current", current_models)
def test_multiple_species_eta(current, electron, proton, caplog):

    args = getargspec(current).args
    if 'eta' in args:

        I = current(Sphere(electron.debye), [electron, proton], eta=1)
        assert(caplog.records[0].levelname == 'ERROR')

@pytest.mark.parametrize("current", current_models)
def test_multiple_species_normalize(current, electron, proton, caplog):

    args = getargspec(current).args
    kwargs = {}
    geometry = Cylinder(r=electron.debye, l=10*electron.debye)
    if 'V' in args: kwargs['V'] = 0.1
    if 'normalization' in args:

        I = current(geometry, [electron, proton],
                    normalization='thmax', **kwargs)
        assert(caplog.records[0].levelname == 'ERROR')

@pytest.mark.parametrize("current", current_models)
def test_multiple_species_zeta(current, electron, proton, caplog):

    args = getargspec(current).args
    kwargs = {}
    geometry = Cylinder(r=electron.debye, l=10*electron.debye)
    if 'V' in args: kwargs['V'] = 0.1
    if 'zeta' in args:

        I = current(geometry, [electron, proton], zeta=39, **kwargs)
        assert(caplog.records[0].levelname == 'ERROR')

@pytest.mark.parametrize("current", current_models)
def test_input_output_format(current, electron):

    geo = Cylinder(r=0.1*electron.debye, l=10*electron.debye)
    # geo = Sphere(electron.debye)

    args = getargspec(current).args
    if 'V' in args:

        I = current(geo, electron, 0.1)
        assert(isinstance(I, float))

        I = current(geo, electron, 1)
        assert(isinstance(I, float))

        I = current(geo, electron, [1, 2])
        assert(isinstance(I, np.ndarray))
        assert(I.dtype==np.float)
        assert(I.shape==(2,))

        I = current(geo, electron, np.array([1, 2]))
        assert(isinstance(I, np.ndarray))
        assert(I.dtype==np.float)
        assert(I.shape==(2,))

@pytest.mark.parametrize("current", current_models)
def test_geometry_error(current, electron):

    args = getargspec(current).args
    kwargs = {}
    if 'V' in args: kwargs['V'] = 0.1

    # Many functions will fail when obtaining the
    # normalization current, which means to test the
    # current function's own fail guard you must not
    # use the normalization current.
    if 'normalization' in args: kwargs['normalization'] = 'thmax'

    with pytest.raises(ValueError):
        current('Bullshit', electron, **kwargs)

@pytest.mark.parametrize("current", current_models)
def test_normalization_error(current, electron):

    args = getargspec(current).args
    kwargs = {}
    if 'eta' in args: kwargs['eta'] = -1
    if 'normalization' in args:
        kwargs['normalization'] = 'Bullshit'

        geo = Cylinder(electron.debye, electron.debye)
        with pytest.raises(ValueError):
            current(geo, electron, **kwargs)
