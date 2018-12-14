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
import inspect
from scipy.constants import value as constants

current_models = [
    thermal_current,
    normalization_current,
    OML_current,
    finite_radius_current
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
    args = inspect.getargspec(current).args
    kwargs = {}
    if 'V' in args: kwargs['V'] = 0.1

    geo = Sphere(electron.debye)
    Ie = current(geo, electron, **kwargs)
    Ii = current(geo, proton, **kwargs)
    I  = current(geo, [electron, proton], **kwargs)

    assert(Ie+Ii == approx(I))

@pytest.mark.parametrize("current", current_models)
def test_normalize(current, electron):

    if current != normalization_current:

        args = inspect.getargspec(current).args
        kwargs = {}
        if 'V' in args: kwargs['V'] = 0.1

        In = current(Sphere(electron.debye), electron, normalize=True, **kwargs)
        I  = current(Sphere(electron.debye), electron, **kwargs)
        I0 = normalization_current(Sphere(electron.debye), electron)

        assert(I0*In == approx(I))

@pytest.mark.parametrize("current", current_models)
def test_eta(current, electron):

    args = inspect.getargspec(current).args
    kwargs = {}
    if 'V' in args or 'eta' in args:

        # Should have both or none
        assert 'V' in args
        assert 'eta' in args

        e = constants('elementary charge')
        k = constants('Boltzmann constant')
        T = 1000
        eta = -10

        I_eta = current(Sphere(electron.debye), electron, eta=eta)
        I_V   = current(Sphere(electron.debye), electron, V=-eta*k*T/e)

        assert I_eta == approx(I_V)

@pytest.mark.parametrize("current", current_models)
def test_multiple_species_eta(current, electron, proton, caplog):

    args = inspect.getargspec(current).args
    if 'eta' in args:

        I = current(Sphere(electron.debye), [electron, proton], eta=1)
        assert(caplog.records[0].levelname == 'ERROR')

@pytest.mark.parametrize("current", current_models)
def test_multiple_species_normalize(current, electron, proton, caplog):

    args = inspect.getargspec(current).args
    kwargs = {}
    if 'V' in args: kwargs['V'] = 0.1
    if 'normalize' in args:

        I = current(Sphere(electron.debye), [electron, proton],
                    normalize=True, **kwargs)
        assert(caplog.records[0].levelname == 'ERROR')

@pytest.mark.parametrize("current", current_models)
def test_input_output_format(current, electron):

    args = inspect.getargspec(current).args
    if 'V' in args:

        I = current(Sphere(electron.debye), electron, 0.1)
        assert(isinstance(I, float))

        I = current(Sphere(electron.debye), electron, 1)
        assert(isinstance(I, float))

        I = current(Sphere(electron.debye), electron, [1, 2])
        assert(isinstance(I, np.ndarray))
        assert(I.dtype==np.float)
        assert(I.shape==(2,))

        I = current(Sphere(electron.debye), electron, np.array([1, 2]))
        assert(isinstance(I, np.ndarray))
        assert(I.dtype==np.float)
        assert(I.shape==(2,))

@pytest.mark.parametrize("current", current_models)
def test_geometry_error(current, electron):

    args = inspect.getargspec(current).args
    kwargs = {}
    if 'V' in args: kwargs['V'] = 0.1

    # Many functions will fail when obtaining the
    # normalization current, which means to test the
    # current function's own fail guard you must not
    # use the normalization current.
    if 'normalize' in args: kwargs['normalize'] = True

    with pytest.raises(ValueError):
        current('Bullshit', electron, **kwargs)
