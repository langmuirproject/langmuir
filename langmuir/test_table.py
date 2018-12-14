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
# These tests are basically duplicates of the ones in test_interpolated.py,
# but these were written first and tests the tables directly, whereas the
# ones in test_interpolated.py tests that the functions using them returns
# the correct values. If the tables are later organized in some other fashion,
# one can probably remove this one.

def test_laframboise_sphere():
    t = get_table('laframboise sphere')
    assert(np.allclose(get_coords_and_value(t, 4, 3), [1, -0.6, 1.595]))

def test_laframboise_cylinder():
    t = get_table('laframboise cylinder')
    assert(np.allclose(get_coords_and_value(t, 5, 6), [3, -2.0, 1.928]))

def test_darian_marholm_uncomplete_sphere():
    t = get_table('darian-marholm uncomplete sphere')

    # Testing one for each alpha/kappa table
    assert(np.allclose(get_coords_and_value(t, 0, 0, 2, 1), [0,   0  , 2, -1, 1.958]))
    assert(np.allclose(get_coords_and_value(t, 0, 1, 2, 1), [0,   0.2, 2, -1, 2.069]))
    assert(np.allclose(get_coords_and_value(t, 1, 0, 2, 1), [1/6, 0  , 2, -1, 1.982]))
    assert(np.allclose(get_coords_and_value(t, 1, 1, 2, 1), [1/6, 0.2, 2, -1, 2.369]))
    assert(np.allclose(get_coords_and_value(t, 2, 0, 2, 1), [1/4, 0  , 2, -1, 1.990]))
    assert(np.allclose(get_coords_and_value(t, 2, 1, 2, 1), [1/4, 0.2, 2, -1, 2.846]))

    # Testing that V=0 is the inaccurate one and that R=0 doesn't exist
    assert(np.allclose(get_coords_and_value(t, 0, 0, 0, 0), [0, 0, 0.2, 0, 0.965]))

def test_darian_marholm_uncomplete_cylinder():
    t = get_table('darian-marholm uncomplete cylinder')

    # Testing one for each alpha/kappa table
    assert(np.allclose(get_coords_and_value(t, 0, 0, 2, 1), [0,   0  , 3, -1, 1.538]))
    assert(np.allclose(get_coords_and_value(t, 0, 1, 2, 1), [0,   0.2, 3, -1, 1.914]))
    assert(np.allclose(get_coords_and_value(t, 1, 0, 2, 1), [1/6, 0  , 3, -1, 1.509]))
    assert(np.allclose(get_coords_and_value(t, 1, 1, 2, 1), [1/6, 0.2, 3, -1, 2.165]))
    assert(np.allclose(get_coords_and_value(t, 2, 0, 2, 1), [1/4, 0  , 3, -1, 1.510]))
    assert(np.allclose(get_coords_and_value(t, 2, 1, 2, 1), [1/4, 0.2, 3, -1, 2.547]))

    # Testing that V=0 is the inaccurate one and that R=0 doesn't exist
    assert(np.allclose(get_coords_and_value(t, 0, 0, 0, 0), [0, 0, 1.0, 0, 0.974]))

def test_darian_marholm_sphere():
    t = get_table('darian-marholm sphere')

    # Testing one from the uncomplete
    assert(np.allclose(get_coords_and_value(t, 1, 0, 3, 1), [1/6, 0  , 2, -1, 1.982]))

    # Testing that V=0 is accurate one and that R=0 exists
    assert(np.allclose(get_coords_and_value(t, 0, 0, 1, 0), [0, 0, 0.2, 0, 1.0]))
    assert(np.allclose(get_coords_and_value(t, 0, 0, 0, 3), [0, 0, 0.0, -3, 4.0]))

def test_darian_marholm_cylinder():
    t = get_table('darian-marholm cylinder')

    # Testing one for each alpha/kappa table
    assert(np.allclose(get_coords_and_value(t, 1, 0, 3, 1), [1/6, 0  , 3, -1, 1.509]))

    # Testing that V=0 is accurate one and that R=0 exists
    assert(np.allclose(get_coords_and_value(t, 0, 0, 1, 0), [0, 0, 1.0, 0, 1.0]))
    assert(np.allclose(get_coords_and_value(t, 0, 0, 0, 3), [0, 0, 0.0, -3, 2.2417], rtol=1e-3))

def test_laframboise_darian_marholm_sphere():
    t = get_table('laframboise-darian-marholm sphere')

    # Testing one from each of Darian-Marholm and Laframboise

    index = map(all, t['points'] == [0.25, 0.2, 3, -5])
    index = np.where(list(index))
    assert(t['values'][index] == approx(4.535))

    index = map(all, t['points'] == [0, 0, 3, -5])
    index = np.where(list(index))
    assert(t['values'][index] == approx(4.640))

def test_laframboise_darian_marholm_cylinder():
    t = get_table('laframboise-darian-marholm cylinder')

    # Testing one from each of Darian-Marholm and Laframboise

    index = map(all, t['points'] == [0.25, 0.2, 3, -5])
    index = np.where(list(index))
    assert(t['values'][index] == approx(3.465))

    index = map(all, t['points'] == [0, 0, 3, -5])
    index = np.where(list(index))
    assert(t['values'][index] == approx(2.701))

def test_invalid_table():
    with pytest.raises(ValueError):
        t = get_table('bullshit')
