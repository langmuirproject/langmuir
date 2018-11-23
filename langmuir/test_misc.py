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

def test_make_array():
    assert(isinstance(make_array(2), np.ndarray))
    assert(isinstance(make_array((2,3)), np.ndarray))
    assert(isinstance(make_array([2,3]), np.ndarray))
    assert(np.allclose(2, make_array(  int(2))))
    assert(np.allclose(2, make_array(float(2))))
    assert(np.allclose([10, 1/10], make_array([10, 1/10])))
