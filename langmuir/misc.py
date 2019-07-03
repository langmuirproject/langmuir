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
import logging
import numpy as np

logger = logging.getLogger('langmuir')
logging.basicConfig()

def make_array(arr):
    """
    Takes a list, tuple, integer or float and returns a numpy array.
    """
    if isinstance(arr, (int, float, np.int64, np.int32, np.float64, np.float32)):
        arr = [arr]
    return np.array(arr, dtype=np.float)

