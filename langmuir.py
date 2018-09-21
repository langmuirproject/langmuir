"""
Copyright 2018 Sigvald Marholm <marholm@marebakken.com>

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

from scipy.interpolate import interp2d

def lafr_func(geometry, kind='linear'):
    """
    Returns normalized current per Laframboise's table.
    """

    if geometry=='Sphere':
        Rs =  [0, 0.2, 0.3, 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 50, 100]

        Vs =  [0.0  , 0.1  , 0.3  , 0.6  , 1.0  , 1.5  , 2.0  , 3.0  , 5.0  , 7.5  , 10.0  , 15.0  , 20.0  , 25.0]
        Is = [[1.000, 1.100, 1.300, 1.600, 2.000, 2.500, 3.000, 4.000, 6.000, 8.500, 11.000, 16.000, 21.000, 26.000],   # R=0
              [1.000, 1.100, 1.300, 1.600, 2.000, 2.500, 3.000, 4.000, 6.000, 8.500, 11.000, 16.000, 21.000, 25.763],   # R=0.2
              [1.000, 1.100, 1.300, 1.600, 2.000, 2.500, 3.000, 4.000, 6.000, 8.500, 11.000, 16.000, 21.000, 25.462],   # R=0.3
              [1.000, 1.100, 1.300, 1.600, 2.000, 2.493, 2.987, 3.970, 5.917, 8.324, 10.704, 15.403, 20.031, 24.607],   # R=0.5
              [1.000, 1.0999,1.299, 1.595, 1.987, 2.469, 2.945, 3.878, 5.687, 7.871,  9.990, 14.085, 18.041, 21.895],   # R=1
              [1.000, 1.0999,1.299, 1.584, 1.955, 2.399, 2.824, 3.632, 5.126, 6.847,  8.460, 11.482, 14.314, 17.018],   # R=2
              [1.000, 1.0999,1.293, 1.572, 1.922, 2.329, 2.709, 3.406, 4.640, 6.007,  7.258,  9.542, 11.636, 13.603],   # R=3
              [1.000, 1.099, 1.288, 1.552, 1.869, 2.219, 2.529, 3.086, 3.957, 4.887,  5.710,  7.167,  8.473,  9.676],   # R=5
              [1.000, 1.099, 1.288, 1.552, 1.869, 2.219, 2.529, 3.086, 3.957, 4.094,  4.658,  5.645,  6.518,  7.318],   # R=7.5
              [1.000, 1.098, 1.280, 1.518, 1.783, 2.050, 2.226, 2.609, 3.119, 3.620,  4.050,  4.796,  5.453,  6.053],   # R=10
              [1.000, 1.098, 1.280, 1.518, 1.783, 2.050, 2.226, 2.609, 3.119, 3.620,  4.050,  4.796,  4.318,  4.719],   # R=15
              [1.000, 1.097, 1.269, 1.481, 1.694, 1.887, 2.030, 2.235, 2.516, 2.779,  3.002,  3.383,  3.716,  4.018],   # R=20
              [1.000, 1.095, 1.255, 1.433, 1.592, 1.719, 1.803, 1.910, 2.037, 2.148,  2.241,  2.397,  2.532,  2.658],   # R=50
              [1.000, 1.094, 1.245, 1.402, 1.534, 1.632, 1.694, 1.762, 1.833, 1.891,  1.938,  2.022,  2.097,  2.166]]   # R=50
    else:
        raise ValueError('Geometry {} not implemented'.format(geometry))

    f = interp2d(Vs, Rs, Is, kind=kind)
    return lambda R, V: float(f(V, R))
