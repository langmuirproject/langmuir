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

class Plane(object):
    """
    A plane with specified area A [m^2]
    """
    def __init__(self, A):
        self.A = A
    def __repr__(self):
        return "Plane(A={})".format(self.A)

class Cylinder(object):
    """
    A cylinder with specified radius r [m] and length l [m]
    """
    def __init__(self, r, l):
        self.r = r
        self.l = l
    def __repr__(self):
        return "Cylinder(r={}, l={})".format(self.r, self.l)

class Sphere(object):
    """
    A sphere with specified radius r [m]
    """
    def __init__(self, r):
        self.r = r
    def __repr__(self):
        return "Sphere(r={})".format(self.r)
