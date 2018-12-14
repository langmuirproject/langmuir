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
from langmuir.tables import *
from langmuir.geometry import *
from langmuir.species import *
from langmuir.misc import *
from scipy.interpolate import griddata
from scipy.constants import value as constants
from copy import deepcopy
import numpy as np

def finite_radius_current(geometry, species, V=None, eta=None, table='laframboise-darian-marholm', normalize=False):

    if isinstance(species, list):
        if normalize == True:
            logger.error('Cannot normalize current to more than one species')
            return None
        if eta is not None:
            logger.error('Cannot normalize voltage to more than one species')
            return None
        I = 0
        for s in species:
            I += finite_radius_current(geometry, s, V, eta, table)
        return I

    q, m, n, T = species.q, species.m, species.n, species.T
    kappa, alpha = species.kappa, species.alpha
    k = constants('Boltzmann constant')

    if V is not None:
        V = make_array(V)
        eta = q*V/(k*T)
    else:
        eta = make_array(eta)

    eta = deepcopy(eta)

    I = np.zeros_like(eta)

    indices_n = np.where(eta > 0)[0]   # indices for repelled particles
    indices_p = np.where(eta <= 0)[0]  # indices for attracted particles

    if normalize:
        I0 = 1
    else:
        I0 = normalization_current(geometry, species)

    if isinstance(geometry, Sphere):
        table += ' sphere'
    elif isinstance(geometry, Cylinder):
        table += ' cylinder'
    else:
        raise ValueError('Geometry not supported: {}'.format(geometry))

    R = geometry.r/species.debye

    if "darian-marholm" in table:
        table = get_table(table)
        pts = table['points']
        vals = table['values'].reshape(-1)
        I[indices_p] = I0*griddata(pts, vals, (1/kappa, alpha, R, eta[indices_p]))

    else:
        table = get_table(table)
        pts = table['points']
        vals = table['values'].reshape(-1)
        I[indices_p] = I0*griddata(pts, vals, (R, eta[indices_p]))
        if(kappa != float('inf') or alpha != 0):
            logger.warning("Using pure Laframboise tables discards spectral indices kappa and alpha")

    if len(indices_n)>0:
        pos_neg = "positive" if q>0 else "negative"
        logger.warning("Only attracted species current is covered by tabulated "
                       "values. Currents due to {} is set to zero for "
                       "{} potentials".format(species, pos_neg))

    if any(np.isnan(I)):
        logger.warning("Data points occurred outside the domain of tabulated values resulting in nan")

    return I[0] if len(I) == 1 else I
