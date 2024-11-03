"""
Copyright 2018
    Sigvald Marholm <marholm@marebakken.com>
    Diako Darian <diako.darian@gmail.com>

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
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator
from scipy.constants import value as constants
from copy import deepcopy
import scipy.special as special
import numpy as np
import os

def finite_radius_current(geometry, species, V=None, eta=None, normalization=None,
                          table=None):
    """
    A current model taking into account the effects of finite radius by
    interpolating between tabulated normalized currents. It does not
    extrapolate, but returns ``nan`` when the input parameters are outside the
    domain of the model. This happens when the normalized potential for any
    given species is more than 25, when kappa is less than 4, when alpha is
    more than 0.2 or when the radius is more than 10 or sometimes all the way
    up towards 100. Normally finite radius effects are negligible for radii
    less than 0.2 Debye lengths (spheres) or 1.0 Debye lengths (cylinders).

    One of several tables can be interpolated from, as decided by the ``table``
    argument. Some tables contain values on a regular grid in the parameter
    space, whereas other tables do not. On a regular grid interpolation will be
    linear along each dimension (using SciPy's RegularGridInterpolation), which
    results in faster computation, and in fact, less jagged and more accurate
    resulting curves. For non-regular grids unstructured interpolation through
    SciPy's LinearNDInterpolation is used. Available tables are listed below.
    Before experimenting with non-default settings the user is encouraged to
    read the relevant background literature, and try to understand how the
    tables are made.

    - ``'laframboise'`` (non-regular)
      The de-facto standard tables for finite radius currents, tables 5c
      and 6c in Laframboise, "Theory of Spherical and Cylindrical Langmuir
      Probes in a Collisionless, Maxwellian Plasma at Rest", PhD Thesis.
      Covers Maxwellian plasmas only, and probe radii ranging from 0 to 100
      Debye lengths. This table is believed to contain the numerically most
      accurate values, but because the table has many gaps, unstructured linear
      interpolation is used, which sometimes causes jagged, and unaccurate
      results. Especially near normalized radii of 7.5 or 15 (spheres) or 30 or
      40 (cylinders). Consider ``'laframboise regularized'`` instead.

    - ``'laframboise regularized'`` (regular)
      This fills in the gaps in the ``'laframboise'`` tables by one-dimensional
      interpolation along the radius axis, and subsequently regular grid
      interpolation is used. This fixes the issues with unstructured
      interpolation. It is faster and more accurate than ``'laframboise'``.
      This table will be used by default for Maxwellian plasmas.

    - ``'darian-marholm uncomplete'`` (regular)
      These are the tables reported in Darian, Marholm, Mortensen, Miloch,
      "Theory and simulations of spherical and cylindrical Langmuir probes in
      non-Maxwellian plasmas", PPCF, 2019. These tables covers Maxwellian,
      Kappa, Cairns and Kappa-Cairns distributions for normalized radii ranging
      from 0.2 Debye lengths (spheres) or 1.0 Debye length (cylinders) up to 10
      Debye lengths. They do not extend downwards to 0 Debye lengths. They are
      not as accurate as ``'laframboise'`` for pure Maxwellian distribution,
      but usually within one or two percent. Consider
      ``'laframboise-darian-marholm'`` instead.

    - ``'darian-marholm'`` (regular)
      Same as above, but this is complemented by adding analytical values from
      OML theory for zero radii, thereby extending the range of valid radii
      down to zero Debye lengths. In addition, the values for zero potential
      are replaced by exact analytical values (i.e. the thermal current). These
      are known to be among the most inaccurate in the original
      ``'darian-marholm'`` tables.

    - ``'laframboise-darian-marholm'`` (non-regular)
      This combines the values from ``'laframboise'`` for Maxwellian
      distribution, with the values from ``'darian-marholm'`` for
      non-Maxwellian distributions (i.e., excluding the Maxwellian values in
      the ``darian-marholm'`` table). As such this table spans the largest
      parameter domain, with normalized radii from 0 to somewhere between 10 or
      100, depending on values of the spectral indices. This table suffers from
      inaccurate extrapolation in the same way as the ``'laframboise'`` table.
      In addition, normalized radii over 10 migh be unaccurate unless the
      distribution is nearly Maxwellian. Consider ``'laframboise-darian-marholm
      regularized'`` instead.

    - ``'laframboise-darian-marholm regularized'`` (regular)
      This combines the values in ``'darian-marholm'`` for non-Maxwellian
      distributions with those in ``'laframboise regularized'`` in such a
      manner that the values form a regular grid. This involves discarding
      values for certain radii from the ``'laframboise regularized'`` tables
      which do not exist for ``'darian-marholm'``. In particular normalized
      radii of 7.5 for spheres, and 1.5, 2.5 and 4 for cylinders. In the
      vicinity of these radii, the ``'laframboise'`` tables will be more
      accurate. In addition, normalized radii above 10 are not supported. This
      is the default tables used for non-Maxwellian plasmas.

    Parameters
    ----------
    geometry: Plane, Cylinder or Sphere
        Probe geometry

    species: Species or array-like of Species
        Species constituting the background plasma

    V: float or array-like of floats
        Probe voltage(s) in [V]. Overrides eta.

    eta: float or array-like of floats
        Probe voltage(s) normalized by k*T/q, where q and T are the species'
        charge and temperature and k is Boltzmann's constant.

    normalization: 'th', 'thmax', 'oml', None
        Wether to normalize the output current by, respectively, the thermal
        current, the Maxwellian thermal current, the OML current, or not at
        all, i.e., current in [A/m].

    table: string
        Which table to use for interpolation. See detailed description above.

    Returns
    -------
    float if voltage is float. array of floats corresponding to voltage if
    voltage is array-like.
    """
    if isinstance(species, list):
        if normalization is not None:
            logger.error('Cannot normalize current to more than one species')
            return None
        if eta is not None:
            logger.error('Cannot normalize voltage to more than one species')
            return None
        I = 0
        for s in species:
            I += finite_radius_current(geometry, s, V, eta, table=table)
        return I

    q, m, n, T = species.q, species.m, species.n, species.T
    kappa, alpha = species.kappa, species.alpha
    k = constants('Boltzmann constant')

    if V is not None:
        V = make_array(V)
        eta = -q*V/(k*T)
    else:
        eta = make_array(eta)

    if table is None:
        if alpha==0 and kappa==np.inf:
            table = 'laframboise regularized'
        else:
            table = 'laframboise-darian-marholm regularized'

    eta = deepcopy(eta)

    I = np.zeros_like(eta)

    indices_n = np.where(eta < 0)[0]   # indices for repelled particles
    indices_p = np.where(eta >= 0)[0]  # indices for attracted particles

    if normalization is None:
        I0 = normalization_current(geometry, species)
    elif normalization.lower() == 'thmax':
        I0 = 1
    elif normalization.lower() == 'th':
        I0 = normalization_current(geometry, species)/\
             thermal_current(geometry, species)
    elif normalization.lower() == 'oml':
        I0 = normalization_current(geometry, species)/\
             OML_current(geometry, species, eta=eta)
    else:
        raise ValueError('Normalization not supported: {}'.format(normalization))

    if isinstance(geometry, Sphere):
        table += ' sphere'
    elif isinstance(geometry, Cylinder):
        table += ' cylinder'
    else:
        raise ValueError('Geometry not supported: {}'.format(geometry))

    R = geometry.r/species.debye

    lerp = LERPERS[table]

    if "darian-marholm" in table:
        I[indices_p] = I0*lerp((1/kappa, alpha, R, eta[indices_p]))
    else:
        I[indices_p] = I0*lerp((R, eta[indices_p]))
        if(kappa != float('inf') or alpha != 0):
            logger.warning("Using pure Laframboise tables discards spectral indices kappa and alpha")

    I[indices_n] = I0*OML_current(geometry, species, eta=eta[indices_n], normalization='thmax')

    if any(np.isnan(I)):
        logger.warning("Data points occurred outside the domain of tabulated values resulting in nan")

    return I[0] if len(I) == 1 else I

def get_fr_lerper(table_name):
    """
    Returns a linear interpolator. For non-regular data, it will return an
    LinearNDInterpolator object (triangulates the parameter space). For regular
    data it will return a RegularGridInterpolator object.
    """

    table = get_table(table_name)
    if table['regular']:
        pts = table['points']
        axes = table['axes']
        vals = table['values']
        lerper = RegularGridInterpolator(axes, vals, bounds_error=False, fill_value=np.nan)
    else:
        pts = table['points']
        vals = table['values'].reshape(-1)
        lerper = LinearNDInterpolator(pts, vals)

    return lerper

# Creating these lerpers as global constant objects means that they will only
# be initialized once, which improves performance.
LERPERS = {table_name: get_fr_lerper(table_name) for table_name in TABLE_NAMES}
