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
from langmuir.species import *
from langmuir.geometry import *
from langmuir.misc import *
from scipy.constants import value as constants
from scipy.special import gamma, erfc, hyp2f1
from copy import deepcopy
import numpy as np

def normalization_current(geometry, species):
    """
    Returns the normalization current for the given species and geometry.
    The normalization current is the current the species would have contributed
    to a probe at zero potential with respect to the background plasma due to
    random thermal movements of particles, if the species had been Maxwellian.

    Parameters
    ----------
    geometry: Plane, Cylinder or Sphere
        Probe geometry

    species: Species or array-like of Species
        Species constituting the background plasma

    Returns
    -------
    float
    """

    if isinstance(species, list):
        I = 0
        for s in species:
            I += normalization_current(geometry, s)
        return I

    q, n, vth = species.q, species.n, species.vth

    if isinstance(geometry, Sphere):
        r = geometry.r
        I0 = 2*np.sqrt(2*np.pi)*r**2*q*n*vth

    elif isinstance(geometry, Cylinder):
        r, l = geometry.r, geometry.l
        I0 = np.sqrt(2*np.pi)*r*l*q*n*vth

    else:
        raise ValueError('Geometry not supported: {}'.format(geometry))

    return I0

def thermal_current(geometry, species, normalization=None):
    """
    Returns the thermal current.

    Parameters
    ----------
    geometry: Plane, Cylinder or Sphere
        Probe geometry

    species: Species or array-like of Species
        Species constituting the background plasma

    Returns
    -------
    float
    """

    if isinstance(species, list):
        if normalization is not None:
            logger.error('Cannot normalize current to more than one species')
            return None
        I = 0
        for s in species:
            I += thermal_current(geometry, s)
        return I

    kappa, alpha = species.kappa, species.alpha

    if normalization is None:
        I0 = normalization_current(geometry, species)
    elif normalization.lower() == 'thmax':
        I0 = 1
    else:
        raise ValueError('Normalization not supported: {}'.format(normalization))

    if isinstance(geometry, Sphere):
        if kappa == float('inf'):
            C = 1.0
            D = (1.+24*alpha)/(1.+15*alpha)
        else:
            C = np.sqrt(kappa-1.5)*gamma(kappa-1.)/gamma(kappa-0.5)
            D = (1.+24*alpha*((kappa-1.5)**2/((kappa-2.)*(kappa-3.))))/(1.+15*alpha*((kappa-1.5)/(kappa-2.5)))
    elif isinstance(geometry, Cylinder):
        if kappa == float('inf'):
            C = 1.0
            D = (1.+3.*alpha)/(1.+15.*alpha)
        else:
            C = np.sqrt(kappa - 1.5) * (kappa - .5) / (kappa - 1.0)
            D = (1. + 3 * alpha * ((kappa - 1.5) / (kappa - 0.5))) / \
                (1. + 15 * alpha * ((kappa - 1.5) / (kappa - 2.5)))
    else:
        raise ValueError('Geometry not supported: {}'.format(geometry))

    return I0*C*D

def OML_current(geometry, species, V=None, eta=None, normalization=None):
    """
    Current collected by a probe according to the Orbital Motion Limited (OML)
    theory. The model assumes a probe of infinitely small radius compared to
    the Debye length, and for a cylindrical probe, that it is infinitely long.
    Probes with radii up to 0.2 Debye lengths (for spherical probes) or 1.0
    Debye lengths (for cylindrical probes) are very well approximated by this
    theory, although the literature is diverse as to how long cylindrical probes
    must be for this theory to be a good approximation.

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

    normalization: 'th', 'thmax', None
        Wether to normalize the output current by, respectively, the thermal
        current, the Maxwellian thermal current, or not at all, i.e., current
        in [A/m].

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
            I += OML_current(geometry, s, V, eta)
        return I

    q, T = species.q, species.T
    kappa, alpha = species.kappa, species.alpha

    k = constants('Boltzmann constant')

    if V is not None:
        eta_isarray = isinstance(V, (np.ndarray, list, tuple))
        V = make_array(V)
        eta = q*V/(k*T)
    else:
        eta_isarray = isinstance(eta, (np.ndarray, list, tuple))
        eta = make_array(eta)

    eta = deepcopy(eta) # Prevents this function from changing eta in caller

    I = np.zeros_like(eta, dtype=float)

    eta[np.where(eta==0)] = np.finfo(float).eps # replace zeros with machine eps
    indices_n = np.where(eta > 0)[0]   # indices for repelled particles
    indices_p = np.where(eta <= 0)[0]  # indices for attracted particles

    if kappa == float('inf'):
        C = 1.0
        D = (1.+24*alpha)/(1.+15*alpha)
        E = 4.*alpha/(1.+24*alpha)
        F = (1.+8*alpha)/(1.+24*alpha)
    else:
        C = np.sqrt(kappa-1.5)*gamma(kappa-1.)/gamma(kappa-0.5)
        D = (1.+24*alpha*((kappa-1.5)**2/((kappa-2.)*(kappa-3.))))/(1.+15*alpha*((kappa-1.5)/(kappa-2.5)))
        E = 4.*alpha*kappa*(kappa-1.)/( (kappa-2.)*(kappa-3.)+24*alpha*(kappa-1.5)**2 )
        F = ((kappa-1.)*(kappa-2.)*(kappa-3.)+8*alpha*(kappa-3.)*(kappa-1.5)**2) /( (kappa-2.)*(kappa-3.)*(kappa-1.5)+24*alpha*(kappa-1.5)**3 )

    if normalization is None:
        I0 = normalization_current(geometry, species)
    elif normalization.lower() == 'thmax':
        I0 = 1
    elif normalization.lower() == 'th':
        I0 = normalization_current(geometry, species)/\
             thermal_current(geometry, species)
    else:
        raise ValueError('Normalization not supported: {}'.format(normalization))

    if isinstance(geometry, Sphere):

        # repelled particles:
        if species.kappa == float('inf'):
            I[indices_n] = I0 * C * D * np.exp(-eta[indices_n]) * (1. + E * eta[indices_n] * (eta[indices_n] + 4.))

        else:
            I[indices_n] = I0 * C * D * (1. + eta[indices_n] / (kappa - 1.5))**(
                1. - kappa) * (1. + E * eta[indices_n] * (eta[indices_n] + 4. * ((kappa - 1.5) / (kappa - 1.))))

        # attracted particles:
        eta[indices_p] = np.abs(eta[indices_p])
        I[indices_p] = I0*C*D*(1.+F*eta[indices_p])

    elif isinstance(geometry, Cylinder):

        # repelled particles:
        if species.kappa == float('inf'):
            I[indices_n] = I0 * C * D * \
                np.exp(-eta[indices_n]) * (1. + E *
                                           eta[indices_n] * (eta[indices_n] + 4.))

        else:
            I[indices_n] = I0 * C * D * (1. + eta[indices_n] / (kappa - 1.5))**(
                1. - kappa) * (1. + E * eta[indices_n] * (eta[indices_n] + 4. * ((kappa - 1.5) / (kappa - 1.))))

        # attracted particles:
        eta[indices_p] = np.abs(eta[indices_p])
        if species.kappa == float('inf'):
            I[indices_p] = I0 * C * D * \
                           ((2. / np.sqrt(np.pi)) * ( 1 - 0.5 * E * eta[indices_p]) * np.sqrt(eta[indices_p]) + \
                           np.exp(eta[indices_p]) * (1. + E * eta[indices_p] * (eta[indices_p] - 4.)) * erfc(np.sqrt(eta[indices_p])))

        else:
            C = np.sqrt(kappa - 1.5) * (kappa - .5) / (kappa - 1.0)
            D = (1. + 3 * alpha * ((kappa - 1.5) / (kappa - 0.5))) / \
                (1. + 15 * alpha * ((kappa - 1.5) / (kappa - 2.5)))
            E = 4. * alpha * kappa * \
                (kappa - 1.) / ((kappa - .5) *
                                (kappa - 1.5) + 3. * alpha * (kappa - 1.5)**2)

            I[indices_p] = (2./np.sqrt(np.pi))*I0 * C * D * (eta[indices_p]/(kappa-1.5))**(1.-kappa) * \
                (((kappa - 1.) / (kappa - 3.)) * E * (eta[indices_p]**2) * hyp2f1(kappa - 3, kappa + .5, kappa - 2., 1. - (kappa - 1.5) / (eta[indices_p])) + \
                ((kappa - 1.5 - 2. * (kappa - 1.) * eta[indices_p]) / (kappa - 2.)) * E * eta[indices_p] * hyp2f1(kappa - 2, kappa + .5, kappa - 1., 1. - (kappa - 1.5) / (eta[indices_p])) +
                (1. + E * eta[indices_p] * (eta[indices_p]-((kappa-1.5)/(kappa-1.)))) * hyp2f1(kappa - 1., kappa + .5, kappa, 1. - (kappa - 1.5) / (eta[indices_p])))

    else:
        raise ValueError('Geometry not supported: {}'.format(geometry))

    return I if eta_isarray else I[0]
