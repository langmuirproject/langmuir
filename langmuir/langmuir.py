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

import numpy as np
import logging
from langmuir.tables import *
from scipy.interpolate import griddata
from scipy.constants import value as constants
from scipy.special import gamma, erfc, hyp2f1
from copy import deepcopy

logger = logging.getLogger('langmuir')
logging.basicConfig()

class Species(object):
    """
    Defines a species using a set of flags and keyword parameters. Maxwellian
    electrons are the default but other distributions can be specified using a
    flag, e.g. 'kappa'. Other parameters can be specified/overridden using the
    keyword parameters. Examples::

        >>> # Maxwellian electrons with density 1e11 m^(-3) and 1000 K
        >>> species = Species(n=1e11, T=1000)

        >>> # Kappa-distributed electrons with kappa=2
        >>> species = Species('kappa', n=1e11, T=1000, kappa=2)

        >>> # Doubly charged Maxwellain Oxygen ions with 0.26 eV
        >>> species = Species(Z=2, amu=16, n=1e11, eV=0.26)

        >>> # Cairns-distributed protons with alpha=0.2
        >>> species = Species('cairns', 'proton', n=1e11, T=1000, alpha=0.2)

    A plasma is fully specified by a list of species. E.g. for a Maxwellian
    electron-proton plasma::

        >>> plasma = []
        >>> plasma.append(Species('electron', n=1e11, T=1000))
        >>> plasma.append(Species('proton',   n=1e11, T=1000))

    Flags:
    ------
    'maxwellian'   : Maxwellian distribution (default)
    'kappa'        : Kappa distribution
    'cairns'       : Cairns distribution
    'kappa-cairns' : Kappa-Cairns distribution
    'electron'     : Elecron species (default)
    'proton'       : Proton species
    'positron'     : Positron species

    Keyword parameters:
    -------------------
    m     : mass [kg]
    amu   : mass [amu]
    q     : charge [C]
    Z     : charge [elementary charges]
    n     : density [m^(-3)]
    T     : temperature [K]
    eV    : temperature [eV]
    vth   : thermal velocity [m/s]
    kappa : spectral index kappa (for Kappa and Kappa-Cairns distributions)
    alpha : spectral index alpha (for Cairns and Kappa-Cairns distributions)
    """

    def __init__(self, *args, **kwargs):

        e    = constants('elementary charge')
        me   = constants('electron mass')
        mp   = constants('proton mass')
        k    = constants('Boltzmann constant')
        eps0 = constants('electric constant')
        amu  = constants('atomic mass constant')

        self.dist = 'maxwellian'
        valid_dists = ['maxwellian', 'kappa', 'cairns', 'kappa-cairns']

        if 'Z' in kwargs:
            self.q = kwargs['Z']*e
        else:
            self.q = kwargs.pop('q', -e)

        if 'amu' in kwargs:
            self.m = kwargs['amu']*amu
        else:
            self.m = kwargs.pop('m', me)

        for arg in args:
            if arg.lower() in valid_dists:
                self.dist = arg.lower()

            if arg.lower()=='proton':
                self.q = e
                self.m = mp

            if arg.lower()=='positron':
                self.q = e
                self.m = me

        self.n = kwargs['n']

        if 'eV' in kwargs:
            self.T = kwargs['eV']*e/k
        elif 'vth' in kwargs:
            self.T = self.m*kwargs['vth']**2/k
        else:
            self.T = kwargs['T']

        self.vth = np.sqrt(k*self.T/self.m)

        if self.dist == 'maxwellian':
            self.alpha = 0
            self.kappa = float('inf')
            self.debye = np.sqrt(eps0*k*self.T/(self.q**2*self.n))

        if self.dist == 'kappa':
            self.alpha = 0
            self.kappa = kwargs['kappa']
            self.debye = np.sqrt(eps0*k*self.T/(self.q**2*self.n))*np.sqrt((self.kappa-1.5)/(self.kappa-0.5))

        if self.dist == 'cairns':
            self.alpha = kwargs['alpha']
            self.kappa = float('inf')
            self.debye = np.sqrt(eps0*k*self.T/(self.q**2*self.n))*np.sqrt((1.0 + 15.0*self.alpha)/(1.0 + 3.0*self.alpha))

        if self.dist == 'kappa-cairns':
            self.alpha = kwargs['alpha']
            self.kappa = kwargs['kappa']
            self.debye = np.sqrt(eps0 * k * self.T / (self.q**2 * self.n)) *\
                         np.sqrt( ((self.kappa - 1.5) / (self.kappa - 0.5)) *\
                                  ( (1.0 + 15.0 * self.alpha * ((self.kappa - 1.5) / (self.kappa - 2.5))) / (1.0 + 3.0 * self.alpha * ((self.kappa - 1.5) / (self.kappa - 0.5)))))

    def __repr__(self):

        s = "Species('{}', q={}, m={}, n={}, T={}".format(self.dist, self.q, self.m, self.n, self.T)

        if(self.dist == 'kappa' or self.dist == 'kappa-cairns'):
            s += ", kappa={}".format(self.kappa)

        if(self.dist == 'cairns' or self.dist == 'kappa-cairns'):
            s += ", alpha={}".format(self.alpha)

        s += ")"

        return s


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

def normalization_current(geometry, species):
    """
    Returns the normalization current for the given species and geometry.
    The normalization current is the current the species would have contributed
    to a neutral probe due to random thermal movements of particles if the species
    had been Maxwellian.

    Parameters:
    -----------
    geometry  (Plane, Cylinder, Sphere)       : geometry of the probe
    species   (Species, list of Species)      : plasma species
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

def thermal_current(geometry, species):
    """
    Returns the thermal current for the given species and geometry. The
    thermal current is the current the species contributes to a neutral
    probe due to random thermal movements of particles.

    Parameters:
    -----------
    geometry  (Plane, Cylinder, Sphere)       : geometry of the probe
    species   (Species, list of Species)      : plasma species
    """

    if isinstance(species, list):
        I = 0
        for s in species:
            I += thermal_current(geometry, s)
        return I

    raise NotImplementedError

def make_array(arr):
    if isinstance(arr, (int, float)):
        arr = np.array([arr], dtype=np.float)
    elif isinstance(arr, (list, tuple)):
        arr = np.array(arr, dtype=np.float)
    return arr

def OML_current(geometry, species, V=None, eta=None, normalize=False):

    """
    Returns the OML current.

    Parameters:
    -----------
    geometry  (Plane, Cylinder, Sphere)       : geometry of the probe
    species   (Species, list of Species)      : plasma species
    V  (int, float, list, tuple, numpy array) : probe's biased voltage
    normalize (bool)                          : normalize current with respect to random/thermal current, a la Laframboise
    """

    if isinstance(species, list):
        if normalize == True:
            logger.error('Cannot normalize current to more than one species')
            return None
        if eta is not None:
            logger.error('Cannot normalize voltage to more than one species')
            return None
        I = 0
        for s in species:
            I += OML_current(geometry, s, V, eta)
        return I

    q, m, n, T = species.q, species.m, species.n, species.T
    kappa, alpha = species.kappa, species.alpha
    vth = species.vth
    k = constants('Boltzmann constant')

    if V is not None:
        V = make_array(V)
        eta = q*V/(k*T)
    else:
        eta = make_array(eta)

    eta = deepcopy(eta) # Prevents this function from changing eta in caller

    I = np.zeros_like(eta)

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

    if normalize:
        I0 = 1
    else:
        I0 = normalization_current(geometry, species)

    if isinstance(geometry, Sphere):
        r = geometry.r

        # repelled particles:
        if species.dist == 'maxwellian' or species.dist == 'cairns':
            I[indices_n] = I0 * C * D * np.exp(-eta[indices_n]) * (1. + E * eta[indices_n] * (eta[indices_n] + 4.))

        elif species.dist == 'kappa' or species.dist == 'kappa-cairns':
            I[indices_n] = I0 * C * D * (1. + eta[indices_n] / (kappa - 1.5))**(
                1. - kappa) * (1. + E * eta[indices_n] * (eta[indices_n] + 4. * ((kappa - 1.5) / (kappa - 1.))))

        # attracted particles:
        eta[indices_p] = np.abs(eta[indices_p])
        I[indices_p] = I0*C*D*(1.+F*eta[indices_p])

    elif isinstance(geometry, Cylinder):
        r, l = geometry.r, geometry.l

        # repelled particles:
        if species.dist == 'maxwellian' or species.dist == 'cairns':
            I[indices_n] = I0 * C * D * \
                np.exp(-eta[indices_n]) * (1. + E *
                                           eta[indices_n] * (eta[indices_n] + 4.))

        elif species.dist == 'kappa' or species.dist == 'kappa-cairns':
            I[indices_n] = I0 * C * D * (1. + eta[indices_n] / (kappa - 1.5))**(
                1. - kappa) * (1. + E * eta[indices_n] * (eta[indices_n] + 4. * ((kappa - 1.5) / (kappa - 1.))))

        # attracted particles:
        eta[indices_p] = np.abs(eta[indices_p])
        if species.dist == 'maxwellian' or species.dist == 'cairns':
            I[indices_p] = I0 * C * D * \
                           ((2. / np.sqrt(np.pi)) * ( 1 - 0.5 * E * eta[indices_p]) * np.sqrt(eta[indices_p]) + \
                           np.exp(eta[indices_p]) * (1. + E * eta[indices_p] * (eta[indices_p] - 4.)) * erfc(np.sqrt(eta[indices_p])))

        elif species.dist == 'kappa' or species.dist == 'kappa-cairns':
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

    return I

def finite_radius_current(geometry, species, V=None, eta=None, table='laframboise-darian-marholm', normalize=False, flip_kappa=False):

    if isinstance(species, list):
        if normalize == True:
            logger.error('Cannot normalize current to more than one species')
            return None
        if eta is not None:
            logger.error('Cannot normalize voltage to more than one species')
            return None
        I = 0
        for s in species:
            I += tabulated_current(geometry, s, V, eta, table)
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

    R = geometry.r/species.debye

    if normalize:
        I0 = 1
    else:
        I0 = normalization_current(geometry, species)

    if isinstance(geometry, Sphere):
        table += ' sphere'
    elif isinstance(geometry, Cylinder):
        table += ' cylinder'

    if "darian-marholm" in table:
        table = get_table(table, flip_kappa=flip_kappa)
        pts = table['points']
        vals = table['values'].reshape(-1)
        if flip_kappa:
            kappa = 1/kappa
        else:
            kappa = np.min([kappa, 1000])
        print(kappa, alpha, R, eta, -eta[indices_p])
        I[indices_p] = I0*griddata(pts, vals, (kappa, alpha, R, -eta[indices_p]))
        print(I)

    else:
        table = get_table(table)
        pts = table['points']
        vals = table['values'].reshape(-1)
        I[indices_p] = I0*griddata(pts, vals, (R, -eta[indices_p]))
        if(kappa != float('inf') or alpha != 0):
            logger.warning("Using pure Laframboise tables discards spectral indices kappa and alpha")

    if len(indices_n)>0:
        pos_neg = "positive" if q>0 else "negative"
        logger.warning("Only attracted species current is covered by tabulated "
                       "values. Currents due to {} is set to zero for "
                       "{} potentials".format(species, pos_neg))

    if any(np.isnan(I)):
        logger.warning("Data points occurred outside the domain of tabulated values resulting in nan")

    return I

# def current(geometry, species, V):

#     if isinstance(species, list):
#         I = 0
#         for s in species:
#             I += current(geometry, species, V)
#         return I

#     """
#     if R=0:
#         call OML function for correct distribution
#         It should take care of geometry, attracted/repelled
#     elif R=inf:
#         thin sheath
#     else:
#         table
#     """

