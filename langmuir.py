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
from scipy.interpolate import interp2d
from scipy.constants import value as constants
from scipy.special import gamma, erfc, hyp2f1

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

        >>> # Doubly charged Maxwellain Oxygen ions # NOT IMPLEMENTED
        >>> # species = Species(Z=2, amu=16, n=1e11, T=1000)

        >>> # Cairns-distributed protons with alpha=0.2
        >>> species = Species('cairns', 'proton', n=1e11, T=1000, alpha=0.2)

    A plasma is fully specified by a list of species. E.g. for an
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

    Keyword parameters:
    -------------------
    m     : mass [kg]
    amu   : mass [amu] # NOT IMPLEMENTED
    q     : charge [C]
    Z     : charge [elementary charges] # NOT IMPLEMENTED
    n     : density [m^{-3}]
    T     : temperature [K]
    eV    : temperature [eV] # NOT IMPLEMENTED
    vth   : thermal velocity [m/s] # NOT IMPLEMENTED
    kappa : kappa (for Kappa and Kappa-Cairns distributions)
    alpha : alpha (for Cairns and Kappa-Cairns distributions)
    """

    def __init__(self, *args, **kwargs):

        dist = 'maxwellian'

        for arg in args:
            if arg.lower() in ['maxwellian',
                               'kappa',
                               'cairns',
                               'kappa-cairns']:
                dist = arg.lower()

            if arg.lower()=='proton':
                kwargs['q']=constants('elementary charge')
                kwargs['m']=constants('proton mass')

        self.dist = dist
        self.q = kwargs.pop('q', -constants('elementary charge'))
        self.m = kwargs.pop('m', constants('electron mass'))
        self.n = kwargs.pop('n')

        if dist == 'maxwellian':
            self.T     = kwargs.pop('T')
            self.alpha = 0
            self.kappa = float('inf')

        if dist == 'kappa':
            self.T     = kwargs.pop('T')
            self.alpha = 0
            self.kappa = kwargs.pop('kappa')

        if dist == 'cairns':
            self.T     = kwargs.pop('T')
            self.alpha = kwargs.pop('alpha')
            self.kappa = float('inf')

        if dist == 'kappa-cairns':
            self.T     = kwargs.pop('T')
            self.alpha = kwargs.pop('alpha')
            self.kappa = kwargs.pop('kappa')

class Geometry(object):
    """
    Specifies a probe geometry. Examples::

        >>> # Plane probe with 0.01 m^2 area
        >>> geometry = Geometry('plane', A=0.01)

        >>> # Cylindrical probe with 1e-3 m radius and 25e-3 m length
        >>> geometry = Geometry('cylinder', r=1e-3, l=25e-3)

        >>> # Spherical probe with 1e-2 m radius
        >>> geometry = Geometry('sphere', r=0.01)

    Flags:
    ------
    'plane'    : Plane probe
    'cylinder' : Cylindrical probe (default)
    'sphere'   : Spherical probe

    Keyword parameters:
    -------------------
    r : radius [m]
    l : length [m]
    A : area [m^2]
    """
    def __init__(self, shape='cylinder', **kwargs):
        self.shape = shape.lower()
        assert self.shape in ['plane', 'sphere', 'cylinder']

        if shape == 'plane':
            self.A = kwargs.pop('A')

        if shape == 'sphere':
            self.r = kwargs.pop('r')

        if shape == 'cylinder':
            self.r = kwargs.pop('r')
            self.l = kwargs.pop('l')

def OML_current(geometry, species, V, normalize=True):

    """
    Returns the OML current.

    Parameters:
            geometry  (Geometry): geometry of the probe (spherical, cylindrical or planar)
            species   (Species) : plasma species (electron or ion)
            V  (int, float, list, tuple, numpy array) : probe's biased voltage
            normalize (bool) : normalize current with respect to random/thermal current, a la Laframboise
    """
    if isinstance(V, (int, float)):
        V = np.array([V], dtype=np.float)
    elif isinstance(V, (list, tuple)):
        V = np.array(V, dtype=np.float)

    I = np.zeros_like(V)

    q, m, n, T = species.q, species.m, species.n, species.T
    kappa, alpha = species.kappa, species.alpha
    k = constants('Boltzmann constant')

    eta = q*V/(k*T)        # Normalized voltage
    vth = np.sqrt(k*T/m)   # Thermal velocity

    indices_n = np.where(eta > 0)   # indices for repelled particles
    indices_p = np.where(eta <= 0)  # indices for attracted particles

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

    if geometry.shape == 'sphere':
        r = geometry.r
        I0 = 2*np.sqrt(2*np.pi)*r**2*q*n*vth # Current due to random particle flux for a spherical probe

        # repelled particles:
        if species.dist == 'maxwellian' or species.dist == 'cairns':
            I[indices_n] = I0 * C * D * np.exp(-eta[indices_n]) * (1. + E * eta[indices_n] * (eta[indices_n] + 4.))

        elif species.dist == 'kappa' or species.dist == 'kappa-cairns':
            I[indices_n] = I0 * C * D * (1. + eta[indices_n] / (kappa - 1.5))**(
                1. - kappa) * (1. + E * eta[indices_n] * (eta[indices_n] + 4. * ((kappa - 1.5) / (kappa - 1.))))
        
        # attracted particles:
        eta[indices_p] = np.abs(eta[indices_p])
        I[indices_p] = I0*C*D*(1.+F*eta[indices_p])

        if normalize:
            I /= I0

    elif geometry.shape == 'cylinder':
        r, l = geometry.r, geometry.l
        I0 = np.sqrt(2*np.pi)*r*l*q*n*vth # Current due to random particle flux for a cylindrical probe

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
            
        if normalize:
            I /= I0 
     
    elif geometry.shape == 'plane':
        raise ValueError('Geometry {} not implemented yet'.format(geometry.shape))

    else:
        raise ValueError('Geometry {} not supported'.format(geometry.shape))
    
    return I 

def current(geometry, species, V):

    if isinstance(species, list):
        I = 0
        for s in species:
            I += current(geometry, species, V)
        return I

    """
    if R=0:
        call OML function for correct distribution
        It should take care of geometry, attracted/repelled
    elif R=inf:
        thin sheath
    else:
        table
    """


def lafr_norm_current(geometry, R, n, T, q=None, m=None):

    if q==None:
        q = -constants('elementary charge')

    if m==None:
        m = constants('electron mass')

    k = constants('Boltzmann constant')

    if geometry.lower()=='sphere':
        I0 = q*n*R**2*np.sqrt(8*np.pi*k*T/m)
    elif geometry.lower()=='cylinder':
        I0 = q*n*R*np.sqrt(2*np.pi*k*T/m)
    else:
        raise ValueError('Geometry {} not supported'.format(geometry))

    return I0

def lafr_attr_current(geometry, kind='linear'):

    if geometry.lower()=='sphere':
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
              [1.000, 1.094, 1.245, 1.402, 1.534, 1.632, 1.694, 1.762, 1.833, 1.891,  1.938,  2.022,  2.097,  2.166]]   # R=100

    elif geometry.lower()=='cylinder':
        Rs =  [0, 1, 1.5, 2, 2.5, 3, 4, 5, 10, 20, 30, 40, 50, 100]

        Vs =  [0.0  , 0.1   , 0.3   , 0.6   , 1.0   , 1.5   , 2.0   , 3.0   , 5.0   , 7.5   , 10.0  , 15.0  , 20.0  , 25.0]
        Is = [[1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.7551, 1.9320, 2.2417, 2.7555, 3.2846, 3.7388, 4.5114, 5.1695, 5.7526],   # R=0
              [1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.7551, 1.9320, 2.2417, 2.7555, 3.2846, 3.7388, 4.5114, 5.1695, 5.7525],   # R=1
              [1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.7551, 1.9320, 2.2417, 2.7555, 3.2846, 3.735 , 4.493 , 5.141 , 5.711 ],   # R=1.5
              [1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.7551, 1.9320, 2.247 , 2.750 , 3.266 , 3.703 , 4.439 , 5.060 , 5.607 ],   # R=2
              [1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.7551, 1.9320, 2.237 , 2.731 , 3.227 , 3.645 , 4.342 , 4.936 , 5.462 ],   # R=2.5
              [1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.754 , 1.928 , 2.226 , 2.701 , 3.174 , 3.567 , 4.235 , 4.789 , 5.291 ],   # R=3
              [1.000, 1.0804, 1.2101, 1.3721, 1.554 , 1.747 , 1.913 , 2.192 , 2.626 , 3.050 , 3.402 , 3.990 , 4.489 , 4.926 ],   # R=4
              [1.000, 1.0804, 1.2100, 1.371 , 1.549 , 1.735 , 1.893 , 2.151 , 2.544 , 2.920 , 3.231 , 3.749 , 4.183 , 4.565 ],   # R=5
              [1.000, 1.0803, 1.208 , 1.362 , 1.523 , 1.677 , 1.798 , 1.98  , 2.22  , 2.442 , 2.622 , 2.919 , 3.166 , 3.384 ],   # R=10
              [1.000, 1.0803, 1.205 , 1.348 , 1.486 , 1.605 , 1.689 , 1.801 , 1.940 , 2.060 , 2.157 , 2.319 , 2.455 , 2.576 ],   # R=20
              [1.000, 1.0803, 1.205 , 1.348 , 1.486 , 1.605 , 1.689 , 1.801 , 1.940 , 2.060 , 2.157 , 2.082 , 2.177 , 2.262 ],   # R=30
              [1.000, 1.0803, 1.205 , 1.348 , 1.486 , 1.605 , 1.689 , 1.801 , 1.940 , 2.060 , 2.157 , 2.082 , 2.025 , 2.092 ],   # R=40
              [1.000, 1.0803, 1.198 , 1.327 , 1.439 , 1.523 , 1.576 , 1.638 , 1.703 , 1.756 , 1.798 , 1.868 , 1.929 , 1.983 ],   # R=50
              [1.000, 1.0803, 1.194 , 1.314 , 1.409 , 1.478 , 1.518 , 1.561 , 1.599 , 1.628 , 1.650 , 1.686 , 1.719 , 1.747 ]]   # R=100

    else:
        raise ValueError('Geometry {} not supported'.format(geometry))

    f = interp2d(Vs, Rs, Is, kind=kind)
    return lambda R, V: float(f(V, R))
