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

from scipy.constants import value as constants
from scipy.constants import epsilon_0
from langmuir.misc import *
import numpy as np

class Species(object):
    """
    Defines a species using a set of flags and keyword parameters. Maxwellian
    electrons are the default but Kappa, Cairns, and Kappa-Cairns distributions
    can be specified by setting the spectral indices kappa and alpha. Other
    parameters can be specified/overridden using the keyword parameters.
    Examples::

        >>> # Maxwellian electrons with density 1e11 m^(-3) and 1000 K
        >>> species = Electron(n=1e11, T=1000)

        >>> # Kappa-distributed electrons with kappa=2
        >>> species = Electron(n=1e11, T=1000, kappa=2)

        >>> # Doubly charged Maxwellain Oxygen ions with 0.26 eV
        >>> species = Species(Z=2, amu=16, n=1e11, eV=0.26)

        >>> # Kappa-Cairns-distributed protons with kappa = 5 and alpha=0.2
        >>> species = Proton(n=1e11, T=1000, kappa=5, alpha=0.2)

    A plasma is fully specified by a list of species. E.g. for a Maxwellian
    electron-proton plasma::

        >>> plasma = []
        >>> plasma.append(Electron(n=1e11, T=1000))
        >>> plasma.append(Proton(n=1e11, T=1000))

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
    kappa : spectral index kappa for Kappa and Kappa-Cairns distribution
    alpha : spectral index alpha for Cairns and Kappa-Cairns distribution
    """

    def_Z   = -1
    def_amu = constants('electron mass')/constants('atomic mass constant')

    def __init__(self, **kwargs):

        e    = constants('elementary charge')
        k    = constants('Boltzmann constant')
        eps0 = epsilon_0 #constants('vacuum electric permittivity')
        amu  = constants('atomic mass constant')

        if 'Z' in kwargs:
            self.q = kwargs['Z']*e
        elif 'q' in kwargs:
            self.q = kwargs['q']
        else:
            self.q = self.def_Z*e

        if 'amu' in kwargs:
            self.m = kwargs['amu']*amu
        elif 'm' in kwargs:
            self.m = kwargs['m']
        else:
            self.m = self.def_amu*amu

        if 'eV' in kwargs:
            self.T = kwargs['eV']*e/k
        elif 'vth' in kwargs:
            self.T = self.m*kwargs['vth']**2/k
        elif 'T' in kwargs:
            self.T = kwargs['T']
        else:
            self.T = 1000

        if self.T<0:
            self.T=0
            logger.warning('Negative temperature interpreted as zero')

        self.n     = kwargs.pop('n', 1e11)
        self.alpha = kwargs.pop('alpha', 0)
        self.kappa = kwargs.pop('kappa', float('inf'))

        self.Z = self.q/e
        self.amu = self.m/amu
        self.eV  = self.T*k/e
        self.vth = np.sqrt(k*self.T/self.m)
        self.omega_p = np.sqrt(self.q**2*self.n/(eps0*self.m))
        self.freq_p = self.omega_p/(2*np.pi)
        self.period_p = 1/self.freq_p

        if self.kappa == float('inf'):
            self.debye = np.sqrt(eps0*k*self.T/(self.q**2*self.n))*np.sqrt((1.0 + 15.0*self.alpha)/(1.0 + 3.0*self.alpha))
        else:
            self.debye = np.sqrt(eps0 * k * self.T / (self.q**2 * self.n)) *\
                         np.sqrt(((self.kappa - 1.5) / (self.kappa - 0.5)) *\
                        ((1.0 + 15.0 * self.alpha * ((self.kappa - 1.5) / (self.kappa - 2.5))) / (1.0 + 3.0 * self.alpha * ((self.kappa - 1.5) / (self.kappa - 0.5)))))

    def omega_c(self, B):
        """
        Returns the angular cyclotron frequency of the species for a given B-field
        """
        return np.abs(B*self.q/self.m)

    def freq_c(self, B):
        """
        Returns the cyclotron frequency of the species for a given B-field
        """
        return self.omega_c(B)/(2*np.pi)

    def period_c(self, B):
        """
        Returns the cyclotron period of the species for a given B-field
        """
        return 1/self.freq_c(B)

    def larmor(self, B):
        """
        Returns the Larmor radius of the species for a given B-field
        """
        return self.vth/self.omega_c(B)

    def __repr__(self):

        s = "Species(q={}, m={}, n={}, T={}".format(self.q, self.m, self.n, self.T)

        if(self.kappa != float('inf')):
            s += ", kappa={}".format(self.kappa)

        if(self.alpha != 0):
            s += ", alpha={}".format(self.alpha)

        s += ")"

        return s

def debye(species):
    """
    Get the total Debye length of a species or a plasma (represented as a list
    of species). If species is a list this includes the Debye length due to
    all species.
    """
    if not isinstance(species, list):
        species = [species]
    return sum([s.debye**(-2) for s in species])**(-0.5)

class Electron(Species): pass
class Proton(Species): def_Z = 1; def_amu = 1
class Positron(Species): def_Z = 1
class Antiproton(Species): def_amu = 1
class Vogon(Species): def_Z = 1j
