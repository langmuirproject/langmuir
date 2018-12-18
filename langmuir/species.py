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
import numpy as np

class Species(object):
    """
    Defines a species using a set of flags and keyword parameters. Maxwellian
    electrons are the default but Kappa, Cairns, and Kappa-Cairns distributions
    can be specified by setting the spectral indices kappa and alpha. Other
    parameters can be specified/overridden using the keyword parameters.
    Examples::

        >>> # Maxwellian electrons with density 1e11 m^(-3) and 1000 K
        >>> species = Species(n=1e11, T=1000)

        >>> # Kappa-distributed electrons with kappa=2
        >>> species = Species(n=1e11, T=1000, kappa=2)

        >>> # Doubly charged Maxwellain Oxygen ions with 0.26 eV
        >>> species = Species(Z=2, amu=16, n=1e11, eV=0.26)

        >>> # Kappa-Cairns-distributed protons with kappa = 5 and alpha=0.2
        >>> species = Species('proton', n=1e11, T=1000, kappa=5, alpha=0.2)

    A plasma is fully specified by a list of species. E.g. for a Maxwellian
    electron-proton plasma::

        >>> plasma = []
        >>> plasma.append(Species('electron', n=1e11, T=1000))
        >>> plasma.append(Species('proton',   n=1e11, T=1000))

    Flags:
    ------
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
    kappa : spectral index kappa for Kappa and Kappa-Cairns distribution
    alpha : spectral index alpha for Cairns and Kappa-Cairns distribution
    """

    def __init__(self, *args, **kwargs):

        e    = constants('elementary charge')
        me   = constants('electron mass')
        mp   = constants('proton mass')
        k    = constants('Boltzmann constant')
        eps0 = constants('electric constant')
        amu  = constants('atomic mass constant')

        if 'Z' in kwargs:
            self.q = kwargs['Z']*e
        else:
            self.q = kwargs.pop('q', -e)

        if 'amu' in kwargs:
            self.m = kwargs['amu']*amu
        else:
            self.m = kwargs.pop('m', me)

        for arg in args:
            if arg.lower()=='proton':
                self.q = e
                self.m = mp

            if arg.lower()=='positron':
                self.q = e
                self.m = me

        self.Z = self.q/e
        self.amu = self.m/amu

        self.n = kwargs['n']

        if 'eV' in kwargs:
            self.T = kwargs['eV']*e/k
        elif 'vth' in kwargs:
            self.T = self.m*kwargs['vth']**2/k
        else:
            self.T = kwargs['T']

        self.vth = np.sqrt(k*self.T/self.m)
        self.eV  = self.T*k/e
        self.omega_p = np.sqrt(self.q**2*self.n/(eps0*self.m))
        self.freq_p = self.omega_p/(2*np.pi)
        self.period_p = 1/self.freq_p

        self.alpha = kwargs.pop('alpha', 0)
        self.kappa = kwargs.pop('kappa', float('inf'))

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
