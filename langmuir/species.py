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
                         np.sqrt(((self.kappa - 1.5) / (self.kappa - 0.5)) *\
                        ((1.0 + 15.0 * self.alpha * ((self.kappa - 1.5) / (self.kappa - 2.5))) / (1.0 + 3.0 * self.alpha * ((self.kappa - 1.5) / (self.kappa - 0.5)))))

    def __repr__(self):

        s = "Species('{}', q={}, m={}, n={}, T={}".format(self.dist, self.q, self.m, self.n, self.T)

        if(self.dist == 'kappa' or self.dist == 'kappa-cairns'):
            s += ", kappa={}".format(self.kappa)

        if(self.dist == 'cairns' or self.dist == 'kappa-cairns'):
            s += ", alpha={}".format(self.alpha)

        s += ")"

        return s
