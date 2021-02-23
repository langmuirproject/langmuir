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

    def_Z   = 1
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

class Electron(Species): def_Z = -1
class Proton(Species): def_amu = 1
class Positron(Species): pass
class Antiproton(Species): def_Z = -1; def_amu = 1
class Vogon(Species): def_Z = 1j
class Hydrogen(Species): def_amu = 1.0080
class Helium(Species): def_amu = 4.00260
class Lithium(Species): def_amu = 7.0
class Beryllium(Species): def_amu = 9.012183
class Boron(Species): def_amu = 10.81
class Carbon(Species): def_amu = 12.011
class Nitrogen(Species): def_amu = 14.007
class Oxygen(Species): def_amu = 15.999
class Fluorine(Species): def_amu = 18.99840316
class Neon(Species): def_amu = 20.180
class Sodium(Species): def_amu = 22.9897693
class Magnesium(Species): def_amu = 24.305
class Aluminum(Species): def_amu = 26.981538
class Silicon(Species): def_amu = 28.085
class Phosphorus(Species): def_amu = 30.97376200
class Sulfur(Species): def_amu = 32.07
class Chlorine(Species): def_amu = 35.45
class Argon(Species): def_amu = 39.9
class Potassium(Species): def_amu = 39.098
class Calcium(Species): def_amu = 40.08
class Scandium(Species): def_amu = 44.95591
class Titanium(Species): def_amu = 47.87
class Vanadium(Species): def_amu = 50.941
class Chromium(Species): def_amu = 51.996
class Manganese(Species): def_amu = 54.93804
class Iron(Species): def_amu = 55.84
class Cobalt(Species): def_amu = 58.93319
class Nickel(Species): def_amu = 58.693
class Copper(Species): def_amu = 63.55
class Zinc(Species): def_amu = 65.4
class Gallium(Species): def_amu = 69.72
class Germanium(Species): def_amu = 72.63
class Arsenic(Species): def_amu = 74.92159
class Selenium(Species): def_amu = 78.97
class Bromine(Species): def_amu = 79.90
class Krypton(Species): def_amu = 83.80
class Rubidium(Species): def_amu = 85.468
class Strontium(Species): def_amu = 87.6
class Yttrium(Species): def_amu = 88.9058
class Zirconium(Species): def_amu = 91.22
class Niobium(Species): def_amu = 92.9064
class Molybdenum(Species): def_amu = 96.0
class Technetium(Species): def_amu = 97.90721
class Ruthenium(Species): def_amu = 101.1
class Rhodium(Species): def_amu = 102.9055
class Palladium(Species): def_amu = 106.4
class Silver(Species): def_amu = 107.868
class Cadmium(Species): def_amu = 112.41
class Indium(Species): def_amu = 114.82
class Tin(Species): def_amu = 118.71
class Antimony(Species): def_amu = 121.76
class Tellurium(Species): def_amu = 127.6
class Iodine(Species): def_amu = 126.9045
class Xenon(Species): def_amu = 131.29
class Cesium(Species): def_amu = 132.9054520
class Barium(Species): def_amu = 137.33
class Lanthanum(Species): def_amu = 138.9055
class Cerium(Species): def_amu = 140.12
class Praseodymium(Species): def_amu = 140.9077
class Neodymium(Species): def_amu = 144.24
class Promethium(Species): def_amu = 144.91276
class Samarium(Species): def_amu = 150.4
class Europium(Species): def_amu = 151.96
class Gadolinium(Species): def_amu = 157.2
class Terbium(Species): def_amu = 158.92535
class Dysprosium(Species): def_amu = 162.50
class Holmium(Species): def_amu = 164.93033
class Erbium(Species): def_amu = 167.26
class Thulium(Species): def_amu = 168.93422
class Ytterbium(Species): def_amu = 173.04
class Lutetium(Species): def_amu = 174.967
class Hafnium(Species): def_amu = 178.5
class Tantalum(Species): def_amu = 180.9479
class Tungsten(Species): def_amu = 183.8
class Rhenium(Species): def_amu = 186.21
class Osmium(Species): def_amu = 190.2
class Iridium(Species): def_amu = 192.22
class Platinum(Species): def_amu = 195.08
class Gold(Species): def_amu = 196.96657
class Mercury(Species): def_amu = 200.59
class Thallium(Species): def_amu = 204.383
class Lead(Species): def_amu = 207
class Bismuth(Species): def_amu = 208.9804
class Polonium(Species): def_amu = 208.98243
class Astatine(Species): def_amu = 209.98715
class Radon(Species): def_amu = 222.01758
class Francium(Species): def_amu = 223.01973
class Radium(Species): def_amu = 226.02541
class Actinium(Species): def_amu = 227.02775
class Thorium(Species): def_amu = 232.038
class Protactinium(Species): def_amu = 231.0359
class Uranium(Species): def_amu = 238.0289
class Neptunium(Species): def_amu = 237.04817
class Plutonium(Species): def_amu = 244.06420
class Americium(Species): def_amu = 243.06138
class Curium(Species): def_amu = 247.07035
class Berkelium(Species): def_amu = 247.07031
class Californium(Species): def_amu = 251.07959
class Einsteinium(Species): def_amu = 252.0830
class Fermium(Species): def_amu = 257.09511
class Mendelevium(Species): def_amu = 258.09843
class Nobelium(Species): def_amu = 259.10100
class Lawrencium(Species): def_amu = 262.110
class Rutherfordium(Species): def_amu = 267.122
class Dubnium(Species): def_amu = 268.126
class Seaborgium(Species): def_amu = 271.134
class Bohrium(Species): def_amu = 274.144
class Hassium(Species): def_amu = 277.152
class Meitnerium(Species): def_amu = 278.156
class Darmstadtium(Species): def_amu = 281.165
class Roentgenium(Species): def_amu = 282.169
class Copernicium(Species): def_amu = 285.177
class Nihonium(Species): def_amu = 286.183
class Flerovium(Species): def_amu = 289.191
class Moscovium(Species): def_amu = 290.196
class Livermorium(Species): def_amu = 293.205
class Tennessine(Species): def_amu = 294.211
class Oganesson(Species): def_amu = 294.214
# Extracted from PubChem
