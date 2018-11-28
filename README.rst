Langmuir
========

.. image:: https://travis-ci.com/langmuirproject/langmuir.svg?branch=master
    :target: https://travis-ci.com/langmuirproject/langmuir

.. image:: https://coveralls.io/repos/github/langmuirproject/langmuir/badge.svg?branch=master
    :target: https://coveralls.io/github/langmuirproject/langmuir?branch=master

.. image:: https://img.shields.io/pypi/pyversions/langmuir.svg
    :target: https://pypi.org/project/langmuir

Programmatically accessible current-voltage characteristics for ideal and non-ideal Langmuir probes.

Installation
------------
Install from PyPI using ``pip`` (preferred method)::

    pip install langmuir

Or download the GitHub repository https://github.com/langmuirproject/langmuir.git and run::

    python setup.py install

Getting Started
---------------
The Langmuir library contains a collection of functions that compute the current collected by a Langmuir probe according to various theories and models. These functions take as arguments the probe geometry, for instance ``Cylinder``, and a plasma, where a plasma may consist of one or more ``Species``.

As an example, consider a 25mm long cylindrical probe with radius 0.255mm. The plasma consists of electrons and singly charged oxygen ions, both with a density of 1e11 and a temperature of 1000K. The current-voltage charactersitic according to OML theory is easily computed and plotted::

    >>> from langmuir import *
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt

    >>> n = 1e11
    >>> T = 1000
    >>> r = 0.255e-3
    >>> l = 25e-3

    >>> plasma = []
    >>> plasma.append(Species('electron' , n=n, T=T))
    >>> plasma.append(Species(amu=16, Z=1, n=n, T=T))

    >>> Vs = np.linspace(-2, 2, 100)
    >>> Is = OML_current(Cylinder(r, l), plasma, Vs)

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)

    >>> ax.plot(Vs, -Is*1e6)
    >>> ax.set_xlabel(r'$V_{\mathrm{p}} [\mathrm{V}]$')
    >>> ax.set_ylabel(r'$I_{\mathrm{p}} [\mathrm{\mu A}]$')
    >>> ax.grid(True)
    >>> plt.show()

.. image:: IV_characteristics.png

Notice that the characteristic includes all regions (ion saturation, electron retardation and electron saturation), and do not rely on approximations to the OML theory requiring the voltage to be within a certain range. What's more, it's easy to take into account for instance finite-radius effects, by replacing ``OML_current()`` with ``finite_radius_current()``.

Specifying the geometry
-----------------------

Specifying the plasma
---------------------

Computing fundamental plasma parameters
---------------------------------------

Models for collected current
----------------------------
- Include example of normalizing

Solving for an unknown voltage
------------------------------

DEPRECATED: Usage of Tables
---------------------------

The tables for attracted-species current for finite-radius probes in an isothermal Maxwellian plasma given by Laframboise is implemented. E.g. to get the normalized current for a spherical probe of 1 Debye length and a normalized potential of 25::

    >> from langmuir import *
    >> R = 1
    >> eV_kT = 25

    >> f = lafr_attr_current('Sphere')
    >> I = f(R, eV_kT)
    >> print("{:.3f}".format(I))
    21.895

The function linearly interpolates between values given in Laframboise's tables.
The argument ``kind`` can be used to change to quadratic interpolation.
To get the current in Ampére's you must find the normalizing current::

    >> n=1e11
    >> T=1e3

    >> I0 = lafr_norm_current('Sphere', R, n, T)
    >> I = I0*f(R, eV_kT)
    >> print("{:.1f}mA".format(I*1e3))
    -216.5mA

Likewise for cylindrical probes. The current is then in Ampère's per meter so
you must multiply by the probe length::

    >> l = 25e-3
    >> f = lafr_attr_current('Cylinder')
    >> I0 = lafr_norm_current('Cylinder', R, n, T)
    >> I = I0*l*f(R, eV_kT)
    >> print("{:.1f}uA".format(I*1e6))
    -711.0uA
