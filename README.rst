Langmuir
========

.. image:: https://travis-ci.com/langmuirproject/langmuir.svg?branch=master
    :target: https://travis-ci.com/langmuirproject/langmuir

.. image:: https://coveralls.io/repos/github/langmuirproject/langmuir/badge.svg?branch=master
    :target: https://coveralls.io/github/langmuirproject/langmuir?branch=master

.. image:: https://img.shields.io/pypi/pyversions/langmuir.svg
    :target: https://pypi.org/project/langmuir

.. image:: https://zenodo.org/badge/149759145.svg
    :target: https://zenodo.org/badge/latestdoi/149759145

.. image:: https://readthedocs.org/projects/langmuir/badge/?version=latest
    :target: https://langmuir.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Programmatically accessible current-voltage characteristics for ideal and non-ideal Langmuir probes. See documentation on ReadTheDocs_.

.. _ReadTheDocs: http://langmuir.readthedocs.io

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

    >>> plasma = []
    >>> plasma.append(Electron(n=1e11, T=1000))
    >>> plasma.append(Proton(  n=1e11, T=1000))

    >>> Vs = np.linspace(-2, 8, 100)
    >>> Is = OML_current(Cylinder(r=1e-3, l=25e-3), plasma, Vs)

    >>> plt.plot(Vs, Is)
    >>> plt.show()

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)

    >>> ax.plot(Vs, -Is*1e6)
    >>> ax.set_xlabel(r'$V_{\mathrm{p}} [\mathrm{V}]$')
    >>> ax.set_ylabel(r'$I_{\mathrm{p}} [\mathrm{\mu A}]$')
    >>> ax.grid(True)
    >>> plt.show()

.. image:: IV_characteristics.png

Notice that the characteristic includes all regions (ion saturation, electron retardation and electron saturation), and do not rely on approximations to the OML theory requiring the voltage to be within a certain range. What's more, it's easy to take into account for instance finite-radius effects, by replacing ``OML_current()`` with ``finite_radius_current()``.
