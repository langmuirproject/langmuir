Langmuir
========

.. image:: https://travis-ci.com/sigvaldm/langmuir.svg?branch=master
    :target: https://travis-ci.com/sigvaldm/langmuir

.. image:: https://coveralls.io/repos/github/sigvaldm/langmuir/badge.svg?branch=master
    :target: https://coveralls.io/github/sigvaldm/langmuir?branch=master

.. image:: https://img.shields.io/pypi/pyversions/langmuir.svg
    :target: https://pypi.org/project/langmuir

Functions for work on Langmuir probe theory

Installation
------------
Install from PyPI using ``pip`` (preferred method)::

    pip install langmuir

Or download the GitHub repository https://github.com/sigvaldm/langmuir.git and run::

    python setup.py install

Usage
-----

The tables for attracted-species current for finite-radius probes in an isothermal Maxwellian plasma given by Laframboise is implemented. E.g. to get the normalized current for a spherical probe of 1 Debye length and a normalized potential of 25::

    >>> from langmuir import *
    >>> R = 1
    >>> eV_kT = 25

    >>> I = lafr_attr_current('Sphere')
    >>> I(R, eV_kT)
    21.895

The function interpolates between values given in Laframboise's tables.
To get the current in Ampére's you must find the normalizing current::

    >>> n=1e11
    >>> T=1e3

    >>> I0 = lafr_norm_current('Sphere', R, n, T)
    >>> I0*I(R, eV_kT)
    -0.21650719872149787

Likewise for cylindrical probes. The current is then in Ampère's per meter so
you must multiply by the probe length::

    >>> l = 25e-3
    >>> I = lafr_attr_current('Cylinder')
    >>> I0 = lafr_norm_current('Cylinder', R, n, T)
    >>> I0*l*I(R, eV_kT)
    -0.0007110399980848462
