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
The Langmuir library contains a collection of functions that compute the current collected by a conductor immersed in a plasma according to various models (characteristics). These functions take as arguments the probe geometry, for instance ``Sphere``, and a plasma, where a plasma may consist of one or more ``Species``.

As an example, consider a spherical Langmuir probe of radius :math:`r=1\,\mathrm{mm}` immersed in a plasma with an electron density of :math:`n=10^{11}\,\mathrm{m^{-3}}` and an electron temperature of :math:`T=1000\,\mathrm{K}`. The electron current collected by this probe when it has a voltage of :math:`V=4.5\,\mathrm{V}` with respect to the background plasma is computed according to *Orbital Motion-Limited* (OML) theory as follows::

    >>> OML_current(Sphere(r=1e-3), Electron(n=1e11, T=1000), V=4.5)
    -5.262656728335636e-07

Let's consider a more complete example. Below we compare the current-voltage characteristics predicted by the OML theory and the *finite-length* (FL) model for a cylindrical probe with an ideal guard on one end.

.. literalinclude:: demo/basic.py
.. image:: docs/source/basic.png

This example demonstrates that accounting for edge effects on a probe of finite length leads to a larger collected current. Also note that the characteristic correctly captures both the electron saturation, the electron retardation, and the ion saturation regions. Beware that the current collected by Langmuir probes (thus going into it) is usually negative. It is common practice, however, to invert it prior to plotting.

More accurate characteristics such as those available in Langmuir allows the study of non-ideal effects, as well as more accurate inference techniques of plasma parameters. See for instance :doc:`examples`.
