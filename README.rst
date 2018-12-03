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

Specifying a species
--------------------
A species is specified using a range of keyword arguments upon initialization of a ``Species`` object.
At the very least, you must specify the density and the temperature (or thermal speed)::

    >>> Species(n=1e11, T=1000)

Since nothing more is specified, this represents Maxwellian-distributed electrons. The full set of keywords is as follows:

+---------------------------+-----------+--------------------------+------------------+
| Quantity                  | Keyword   | Units                    | Default value    |
+===========================+===========+==========================+==================+
| Charge                    | ``q``     | Coulombs                 | electron         |
|                           +-----------+--------------------------+ charge           |
|                           | ``Z``     | elementary charges       |                  |
+---------------------------+-----------+--------------------------+------------------+
| Mass                      | ``m``     | kilograms                | electron         |
|                           +-----------+--------------------------+ mass             |
|                           | ``amu``   | atomic mass units        |                  |
+---------------------------+-----------+--------------------------+------------------+
| Density                   | ``n``     | partices per cubic meter |                  |
+---------------------------+-----------+--------------------------+------------------+
| Temperature/thermal speed | ``T``     | Kelvin                   |                  |
|                           +-----------+--------------------------+                  |
|                           | ``eV``    | electron-volts           |                  |
|                           +-----------+--------------------------+                  |
|                           | ``vth``   | meters per second        |                  |
+---------------------------+-----------+--------------------------+------------------+
| Spectral index kappa      | ``kappa`` | dimensionless            | ``float('inf')`` |
+---------------------------+-----------+--------------------------+------------------+
| Spectral index alpha      | ``alpha`` | dimensionless            | ``0``            |
+---------------------------+-----------+--------------------------+------------------+

The Langmuir library supports the general Kappa-Cairns velocity distribution, but defaults to Maxwellian, which is a special case of the Kappa-Cairns distribution.Kappa, Cairns or Kappa-Cairns distributed particles are obtained by specifying ``kappa``, ``alpha`` or both, respectively.

The charge and mass can also be specified using one of three self-descriptive flag:

- ``'electron'``
- ``'positron'``
- ``'proton'``

This must come before the keyword arguments, for instance::

    >>> Species('proton', n=1e11, eV=0.1)

Computing fundamental plasma parameters
---------------------------------------
The Langmuir library is also convenient to use for quick computations of fundamental plasma parameters, since ``Species`` computes these upon initialization. For instance, to get the electron Debye length of a plasma with a certain density and temperature::

    >>> Species(n=1e11, eV=0.1).debye
    0.007433942027347403

The following member variables and methods are accessible in ``Species``:

+-----------------+---------------------------------+
| Member          | Description                     |
+=================+=================================+
| ``debye``       | The Debye length                |
+-----------------+---------------------------------+
| ``omega_p``     | The angular plasma frequency    |
+-----------------+---------------------------------+
| ``freq_p``      | The linear plasma frequency     |
+-----------------+---------------------------------+
| ``period_p``    | The plasma period               |
+-----------------+---------------------------------+
| ``omega_c(B)``  | The angular cyclotron frequency |
+-----------------+---------------------------------+
| ``freq_c(B)``   | The linear cyclotron frequency  |
+-----------------+---------------------------------+
| ``period_c(B)`` | The cyclotron period            |
+-----------------+---------------------------------+
| ``larmor(B)``   | The larmor radius               |
+-----------------+---------------------------------+

The latter four members are methods which take the magnitude of the magnetic flux density as an argument. In addition, every valid keyword argument is also a valid member variable::

    >>> Species(n=1e11, T=1000).eV
    0.08617330337217212

Finally, the total Debye length of a plasma consisting of multiple species can be obtained using the ``debye()`` function::

    >>> plasma = []
    >>> plasma.append(Species('electron' , n=1e11, T=1000))
    >>> plasma.append(Species(amu=16, Z=1, n=1e11, T=1000))
    >>> debye(plasma)
    0.004879671013271479

Specifying the geometry
-----------------------
Langmuir supports two probe geometries, with self-descriptive names and the following signatures:

- ``Sphere(r)``
- ``Cylinder(r, l)``

``r`` and ``l`` representes the radius and length, respectively, of the geometry.

Models for collected current
----------------------------
Langmuir comes with several models for the collected current. Each model is represented by a function which takes a ``geometry`` and a ``species`` argument. The ``geometry`` is one of the above probe geometries, and the ``species`` parameters is either a single ``Species`` object or a list of such if it is desirable to take into account the effect of all species in a plasma. Most models also depend on the potential of the probe with respect to the background plasma. The potential can either be specified in volts using th ``V`` argument, or in terms of normalized voltage e*V/(k*T) using the ``eta`` argument. If these are Numpy arrays, the output will be a Numpy array of collected currents. The ``normalize`` argument can be set to ``True`` to return currents in normalized quantities rather than Ampére. Below is a description of all models:

- ``OML_current(geometry, species, V=None, eta=None, normalize=False)``
  Current according to the "Orbital Motion Limited" theory. This assumes an infinitesimal probe radius compared to the Debye length, and for cylinders, infinite length. Spheres with radii less than 0.2 Debye lengths, or Cylinders with radii less than 1 Debye length is usually very well approximated as having an infinitesimal radius.
- ``finite_radius_current(geometry, species, V=None, eta=None, table='laframboise-darian-marholm', normalize=False)``
  This model interpolates between tabulated values by Laframboise and Darian et al. which takes into account the effect of finite radius. If radius, voltage or spectral indices are outside the convex hull of the tabulated values, `nan` are returned. Valid radii are 0-10 debye lengths (100 for Maxwellian), normalized voltages must be between -25 and 0 (repelled species are neglected), alpha between 0 and 0.2, and kappa no less than 4.
- ``thermal_current(geometry, species, normalize=False)``
  This is the current absorbed by an object fixed at zero potential due to random thermal motion of particles.
- ``normalization_current(geometry, species)``
  This is the current which is used for normalizing the current from the other models. For Maxwellian this coincides with the thermal current, whereas for other distributions it is what would have been the thermal current if the particles were Maxwellian.

As an example, the following snippet computes the normalized electron current of a probe of 3 Debye lengths radius and normalized voltage of -10::

    >>> sp  = Species(n=1e11, T=1000)
    >>> geo = Cylinder(r=3*sp.debye, l=1)
    >>> I = finite_radius_current(geo, sp, eta=-10, normalize=True)

Notice that setting ``l==1`` means you get the current per unit length.

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
