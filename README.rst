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
  This model interpolates between tabulated values by Laframboise and Darian et al. which takes into account the effect of finite radius. If radius, voltage or spectral indices are outside the convex hull of the tabulated values, ``nan`` are returned. Valid radii are 0-10 debye lengths (100 for Maxwellian), normalized voltages must be between -25 and 0 (repelled species are neglected), alpha between 0 and 0.2, and kappa no less than 4.
- ``thermal_current(geometry, species, normalize=False)``
  This is the current absorbed by an object fixed at zero potential due to random thermal motion of particles.
- ``normalization_current(geometry, species)``
  This is the current which is used for normalizing the current from the other models. For Maxwellian this coincides with the thermal current, whereas for other distributions it is what would have been the thermal current if the particles were Maxwellian.

As an example, the following snippet computes the normalized electron current of a probe of 3 Debye lengths radius and normalized voltage of -10::

    >>> sp  = Species(n=1e11, T=1000)
    >>> geo = Cylinder(r=3*sp.debye, l=1)
    >>> I = finite_radius_current(geo, sp, eta=-10, normalize=True)

Notice that setting ``l==1`` means you get the current per unit length.

Inverse problems
----------------
Sometimes the collected current of one or more probes is known and one would like to solve for one or more other parameters. The Langmuir library do not address this analytically in part due to the vast number of such inverse problems, and in part due to some characteristics not being invertible (for instance those who are of tabulated values). However, it is in principle possible to apply numerical methods of root solving, least squares, etc. along with the models in Langmuir.

Consider a cylindrical probe with known dimensions and a positive but unknown voltage collecting a current of -0.4uA in a Maxwellian plasma with known density and temperature. What is the voltage? We shall neglect the current due to ions, and define a residual function. This residual is the difference between the current collected by a probe at a given potential, and the actual collected current, and it is used by a least squares algorithm to compute the voltage::

    >>> from langmuir import *
    >>> from scipy.optimize import leastsq

    >>> sp = Species(n=1e11, T=1000)
    >>> geo = Cylinder(1e-3, 25e-3)
    >>> I = -0.4e-6

    >>> def residual(V):
    >>>     return finite_radius_current(geo, sp, V) - I

    >>> x, c = leastsq(residual, 0)
    >>> print(x[0])
    0.6265540484991013

The reader may verify that this voltage indeed results in the correct current. Notice also that we were in fact able to invert the model ``finite_radius_current``, which consists of tabulated values.

A slightly more interesting inversion problem, is that of determining the ionospheric density from four cylindrical Langmuir probes with known bias voltages with respect to a spacecraft, but an unknown floating potential ``V0`` of the spacecraft with respect to the plasma. We shall assume the bias voltages to be 2.5, 4.0, 5.5 and 7.0 volts. In the below example, we first construct the currents for such probes by assuming a floating potential and a set of plasma parameters, but we do not use this knowledge in the inversion. We do, however, make an initial guess ``x0`` which we believe are somewhat close to the answer::

    >>> from langmuir import *
    >>> from scipy.optimize import leastsq

    >>> geo = Cylinder(1e-3, 25e-3)
    >>> V0 = -0.5
    >>> V = np.array([2.5, 4.0, 5.5, 7.0])
    >>> I = OML_current(geo, Species(n=120e10, T=1000), V+V0)

    >>> def residual(x):
    >>>     n, V0 = x
    >>>     return OML_current(geo, Species(n=n, T=1500), V+V0) - I

    >>> x0 = [10e10, -0.3]
    >>> x, c = leastsq(residual, x0)
    >>> n, V0 = x

    >>> print(n)
    1199899818493.931

    >>> print(V0)
    -0.5417515655165968

The method correctly determined the density to be 120e10. However, the floating potential ``V0`` is off by almost ten percent. The reason is that the temperature is considered unknown, and assumed to be 1500K when solving the problem, while it is actually 1000K. Since we have four measurements (four equations) and only two unknowns, it is tempting to also include the temperature as an unknown parameter and try to solve for it. However, if this is done the least squares algorithm will fail miserably. The reason is that the set of equations arising for the attracted-species current of cylindrical probes are singular and cannot be solved for even analytically. Fortunately, both the temperature and floating potential can be eliminated from the equation when analytically solving for the density, and similarly it also works to obtain the density from the least squares algorithm. Since the floating potential and temperature represent a coupled unknown which cannot be solved for, an error in assuming one is reflected as an error in the other.

This demonstrates the usefulness as well as challenges and subtleties of solving inverse Langmuir problems.
