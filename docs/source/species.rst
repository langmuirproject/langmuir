Specifying Species
==================
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
