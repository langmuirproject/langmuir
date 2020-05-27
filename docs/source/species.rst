Specifying the plasma
=====================
A species is specified using a range of keyword arguments upon initialization of a ``Species`` object, for instance::

    Species(Z=1, amu=16, n=1e11, T=1000)

This is a singly charged ion (:math:`Z=1`) with a mass of :math:`16\,\mathrm{AMU}` (atomic mass units), which corresponds to the most abundant oxygen isotope.
The density is :math:`n=10^{11}\,\mathrm{m^{-3}}`, and it is Maxwellian distributed with a temperature of :math:`T=1000\,\mathrm{K}` and no drift-velocity.
The full set of keywords is given below:

+---------------------------+-----------+--------------------------+----------------------------------+
| Quantity                  | Keyword   | Units                    | Default value                    |
+===========================+===========+==========================+==================================+
| Charge                    | ``q``     | Coulombs                 | electron                         |
|                           +-----------+--------------------------+ charge                           |
|                           | ``Z``     | elementary charges       |                                  |
+---------------------------+-----------+--------------------------+----------------------------------+
| Mass                      | ``m``     | kilograms                | electron                         |
|                           +-----------+--------------------------+ mass                             |
|                           | ``amu``   | atomic mass units        |                                  |
+---------------------------+-----------+--------------------------+----------------------------------+
| Density                   | ``n``     | partices per cubic meter | :math:`10^{11}\,\mathrm{m^{-3}}` |
+---------------------------+-----------+--------------------------+----------------------------------+
| Temperature/thermal speed | ``T``     | Kelvin                   | :math:`1000\,\mathrm{K}`         |
|                           +-----------+--------------------------+                                  |
|                           | ``eV``    | electron-volts           |                                  |
|                           +-----------+--------------------------+                                  |
|                           | ``vth``   | meters per second        |                                  |
+---------------------------+-----------+--------------------------+----------------------------------+
| Spectral index kappa      | ``kappa`` | dimensionless            | ``float('inf')``                 |
+---------------------------+-----------+--------------------------+----------------------------------+
| Spectral index alpha      | ``alpha`` | dimensionless            | ``0``                            |
+---------------------------+-----------+--------------------------+----------------------------------+

Kappa, Cairns or Kappa-Cairns distributed particles are obtained by specifying ``kappa``, ``alpha``, or both, respectively.

For convenience the following subclasses of ``Species`` exist, the only difference being that the default charge and mass is according to the particle type the class is named after:

- ``Electron``
- ``Proton``
- ``Positron``
- ``Antiproton``
