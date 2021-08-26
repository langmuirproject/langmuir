Specifying the plasma
=====================
A species is specified using a range of keyword arguments upon initialization of a ``Species`` object, for instance::

    Species(Z=1, amu=16, n=1e11, T=1000)

This is a singly charged ion (:math:`Z=1`) with a mass of :math:`16\,\mathrm{AMU}` (atomic mass units), which corresponds to the most abundant oxygen isotope.
The density is :math:`n=10^{11}\,\mathrm{m^{-3}}`, and it is Maxwellian distributed with a temperature of :math:`T=1000\,\mathrm{K}` and no drift-velocity.
The full set of keywords is given below:

+---------------------------+-----------+--------------------------+-----------------------------------+
| Quantity                  | Keyword   | Units                    | Default value                     |
+===========================+===========+==========================+===================================+
| Charge                    | ``q``     | Coulombs                 | elementary                        |
|                           +-----------+--------------------------+ charge                            |
|                           | ``Z``     | elementary charges       |                                   |
+---------------------------+-----------+--------------------------+-----------------------------------+
| Mass                      | ``m``     | kilograms                | electron                          |
|                           +-----------+--------------------------+ mass                              |
|                           | ``amu``   | atomic mass units        |                                   |
+---------------------------+-----------+--------------------------+-----------------------------------+
| Density                   | ``n``     | partices per cubic meter | :math:`10^{11}\,\mathrm{m^{-3}}`  |
+---------------------------+-----------+--------------------------+-----------------------------------+
| Temperature/thermal speed | ``T``     | Kelvin                   | :math:`1000\,\mathrm{K}`          |
|                           +-----------+--------------------------+                                   |
|                           | ``eV``    | electron-volts           |                                   |
|                           +-----------+--------------------------+                                   |
|                           | ``vth``   | meters per second        |                                   |
+---------------------------+-----------+--------------------------+-----------------------------------+
| Spectral index kappa      | ``kappa`` | dimensionless            | :math:`\infty` (``float('inf')``) |
+---------------------------+-----------+--------------------------+-----------------------------------+
| Spectral index alpha      | ``alpha`` | dimensionless            | :math:`0`                         |
+---------------------------+-----------+--------------------------+-----------------------------------+

Kappa, Cairns or Kappa-Cairns distributed particles are obtained by specifying ``kappa``, ``alpha``, or both, respectively. Maxwell, Kappa and Cairns distributions are all limiting cases of the Kappa-Cairns distribution [Darian]_.

For convenience, the following subclasses of ``Species`` exist:

- ``Electron``
- ``Proton``
- ``Positron``
- ``Antiproton``
- ``Vogon``
- ``Hydrogen``
- ``Helium``
- ``Lithium``
- and similar for the rest of the periodic table

The only difference from ``Species`` is that these have different default charge and mass to represent the named element. ``Oxygen`` for instance, defaults to singly charged oxygen of the most abundant isotope. Hence the example in the beginning of this section could also have been written::

    Oxygen(n=1e11, T=1000)

Finally, a multi-species plasma is represented as a list of its constituents. A typical Oxygen plasma would for instance be::

    plasma = [
        Electron(n=1e11, T=1000),
        Oxygen(n=1e11, T=1000)
        ]

Quick computations of plasma parameters
---------------------------------------
The Langmuir library is also suitable for doing quick computations of fundamental plasma parameters, since ``Species`` computes these upon initialization. One example is calculating the electron Debye length of a certain plasma::

    >>> Electron(n=1e11, T=1000).debye
    0.00690089806774598

The ``Species`` class defines the following useful members:

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

The latter four members are methods which take the magnitude of the magnetic flux density as an argument. In addition, every valid keyword argument of the constructor is also a valid member. This may conveniently be used for instance to convert a temperature from electron-volts to Kelvin::

    >>> Species(eV=0.2).T
    2320.9036243100163

In this case we only specified the input ``eV``, since we know that temperature do not depend on density.

Finally, the total Debye length of a plasma consisting of multiple species can be obtained using the ``debye()`` function. For the oxygen plasma mentioned previously::

    >>> debye(plasma)
    0.004879671013271479

Note that Langmuir uses the correct expression for the Debye length also for general Kappa-Cairns distributed plasmas [Darian]_.
