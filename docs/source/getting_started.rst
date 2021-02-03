Getting Started
===============
The Langmuir library contains a collection of functions that compute the current collected by a conductor immersed in a plasma according to various models. These functions take as arguments the probe geometry, for instance ``Sphere``, and a plasma, where a plasma may consist of one or more ``Species``.

As an example, consider a spherical Langmuir probe of radius :math:`r=1\,\mathrm{mm}` immersed in a plasma with an electron density of :math:`n=10^{11}\,\mathrm{m^{-3}}` and an electron temperature of :math:`T=1000\,\mathrm{K}`. The electron current collected by this probe when it has a voltage of :math:`V=4.5\,\mathrm{V}` with respect to the background plasma is computed according to *Orbital Motion-Limited* (OML) theory as follows::

    >>> OML_current(Sphere(r=1e-3), Electron(n=1e11, T=1000), V=4.5)
    -5.262656728335636e-07

Let's consider a more complete example. Below we compare the current-voltage characteristics predicted by the OML theory and the *finite-length* (FL) model for a cylindrical probe with an ideal guard on one end.

.. literalinclude:: ../../demo/basic.py
.. image:: basic.png

This example demonstrates that accounting for edge effects on a probe of finite length leads to a larger collected current. Also note that the characteristic correctly captures both the electron saturation, the electron retardation, and the ion saturation regions.
