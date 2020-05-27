Models for the collected current
================================

Langmuir comes with several models for the characteristic :math:`I(V)` of Langmuir probes. Each model is represented by a function which, in addition to the voltage, take a ``geometry`` and a ``species`` argument. ``geometry`` is one of the previously mentioned probe geometries (``Plane``, ``Cylinder`` or ``Sphere``), and the ``species`` argument is either a single ``Species`` object or a list of such can be used to account for multiple species. As in OML theory, it is assumed that the individual species' contribution to the current can be linearly superposed. In other words, Langmuir do not at the time account for a possible non-linear coupling between the species that may occur in the sheath or in a wake.

The voltage :math:`V` may either be specified directly through the argument ``V``, or the argument ``eta`` may be used instead. ``eta`` represents the  normalized voltage :math:`\eta=-qV/kT`, where :math:`q` and :math:`T` are the charge and temperature of the species, respectively, and :math:`k` is Boltzmann's constant. If ``V`` or ``eta`` is a Numpy array, the resulting current will also be a Numpy array. The ``normalization`` argument can be set to return normalized currents :math:`I/I_0` rather than currents in Ampére. The following normalization currents :math:`I_0` can be specified:

- ``'oml'``: The current according to OML theory. This is useful for comparing other models with OML theory.
- ``'th'``: The current collected due to the thermal motion of particles by a neutral probe, i.e., what OML theory predicts for :math:`\eta=0`. This is used for instance by Laframboise.
- ``'thmax'``: Same as ``th``, but this normalization current is for a Maxwellian plasma regardless of what the actual distribution specified in ``Species`` is.

In the case that ``species`` represents multiple species, it is unknown which species to normalize with respect to. Consequentially, the arguments ``V`` and ``normalization`` cannot be used in this case (it is of course possible to divide the result by :math:`I_0` for a chosen species manually).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   analytical
   finite_radius
   finite_length
