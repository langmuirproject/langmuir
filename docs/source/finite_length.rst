The finite length model
=======================

- ``finite_length_current(geometry, species, V=None, eta=None, normalization=None)``
  The Marholm-Marchand model for finite-length probes. Works for normalized
  voltages up to 100.

- ``finite_length_current_density(geometry, species, V=None, eta=None, z=None, zeta=None, normalization=None)``
  The current per unit length according to the Marholm-Marchand model for
  finite-length probes. Works for normalized voltages up to 100. ``z`` is
  position on the probe, and ``zeta`` is position normalized by the Debye
  length.
