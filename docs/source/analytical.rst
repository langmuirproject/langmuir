Analytical models
=================
- ``OML_current(geometry, species, V=None, eta=None, normalization=None)``
  Current collected by a probe according to the Orbital Motion Limited (OML)
  theory. The model assumes a probe of infinitely small radius compared to
  the Debye length, and for a cylindrical probe, that it is infinitely long.
  Probes with radii up to 0.2 Debye lengths (for spherical probes) or 1.0
  Debye lengths (for cylindrical probes) are very well approximated by this
  theory, although the literature is diverse as to how long cylindrical probes
  must be for this theory to be a good approximation.

- ``thermal_current(geometry, species, normalization=None)``
  Returns the thermal current for the given species and geometry. The
  thermal current is the current the species contributes to a probe at zero
  potential with respect to the background plasma due to random thermal
  movements of particles.

- ``normalization_current(geometry, species)``
  Returns the normalization current for the given species and geometry.
  The normalization current is the current the species would have contributed
  to a probe at zero potential with respect to the background plasma due to
  random thermal movements of particles, if the species had been Maxwellian.

