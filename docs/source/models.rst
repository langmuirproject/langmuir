Models for collected current
============================
Langmuir comes with several models for the collected current. Each model is represented by a function which takes a ``geometry`` and a ``species`` argument. The ``geometry`` is one of the above probe geometries, and the ``species`` parameters is either a single ``Species`` object or a list of such if it is desirable to take into account the effect of all species in a plasma. Most models also depend on the potential of the probe with respect to the background plasma. The potential can either be specified in volts using th ``V`` argument, or in terms of normalized voltage e*V/(k*T) using the ``eta`` argument. If these are Numpy arrays, the output will be a Numpy array of collected currents. The ``normalization`` argument can be set to ``'th'``, ``'thmax'``, ``'oml'`` to return currents normalized by the thermal current, Maxwellian thermal current, or OML current, respectively, rather than Ampére. Below is a description of all models:

- ``OML_current(geometry, species, V=None, eta=None, normalization=None)``
  Current collected by a probe according to the Orbital Motion Limited (OML)
  theory. The model assumes a probe of infinitely small radius compared to
  the Debye length, and for a cylindrical probe, that it is infinitely long.
  Probes with radii up to 0.2 Debye lengths (for spherical probes) or 1.0
  Debye lengths (for cylindrical probes) are very well approximated by this
  theory, although the literature is diverse as to how long cylindrical probes
  must be for this theory to be a good approximation.

- ``finite_radius_current(geometry, species, V=None, eta=None, table='laframboise-darian-marholm', normalization=None)``
  A current model taking into account the effects of finite radius by
  interpolating between tabulated normalized currents. The model only
  accounts for the attracted-species currents (for which eta<0). It does
  not extrapolate, but returns ``nan`` when the input parameters are outside
  the domain of the model. This happens when the normalized potential for any
  given species is less than -25, when kappa is less than 4, when alpha is
  more than 0.2 or when the radius is more than 10 or sometimes all the way
  up towards 100 (as the distribution approaches Maxwellian). Normally finite
  radius effects are negligible for radii less than 0.2 Debye lengths (spheres)
  or 1.0 Debye lengths (cylinders).

- ``finite_length_current(geometry, species, V=None, eta=None, normalization=None)``
  The Marholm-Marchand model for finite-length probes. Works for normalized
  voltages up to 100.

- ``finite_length_current_density(geometry, species, V=None, eta=None, z=None, zeta=None, normalization=None)``
  The current per unit length according to the Marholm-Marchand model for
  finite-length probes. Works for normalized voltages up to 100. ``z`` is
  position on the probe, and ``zeta`` is position normalized by the Debye
  length.

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

As an example, the following snippet computes the normalized electron current of a probe of 3 Debye lengths radius and normalized voltage of -10::

    >>> sp  = Species(n=1e11, T=1000)
    >>> geo = Cylinder(r=3*sp.debye, l=1)
    >>> I = finite_radius_current(geo, sp, eta=-10, normalization='th')

Notice that setting ``l==1`` means you get the current per unit length.
