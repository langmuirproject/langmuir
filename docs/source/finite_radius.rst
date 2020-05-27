Finite radius models
====================
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


