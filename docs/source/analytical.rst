Analytical theories
-------------------

.. _OML:

The orbital motion-limited (OML) theory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The *orbital motion-limited* (OML) theory is the de-facto standard theory for current collection of a Langmuir probe, and is implemented as the following function::

    OML_current(geometry, species, V=None, eta=None, normalization=None)

Assuming a small probe radius compared to the Debye shielding distance, i.e., the *thick sheath approximation*, the current collected by a probe in a plasma can be calculated based on orbital trajectories within the sheath. The derivations of the implemented analytical expressions assumes:

1. a collisionless plasma
2. a non-drifting Kappa-Cairns-distributed plasma
3. a non-magnetized plasma
4. small probe radius (if spherical or cylindrical) compared to the Debye length
5. large probe length (if cylindrical), or large extent (if planar) compared to the Debye length

It is worth remarking that the Maxwell, Kappa and Cairns distributions are special cases of the Kappa-Cairns distributions and are therefore also covered [Darian]_. Moreover, no assumptions are made on the voltage range.

Radii up to 0.2 Debye lengths (for spherical probes) or 1.0 Debye lengths (for cylindrical probes) typically satisfy the small radius criterion well [Laframboise]_. Cylindrical probes need to be very long, however, to satisfy the long probe criterion [Marholm]_. Other models in Langmuir overcome some of the limitations of the OML theory.

For more information on OML theory, see [MottSmith]_ for the original derivation, or [Darian]_ for the more general derivation with Kappa-Cairns distribution.

.. _thermal-current:

Thermal currents and the sheath-limited (SL) theory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A probe at zero voltage with respect to the background plasma do not accelerate particles, but nevertheless collects a current due to the random thermal motion of particles through the surface of the probe. This is implemented in the following function::

    thermal_current(geometry, species, normalization=None)

This coincides with the OML theory not only for zero voltage, but also for planar probes. This can be understood from the fact that for planar probes, the effective collection area (the sheath boundary) do not increase as the voltage increases and the sheath thickness increases, but remains of the same area.

The *sheath-limited* (SL) theory makes the same assumptions as the OML theory, except that the probe radius is large compared to the Debye length (*thin sheath approximation*). With these approximations, the current collected is :math:`a/r` times the thermal current, where :math:`a` is the distance from the center of the probe to the edge (in some sense) of the sheath surrounding the probe, and :math:`r` is the probe radius. For spheres and cylinders an increase in voltage causes an increase in :math:`a` and hence in collected current. Nevertheless, for very large probes, :math:`a/r` approches unity, and `thermal_current` may be used as an approximation. In this case, the probe can be understod to be locally flat.
