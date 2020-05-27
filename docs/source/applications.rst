Examples
========

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   normalization
   beta
   floating_potential
   inverse


As an example, the following snippet computes the normalized electron current of a probe of 3 Debye lengths radius and normalized voltage of -10::

    >>> sp  = Species(n=1e11, T=1000)
    >>> geo = Cylinder(r=3*sp.debye, l=1)
    >>> I = finite_radius_current(geo, sp, eta=-10, normalization='th')

Notice that setting ``l==1`` means you get the current per unit length.
