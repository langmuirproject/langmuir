Normalization
-------------

The characteristics of Langmuir probes do not depend on every conceivable parameter such as density, temperature, probe length, etc. independently, but instead upon a smaller set of non-dimensional groups of such parameters [Laframboise]_, [Marholm]_. One such group is the normalized voltage:

.. math::

    \eta = -\frac{qV}{kT}
    
where :math:`q` and :math:`T` are the charge and temperature of the collected species, :math:`V` is the voltage, and :math:`k` is Boltzmann's constant (for :math:`\eta>0` the species is attracted, and for :math:`\eta<0` it is repelled). It may therefore be of interest to study an entire class of problems with the same normalized parameters instead of a specific case. In Langmuir one may specify normalized voltages by using the argument ``eta`` instead of ``V`` in models such as ``OML_current``.

Similarly, a normalized current :math:`I/I_0` may be defined, where :math:`I` is the collected current and :math:`I_0` is a characteristic current. The models in Langmuir may return currents normalized by one of several possible characteristic currents depending on the value of the ``normalization`` argument:

- ``'th'``: Normalized by the :ref:`thermal current <thermal-current>`. Often the most natural choice.
- ``'oml'``: Normalized by the current according to :ref:`the OML theory <OML>`. This is useful for comparing other models with OML theory.
- ``'thmax'``:  Normalized by the thermal current of a Maxwellian plasma regardless of what the distribution actually is. Probably only useful for comparison with [Darian]_.

Finally, all lengths are normalized by the Debye length :math:`\lambda_D`. Below is a complete example of obtaining the normalized current for a cylindrical probe of radius :math:`0.2\lambda_D` and length :math:`10\lambda_D` with :math:`\eta=10`::

    >>> sp = Species()
    >>> geometry = Cylinder(r=0.2*sp.debye, l=10*sp.debye)
    >>> OML_current(geometry, sp, eta=10, normalization='th')
    3.7388259506315147
  
Since only the non-dimensional groups determine the normalized collected current, we do not care about the exact parameters of the species, but leave them at the default. Note that in this case it actually does not matter what the probe size is, because the thermal current depends on the probe size in the same way as the current predicted by OML theory. This may differ for other models, however. See also :ref:`example_normalization`.

Note that the non-dimensional groups is specific to each species, e.g., the voltage normalized with respect to electrons is not the same as with respect to ions. If a multi-species plasma is specified, the normalization used will be with respect to the first species in the list. E.g., if electrons are the first element in the list, :math:`\eta=eV/kT_e` where :math:`e` is the elementary charge and :math:`T_e` is the electron temperature.
