Characteristic models
=====================

Langmuir comes with several models for the characteristic :math:`I(V)` of Langmuir probes. Each model is represented by a function which, in addition to the voltage ``V``, take a ``geometry`` and a ``species`` argument. The models have similar signatures and can often be interchanged. ``geometry`` is one of the previously mentioned probe geometries (``Plane``, ``Cylinder`` or ``Sphere``), and the ``species`` argument is either a single ``Species`` object, or a list of ``Species``. If it is a single species, the return value is the current collected from that species alone. This may be sufficient in regimes where the current is dominated by the collection of one species. If the ``species`` argument it is a list of species, it is the sum of the currents from all species (there is currently no model which accounts for non-linear coupling between species). If both ions and electrons are included, the characteristics captures both the electron saturation, electron retardation and ion saturation regimes. The model functions also have arguments ``eta`` and ``normalization`` which will be covered in :doc:`normalization`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   analytical
   normalization
   finite_length
   finite_radius
