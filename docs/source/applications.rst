Examples
========

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   beta
   floating_potential
   inverse

Another example is processing 4-Needle Langmuir Probe (4-NLP) measurements using the Jacobsen-Bekkeng method to get density. If the probe biases are 2, 3, 4, and 5V, and the N measurements samples is contained in an Nx4 array ``I``, the densities can be found as follows::

    n = jacobsen_density(Cylinder(r=0.255e-3, l=25e-3), [2,3,4,5], I)
