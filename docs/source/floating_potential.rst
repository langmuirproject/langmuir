Numerically solving for the floating potential
----------------------------------------------
In this example we consider how to numerically obtain the floating potential by means of a root solver. We shall consider a sphere situated in a hydrogen plasma, which is known to have a floating potential of :math:`-2.5kT/e`, where :math:`k` is Boltzmann's constant, :math:`T` is the temperature, and :math:`e` is the elementary charge [Whipple]_. The code for the example is given below:

.. literalinclude:: ../../demo/floating_potential.py

The requirement for a steady floating potential is that the net current into the object, by all species, is zero. Otherwise the potential would increase or decrease. Thus, in our case we seek to find the root of the characteristic, which we take as ``OML_current`` in the example above. The function ``root_scalar`` from SciPy iterates on the value ``V`` of ``res`` within the range :math:`(-3,0)`, until it returns a value that is sufficiently close to zero, i.e., the root of ``OML_current``. The script returns::

    -0.2157892685234552
    -0.21543333155362945

and hence is in excellent agreement with previous results.
