Inferring plasma parameters from measurements
---------------------------------------------
The purpose of Langmuir probes is to measure plasma parameters, such as the electron and ion densities, and the electron temperature. The traditional techniques rely on OML theory, which predicts that the characteristic behaves differently in different regions of the probe voltage. For voltages between the floating potential and the background plasma potential (here taken to be zero), for instance, the ion current can be neglected and the OML theory then predicts a slope :math:`\mathrm{d}(\ln I) /\mathrm{d}V` depending only on the electron temperature, as well as some physical constants. Doing a voltage sweep across this region therefore allows the determination of the electron temperature. Further on, in the electron (ion) saturation region, the ions (electrons) can be neglected, and the analytical expressions of the remaining part allows determination of the electron (ion) density, once the electron temperature is known [Marholm2]_, [Bekkeng]_, [MottSmith]_, [Bittencourt]_.

.. image:: OML_regions.png
   :width: 400px

Another technique is that of Jacobsen and Bekkeng for the multineedle Langmuir probe (m-NLP) instrument [Jacobsen]_. The m-NLP instrument consists of at least two (typically four) cylindrical Langmuir probes biased at different fixed positive voltages with respect to a spacecraft. The OML theory predicts that the slope :math:`\mathrm{d}I^2/\mathrm{d}V` depends only on the electron density, except for known physical constants. The m-NLP instrument thus allows inferring the electron density without sweeping the voltage, which gives the m-NLP  instrument faster sampling times and thus higher spatial resolution while spaceborn than swept Langmuir probes. A fast implementation of this density inference technique is readily available in Langmuir. Given an :math:`N\times 4` array ``I``, where each row corresponds to the currents measured at some time instant by probes biased at, say, 2, 3, 4, and 5 volts with respect to some common reference voltage, the densities can be inferred as follows::

    n = jacobsen_density(Cylinder(r=0.255e-3, l=25e-3), [2,3,4,5], I)

The main problem with these approaches is that they rely upon specific analytic expressions for the characteristic, which may not hold for non-ideal cases (finite length, finite radius, collisional or non-Maxwellian plasmas). Another problem is in identifying the different regions. In spaceborn Langmuir probes, for instance, the probe voltage is only known with respect to the spacecraft, and not with respect to the background plasma.

A general formulation
~~~~~~~~~~~~~~~~~~~~~
In relation to the Langmuir software we take a more general point of view [Marholm2]_. Consider that a set of currents :math:`\{\hat I_p\}_{p=1}^N` have been measured in a plasma. The currents may for instance be the currents corresponding to different voltages of a swept Langmuir probe, or it may be the currents collected by different probes, such as in the m-NLP instrument. For the sake of generality, we allow each measurement to obey a different characteristic function, which we denote :math:`I_p(V_p; \mathbf P)`. :math:`V_p` is the probe voltage at which the current was measured, and :math:`\mathbf P` is a vector of other parameters that the characteristics depend upon, such as electron temperature and density. The probe voltage is often not known with respect to the background plasma, but instead with respect to a common reference voltage :math:`V_0`. For spaceborn instruments this is the spacecraft floating potential. With this, we may write :math:`V_p=V_0+V_{0p}`, and the measurements form the following system of equations:

.. math::

    \hat I_p = I_p (V_0 + V_{0p}; \mathbf P),\quad p=1,...,N

This set of equations may be solved for the unknowns :math:`(V_0,\mathbf P)` by any suitable numerical method (insofar as it is well-posed), and there is a wide range of free software available depending on how this system is to be solved. However, it requires programmatic access to the characteristics :math:`I_p(V_p; \mathbf P)`. Computing the currents for given physical parameters may be considered the *forward problem*, and it is a prerequisite for solving the *inverse problem*, namely inferring physical parameters from measured currents. Langmuir focuses on the forward problem. In the following, however, we give a few examples of attacking the inverse problem.

Fitting a power law
~~~~~~~~~~~~~~~~~~~

Fitting the finite-length model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Machine learning approaches
~~~~~~~~~~~~~~~~~~~~~~~~~~~
In progress
