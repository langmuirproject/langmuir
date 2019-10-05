Inverse problems
================
Sometimes the collected current of one or more probes is known and one would like to solve for one or more other parameters. The Langmuir library do not address this analytically in part due to the vast number of such inverse problems, and in part due to some characteristics not being invertible (for instance those who are of tabulated values). However, it is in principle possible to apply numerical methods of root solving, least squares, etc. along with the models in Langmuir.

Consider a cylindrical probe with known dimensions and a positive but unknown voltage collecting a current of -0.4uA in a Maxwellian plasma with known density and temperature. What is the voltage? We shall neglect the current due to ions, and define a residual function. This residual is the difference between the current collected by a probe at a given potential, and the actual collected current, and it is used by a least squares algorithm to compute the voltage::

    >>> from langmuir import *
    >>> from scipy.optimize import leastsq

    >>> sp = Species(n=1e11, T=1000)
    >>> geo = Cylinder(1e-3, 25e-3)
    >>> I = -0.4e-6

    >>> def residual(V):
    >>>     return finite_radius_current(geo, sp, V) - I

    >>> x, c = leastsq(residual, 0)
    >>> print(x[0])
    0.6265540484991013

The reader may verify that this voltage indeed results in the correct current. Notice also that we were in fact able to invert the model ``finite_radius_current``, which consists of tabulated values.

A slightly more interesting inversion problem, is that of determining the ionospheric density from four cylindrical Langmuir probes with known bias voltages with respect to a spacecraft, but an unknown floating potential ``V0`` of the spacecraft with respect to the plasma. We shall assume the bias voltages to be 2.5, 4.0, 5.5 and 7.0 volts. In the below example, we first construct the currents for such probes by assuming a floating potential and a set of plasma parameters, but we do not use this knowledge in the inversion. We do, however, make an initial guess ``x0`` which we believe are somewhat close to the answer::

    >>> from langmuir import *
    >>> from scipy.optimize import leastsq

    >>> geo = Cylinder(1e-3, 25e-3)
    >>> V0 = -0.5
    >>> V = np.array([2.5, 4.0, 5.5, 7.0])
    >>> I = OML_current(geo, Species(n=120e10, T=1000), V+V0)

    >>> def residual(x):
    >>>     n, V0 = x
    >>>     return OML_current(geo, Species(n=n, T=1500), V+V0) - I

    >>> x0 = [10e10, -0.3]
    >>> x, c = leastsq(residual, x0)
    >>> n, V0 = x

    >>> print(n)
    1199899818493.931

    >>> print(V0)
    -0.5417515655165968

The method correctly determined the density to be 120e10. However, the floating potential ``V0`` is off by almost ten percent. The reason is that the temperature is considered unknown, and assumed to be 1500K when solving the problem, while it is actually 1000K. Since we have four measurements (four equations) and only two unknowns, it is tempting to also include the temperature as an unknown parameter and try to solve for it. However, if this is done the least squares algorithm will fail miserably. The reason is that the set of equations arising for the attracted-species current of cylindrical probes are singular and cannot be solved for even analytically. Fortunately, both the temperature and floating potential can be eliminated from the equation when analytically solving for the density, and similarly it also works to obtain the density from the least squares algorithm. Since the floating potential and temperature represent a coupled unknown which cannot be solved for, an error in assuming one is reflected as an error in the other.

This demonstrates the usefulness as well as challenges and subtleties of solving inverse Langmuir problems.

