Computing fundamental plasma parameters
=======================================
The Langmuir library is also convenient to use for quick computations of fundamental plasma parameters, since ``Species`` computes these upon initialization. For instance, to get the electron Debye length of a plasma with a certain density and temperature::

    >>> Species(n=1e11, eV=0.1).debye
    0.007433942027347403

The following member variables and methods are accessible in ``Species``:

+-----------------+---------------------------------+
| Member          | Description                     |
+=================+=================================+
| ``debye``       | The Debye length                |
+-----------------+---------------------------------+
| ``omega_p``     | The angular plasma frequency    |
+-----------------+---------------------------------+
| ``freq_p``      | The linear plasma frequency     |
+-----------------+---------------------------------+
| ``period_p``    | The plasma period               |
+-----------------+---------------------------------+
| ``omega_c(B)``  | The angular cyclotron frequency |
+-----------------+---------------------------------+
| ``freq_c(B)``   | The linear cyclotron frequency  |
+-----------------+---------------------------------+
| ``period_c(B)`` | The cyclotron period            |
+-----------------+---------------------------------+
| ``larmor(B)``   | The larmor radius               |
+-----------------+---------------------------------+

The latter four members are methods which take the magnitude of the magnetic flux density as an argument. In addition, every valid keyword argument is also a valid member variable::

    >>> Species(n=1e11, T=1000).eV
    0.08617330337217212

Finally, the total Debye length of a plasma consisting of multiple species can be obtained using the ``debye()`` function::

    >>> plasma = []
    >>> plasma.append(Species('electron' , n=1e11, T=1000))
    >>> plasma.append(Species(amu=16, Z=1, n=1e11, T=1000))
    >>> debye(plasma)
    0.004879671013271479
