Specifying the geometry
=======================
Langmuir supports three probe geometries, with self-descriptive names:

- ``Plane(A)`` represents a planar probe with surface area ``A``.
- ``Cylinder(r, l, lguard=0, rguard=0)`` represents a cylindrical probe of radius ``r`` and length ``l``. The optional arguments ``lguard`` and ``rguard`` may be used to specify the length of the left and right guards, respectively. Setting them to ``float('inf')`` or ``True`` means that there is an ideal guard.
- ``Sphere(r)`` represents a spherical probe of radius ``r``.

All dimensions are in SI units. Not all models may be able to support all geometries, and only the finite-length model makes use of the guard feature.
