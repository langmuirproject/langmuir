---
title: 'Langmuir: Programmatically accessible current--voltage characteristics for ideal and non-ideal Langmuir probes'
tags:
  - Python
  - Langmuir probes
  - Plasma diagnostics
  - Computational physics
  - Model order reduction
authors:
  - name: Sigvald Marholm
    orcid: 
    affiliation: 
  - name: Diako Darian
    orcid:
    affiliation: 
affiliations:
 - name: 
   index: 1
 - name: 
   index: 2
date: 
bibliography: paper.bib
---

# Summary

Langmuir probes are among the most common instruments for measuring plasma
parameters such as the electron density and the electron temperature in space
and laboratory plasma devices [@tbd:3000]. The instrument works by immersing a
conductor of a certain voltage in a plasma (e.g., the ionosphere), and
measuring the current (charged particles) collected from the plasma. Since the
current depend on the plasma parameters, it is in principle possible to infer
plasma parameters from a set of such measurements at different voltages.

Traditionally, the current--voltage characteristic $I(V)$ (current as function
of voltage) for such instruments is derived by means of analytical theories
such as the *orbital motion-limited* (OML) theory [@tbd:3000]. Such theories
make several simplifying assumptions, for instance that cylindrical probes be
infinitely long. With computer simulations it is possible to simulate
characteristics with fewer restrictions, but simulation results are rarely in
the form of an analytical expression $I(V)$ that can easily be implemented
within a computer program, and made use of by others.

Langmuir is a library of Python functions that act as characteristics $I(V)$,
some of which are merely idealized, analytical expressions, while others cover
non-ideal cases based on high-fidelity computer simulation results. The latter
class incorporate appropriate interpolation schemes which make them continuous
within their domain, as well as scaling based on dimensional analyses
(Buckingham's pi theorem [@tbd:3000]) which make the results scale to a wide
range of physical scenarios. Evaluating a characteristic for a set of
parameters is near-instantaneous compared to running a high-fidelity
simulation, which can take hours for each data point. As such these functions
may be seen as simple reduced order models [@tbd:3000].

# Statement of need

When parameters are inferred using equations from these theories, 

# Acknowledgements

This work is partly funded by ..., and is partly not funded.

# References
