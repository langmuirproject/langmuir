---
title: 'Langmuir: Programmatically accessible current--voltage characteristics for ideal and non-ideal Langmuir probes'
tags:
  - Python
  - Langmuir probes
  - Plasma diagnostics
  - Computational physics
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
and laboratory plasma devices. The instrument works by immersing a
conductor of a certain voltage in a plasma (e.g., the ionosphere), and
measuring the current (charged particles) collected from the plasma. Since the
current depend on the plasma parameters, it is in principle possible to infer
plasma parameters from a set of such measurements at different voltages.

Traditionally, the current--voltage characteristic $I(V)$ (current as function
of voltage) for such instruments is derived by means of analytical theories
such as the *orbital motion-limited* (OML) theory [@mott-smith]. Such theories
make several simplifying assumptions, for instance that cylindrical probes be
infinitely long and very thin. With computer simulations it is possible to
simulate current collection with fewer restrictions. However, simulations
rarely result in analytical expressions $I(V)$, and this make them inconvenient
to apply for other researchers.

Langmuir is a library of Python functions that act as characteristics $I(V)$,
some of which are merely analytical expressions for idealized cases, while
others cover non-ideal cases (finite-length probes, non-Maxwellian
distributions, etc.) based on high-fidelity simulation results. The latter
class incorporate appropriate interpolation schemes which make them continuous
within their domain, as well as scaling based on dimensional analyses
(Buckingham's pi theorem [@buckingham]) which make the results scale to a
different orders of magnitudes of physical parameters than the originating
simulations. Evaluating a characteristic for a set of parameters is
near-instantaneous compared to running a high-fidelity simulation, which can
take hours for each data point. As such, these functions may be seen as a
simple kind of reduced-order models.

# Statement of need

Obtaining the current $I(V; \mathbf P)$ for a probe voltage $V$ and a given set
of plasma parameters contained in the vector $\mathbf P$ may be seen as the
*forward problem*, and is interesting in its own right. It can for instance be
used to study how probes should be designed in order to act in a certain way.
More importantly, however, Langmuir enables the use of computational methods to
infer plasma parameters.

Conventionally, plasma parameters are inferred using expressions derived from
analytical theories. The OML theory for cylindrical probes, for example,
predict that for probes with large enough voltage, the slope
$\mathrm{d}I^2/\mathrm{d}V=Cn_e$ where $C$ is a known constant and $n_e$ is the
electron density. Two or more probes can thus be used to calculate the slope
and in turn the density [@jacobsen]. The accuracy of such methods is of course
limited by the assumptions made in the analytical derivations.

Instead, we take a more general point of view, as proposed by @marholm. A set
of measurements $\{\hat I_p\}_{p=1}^N$ taken at voltages $\{V_p\}_{p=1}^N$ by
probe(s) with a characteristic $I(V; \mathbf P)$, should follow this system of
equations:
$$
    \hat I_p = I(V_p; \mathbf P),\quad p=1,...,N
$$
Finding the unknowns (typically $\mathbf P$) such that this system is satisfied
is referred to as the *inverse problem*. The inverse problem can in principle
be solved by any suitable computational method, which is already available in
other software packages, but programmatic access to the characteristic is a
prerequisite. This is where Langmuir comes in, as an enabling technology. In
addition, the tutorial-like documentation is an important part of the
contribution, as it serves to demonstrate how Langmuir may be used together
with third party software to solve the inverse problem.

Langmuir has already been used to solve real-world problems [@marholm,
@guthrie], and continue to be used in ongoing research.

# Acknowledgements

This work is partly funded by ..., and is partly not funded.

# References
