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
rarely result in analytical expressions $I(V)$, and this makes them
inconvenient to apply for other researchers.

Langmuir is a library of Python functions acting as models of the
characteristics $I(V)$, some of which are merely analytical expressions for
idealized cases, while others are based on high-fidelity simulation results,
and cover non-ideal cases such as probes of finite length [@marholm] or
non-Maxwellian plasmas [@darian]. The simulation-based models incorporate
appropriate interpolation schemes between data points which make them
continuous within their domain. They are not limited to physical parameters of
a particular range, but are normalized internally, and scaled according to the
$\pi$-theorem of dimensional analysis [@buckingham]. Evaluating a
characteristic for a set of parameters is near-instantaneous compared to
running a high-fidelity simulation, which can take hours for each input. As
such, the simulation-based models may be seen as a simple kind of reduced-order
models.

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

Instead, we take a more general point of view, as suggested by e.g. @marholm. A
set of measurements $\{\hat I_p\}_{p=1}^N$ taken at voltages $\{V_p\}_{p=1}^N$
by probe(s) with a characteristic $I(V; \mathbf P)$, should follow this system
of equations:
$$
    \hat I_p = I(V_p; \mathbf P),\quad p=1,...,N
$$
Finding the unknowns (e.g., $\mathbf P$) such that this system is satisfied
is referred to as the *inverse problem*. The inverse problem can in principle
be solved by any suitable computational method, which is already available in
other software packages, but programmatic access to the characteristic is a
prerequisite. Langmuir enables such computational approaches by providing the
characteristics. In addition, the tutorial-like documentation of Langmuir is an
important part of the contribution, as it not only documents the features of
Langmuir, but also gives examples of how to solve the inverse problem.

Langmuir has already been used to solve real-world problems [@marholm;
@guthrie], and continue to be used in ongoing research.

# Acknowledgements

This work is partly funded by ..., and is partly not funded.

# References
