# A TRiSK-based Poisson solver for MPAS-A

**Atmospheric electrostatics on the unstructured Voronoi mesh**

*David Fillmore*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19637636.svg)](https://doi.org/10.5281/zenodo.19637636)

---

This is the public design narrative for an elliptic Poisson solver
on the unstructured Voronoi mesh of MPAS-Atmosphere (MPAS-A). The
discrete operator is assembled from the TRiSK div-of-grad
construction
({cite:t}`thuburn2009numerical,ringler2010unified,skamarock2012multiscale`)
in the horizontal and a centered-difference stencil in the
vertical; the resulting sparse symmetric-positive-definite linear
system is solved by a preconditioned conjugate-gradient iteration
with no external linear-algebra dependency. The narrative presents
the mathematical formulation, the algorithm, and the structure of
the five-tier verification programme that will characterize its
accuracy and scaling. Numerical results are reserved for a
follow-on publication.

This narrative is a living document. Major revisions are tagged on
the underlying git repository at
[`github.com/davidfillmore/MPAS`](https://github.com/davidfillmore/MPAS)
and mint a Zenodo DOI via the GitHub↔Zenodo integration. The
formulation-of-record DOI for this first tagged version is
[10.5281/zenodo.19637636](https://doi.org/10.5281/zenodo.19637636)
(`design-poisson-v0.1`); the companion LaTeX manuscript targeting
*Geoscientific Model Development* will cite this DOI.

```{toctree}
:maxdepth: 2
:caption: Contents

introduction
governing-equations
discretization
solver
verification-plan/index
discussion
conclusions
appendices/index
references
```

:::{note}
This site is under active revision. The main-text narrative and
the five-tier verification plan are populated; extended step-by-step
mathematical derivations will be added as further appendix pages
in a subsequent revision.
:::
