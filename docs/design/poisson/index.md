# A TRiSK-based Poisson solver for MPAS-A

**Atmospheric electrostatics on the unstructured Voronoi mesh**

*David Fillmore*

---

This is the public design narrative for an elliptic Poisson solver
on the unstructured Voronoi mesh of MPAS-Atmosphere (MPAS-A). The
discrete operator is assembled from the TRiSK div-of-grad
construction
({cite:t}`thuburn2009numerical,ringler2010unified,skamarock2012multiscale`)
in the horizontal and a centered-difference stencil in the
vertical; the resulting sparse symmetric positive definite (SPD)
linear system is solved by a preconditioned conjugate-gradient
(PCG) iteration with no external linear-algebra dependency. The
narrative presents the mathematical formulation, the algorithm,
and Phase 1 results from method of manufactured solutions (MMS)
verification, idealized Gaussian charge diagnostics, and a
one-way charge-coupled supercell demonstration. Mesh-distortion
sensitivity, variable-resolution transition behavior, and
large-scale parallel performance remain planned follow-on studies.

This narrative is a living document. Major revisions are tagged on
the underlying git repository at
[`github.com/davidfillmore/MPAS`](https://github.com/davidfillmore/MPAS)
and mint a Zenodo DOI via the GitHub↔Zenodo integration. The
companion LaTeX manuscript targeting *Geoscientific Model
Development* will cite the first tagged DOI as the citable
formulation-of-record.

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
This site is under active revision. The main-text narrative, the
five-tier verification plan, and the extended mathematical
derivation appendices are populated. Phase 1 results are reported
for Tiers A, B, and E; Tiers C and D remain planned follow-on
studies.
:::
