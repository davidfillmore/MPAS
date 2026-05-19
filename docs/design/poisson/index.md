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
and mint a Zenodo DOI via the GitHub↔Zenodo integration. The current
RTD stable documentation snapshot is release `26.05.1`, archived at
`10.5281/zenodo.20278575` under Zenodo concept DOI
`10.5281/zenodo.19637632`. This snapshot contains the public Sphinx
documentation and committed figure assets; the active solver source
remains on the `feature/poisson` branch pending a later code-and-data
release.

The public site intentionally follows the paper closely for the
mathematical formulation, Phase 1 results, figures, caveats, and
conclusions. It also includes implementation-roadmap material that is
useful for readers tracking active development but too operational for
the journal manuscript: release provenance, Zenodo versioning, planned
Tier C/D work, and Phase 2 extensions.

```{toctree}
:maxdepth: 2
:caption: Contents

introduction
governing-equations
discretization
solver
verification-plan/index
living-document
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
