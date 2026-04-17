# Embedding SHDOM in MPAS-A

**3D radiative transfer via conservative coupling of a Cartesian sub-domain to an unstructured Voronoi mesh**

*David Fillmore*

---

This is the public design narrative for the SHDOM–MPAS-A coupling
framework: a hybrid 1D/3D radiative-transfer architecture that
embeds the Spherical Harmonic Discrete Ordinate Method (SHDOM,
{cite:t}`evans1998shdom`) on a Cartesian sub-domain within the
unstructured Voronoi host mesh of MPAS-Atmosphere (MPAS-A). The
document describes the mathematical formulation and the algorithmic
architecture of the framework; the implementation plan, validation
plan, and week-to-week working notes are kept in a separate
planning repository and are not part of this public narrative.

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
methods/index
discussion
conclusions
appendices/index
references
```

:::{note}
This site is under active revision. The main-text narrative is
populated; the extended Evans-style SHDOM derivations will be added
as further appendix pages in a subsequent revision.
:::
