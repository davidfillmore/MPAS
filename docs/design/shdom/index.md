# Embedding SHDOM in MPAS-A

**3D radiative transfer via conservative coupling of a Cartesian sub-domain to an unstructured Voronoi mesh**

*David Fillmore*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19637633.svg)](https://doi.org/10.5281/zenodo.19637633)

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
formulation-of-record DOI for this first tagged version is
[10.5281/zenodo.19637633](https://doi.org/10.5281/zenodo.19637633)
(`design-shdom-v0.1`); the companion LaTeX manuscript targeting
*Geoscientific Model Development* will cite this DOI.

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
