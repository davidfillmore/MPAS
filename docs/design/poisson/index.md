# A finite-volume Poisson solver for atmospheric electrostatics in MPAS-A

**Second-order TRiSK discretization on uniform hexagonal meshes with a
conservative terrain remap**

*David Fillmore*

---

This is the public narrative for the first elliptic Poisson solver
designed for the unstructured Voronoi mesh of the Model for Prediction
Across Scales Atmosphere (MPAS-A). The solver discretizes the Poisson
equation $\nabla \cdot (\varepsilon \nabla \varphi) = -\rho$ using a
TRiSK div-of-grad operator that reuses the existing C-grid mesh metrics,
combined with a centered-difference vertical stencil. After integration
by cell volume the resulting 3D discrete operator is sparse and
symmetric in the standard $\ell^2$ inner product, and its negation is
positive definite (SPD). The linear system is solved by a preconditioned
conjugate-gradient (PCG) iteration with a Jacobi preconditioner,
implemented in-tree and requiring no external linear-algebra library.

Results are reported for method-of-manufactured-solutions (MMS)
verification in planar and spherical geometries, idealized Gaussian
charge distributions, a one-way charge-coupled supercell demonstration,
and terrain-following meshes handled by a conservative altitude-grid
remap. The Cartesian MMS confirms second-order horizontal convergence on
a regular hexagonal mesh. On the sphere the operator is second order
only in the regular-hexagonal bulk: the quasi-uniform icosahedral SCVT
mesh exhibits a localized pentagon-defect limitation (finest-pair $L^2$
slope $1.38$), and a controlled mesh-sensitivity study shows that on the
irregular, general centroidal-Voronoi tessellations that
variable-resolution configurations require, the solved-potential
convergence collapses to a near-stalled rate of $\sim\!0.35$. This is
not caused by degenerate cells — every mesh is well-centered — but
reflects that the operator's second-order behavior depends on a
near-uniform hexagonal regularity that only the special icosahedral
family provides. Recovering it requires a higher-order or finite-element
reconstruction rather than a different finite-volume weight. Idealized
point-charge and tripole diagnostics recover the expected field
topology, while a calibrated supercell stub produces a peak
electric-field magnitude of $4.17 \times 10^4$ V m$^{-1}$ without
two-way feedback. A terrain-following extension solves on a uniform
altitude grid with shaved cut cells and grounded walls, preserving the
SPD operator and, for a gentle hill, the second-order interior rate; the
same supercell repeated over a $1$ km hill gives consistent fields.
Large-scale parallel performance remains a follow-on study. The operator
and solver thus developed provide a foundation for lightning
parameterizations, lightning NO$_x$ emission studies, and a future
global-atmospheric-electric-circuit solver that extends the same
elliptic operator by substituting conductivity $\sigma(\boldsymbol{x})$
for permittivity $\varepsilon$.

:::{admonition} Cite as
:class: note

David Fillmore, *A finite-volume Poisson solver for atmospheric
electrostatics in MPAS-A v8.4: second-order TRiSK discretization on
uniform hexagonal meshes with a conservative terrain remap.*
Preprint, 2026.
[doi:10.5281/zenodo.21173830](https://doi.org/10.5281/zenodo.21173830)
(concept DOI, resolves to the latest version;
[10.5281/zenodo.21173831](https://doi.org/10.5281/zenodo.21173831)
pins version 1). See [](code-and-data.md) for source and archive
details.
:::

```{toctree}
:maxdepth: 2
:caption: Contents

introduction
governing-equations
discretization
solver
verification
idealized-applications
mesh-sensitivity
parallel-performance
end-to-end-supercell
terrain-following
discussion
conclusions
code-and-data
appendices/index
references
```
