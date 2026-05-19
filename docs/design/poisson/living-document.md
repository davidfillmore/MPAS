(sec-living-document)=
# Living document and release tracking

This Read the Docs site is the public, living companion to the
Poisson-solver manuscript. The manuscript is the compact journal
record: it states the formulation, reports the bounded Phase 1 results,
and identifies the limitations. This site mirrors that scientific
content, including the Phase 1 figures and caveats, but also keeps the
development record that would be too operational for the paper.

## Relationship to the paper

The paper and this site use the same core scope:

- a TRiSK horizontal div-of-grad operator and centered vertical stencil
  for $\nabla \cdot (\varepsilon \nabla \phi) = -\rho$;
- an in-tree symmetric positive definite (SPD) Poisson operator solved
  by preconditioned conjugate gradients (PCG);
- completed Phase 1 results for Tier A method-of-manufactured-solutions
  (MMS) verification, Tier B idealized charge diagnostics, and Tier E
  one-way charge-coupled supercell output;
- planned follow-on Tier C mesh-sensitivity and Tier D scaling studies;
- Phase 2 deferral of terrain-following coordinates, defect-aware
  spherical correction, physical electrification, lightning discharge,
  and two-way feedback.

The site may be updated more frequently than the manuscript while the
solver matures. When that happens, the versioned Zenodo record is the
stable reference point for a particular public state of the method.

## Zenodo tracking

Tagged releases of the public MPAS fork are archived through the
GitHub-Zenodo integration. Each tagged release receives a citable DOI.
Zenodo creates one concept record per GitHub repository, so the SHDOM
and Poisson design releases are versioned snapshots under the shared
MPAS-fork Zenodo concept. The initial Poisson design release is:

- v0.1 tag: `design-poisson-v0.1`;
- v0.1 DOI: `10.5281/zenodo.19637636`;
- RTD project: `mpas-poisson`;
- source repository: [`davidfillmore/MPAS`](https://github.com/davidfillmore/MPAS).

The manuscript DOI, when minted, should be cited for the peer-reviewed
paper. The version-specific Zenodo DOI should be cited when a reader
needs the exact public source/documentation snapshot associated with a
solver state, benchmark figure, or RTD release.

## Planning information retained here

The journal paper deliberately separates completed Phase 1 results from
future work. This site keeps the same distinction but gives additional
planning detail:

- **Tier C** will quantify sensitivity to perturbed quasi-uniform
  meshes, variable-resolution transition bands, and the documented
  spherical pentagon neighborhoods.
- **Tier D** will measure strong and weak scaling, PCG iteration counts,
  communication cost, and preconditioner behavior.
- **Phase 2 numerical methods** will revisit pentagon defects with a
  correction that preserves conservation, symmetry, and solver
  compatibility.
- **Phase 2 physics and geometry** will add terrain-following metric
  terms, physical electrification, lightning discharge, charge leakage,
  and two-way coupling only after the Phase 1 operator and diagnostic
  path are stable.

Those planning notes are part of the public development record, not
claims of completed verification.
