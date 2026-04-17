(sec-conclusions)=
# Conclusions

This narrative has described a hybrid 1D/3D radiative-transfer
framework for the MPAS-A unstructured-mesh atmospheric model, in
which a 3D solver (SHDOM) operates on a Cartesian sub-domain
embedded within the Voronoi host mesh while a 1D solver (RRTMG)
continues to service the remainder of the domain. The framework
is organized around four architectural invariants —
per-(cell, band) authoritative solver, domain-boundary 1D
treatment, persistent diagnostic output, and clean module
interfaces — that were chosen so that a scoped five-phase roadmap,
from the Phase 1 single-band full-domain diagnostic through to
Phase 5 adaptive cloud-aware 3D patches, proceeds by extension
rather than by refactoring.

The numerical content of the framework is concentrated in the
Voronoi-to-Cartesian mesh mapping of
[](methods/mesh-mapping.md), whose area-weighted conservative
remapping formula preserves the physical bounds on
single-scattering albedo and asymmetry parameter by construction
rather than by post hoc clipping. The RRTMG patching of
[](methods/optical-properties.md) exposes optical properties and
per-band fluxes through add-only entry points that leave every
existing upstream routine byte-identical, and the MPI
redistribution of [](methods/mpi.md) resolves SHDOM's Cartesian
pencil decomposition against MPAS-A's METIS graph partition
through a precomputed `MPI_Alltoallv` plan.

The broader significance of the framework lies beyond its
first-phase result itself. The architectural separation between
unstructured host and structured 3D solver developed here is, to
the author's knowledge, the first such coupling to be formulated
for a production atmospheric model with unstructured horizontal
discretization; the components are interchangeable enough to
accommodate alternative 3D solvers (SPARTACUS, TenStream, MYSTIC)
with modest specialization. The roadmap through Phase 5 terminates
at adaptive cloud-aware 3D RT on variable-resolution global meshes
— a regime in which the inclusion of 3D radiative effects is
beginning to have measurable predictive consequence, and in which
the engineering cost of 3D RT must be paid only where clouds
demand it.
