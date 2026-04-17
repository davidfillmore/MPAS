(sec-discussion)=
# Discussion

## Framework extensibility across phases

The Phase 1 configuration described above is the smallest complete
exercise of the full coupling pipeline. The architectural
invariants of [](methods/framework.md) were chosen so that four
additional scoped phases extend the framework incrementally
without refactoring.

**Phase 2 — multi-band diagnostic.** The pipeline generalizes to
$N_b$ bands by replacing the single-band Phase 1 diagnostic fields
with a band-indexed array of the same diagnostic quantities, by
extending the RRTMG optical-property entry point to return all
active bands rather than the single visible band, and by extending
the SHDOM interface to loop over bands per dispatch. The band loop
can be serial or, if profiling motivates it, parallelized. The
mesh mapper, the heterogeneous MPI redistribution, and the
dispatcher structure are unchanged.

**Phase 3 — 3D/1D interior interface.** The Phase 1 cell mask
takes the trivial value "all cells". Phase 3 flips an arbitrary
sub-domain of the MPAS-A mesh to SHDOM for selected bands and
runs 1D RT on the remainder, introducing the framework's first
interior 3D/1D interface. The outstanding design question is the
prescription of the SHDOM edge-radiance boundary condition from
the adjacent 1D columns: in the two-stream RRTMG view the lateral
radiance field is formally absent, and a closure is required to
reconstruct the directional radiance SHDOM needs from the scalar
flux the 1D solver provides. A conservative default is to impose
the plane-parallel radiance implied by the adjacent 1D column's
fluxes under the Henyey–Greenstein phase function assumption; more
sophisticated closures using spatial information from neighbouring
1D columns are possible and will be evaluated against Phase 3
validation runs.

**Phase 4 — variable-resolution MPAS-A mesh.** On a
variable-resolution MPAS-A mesh, a refined sub-domain runs 3D RT
while the coarser surrounding mesh continues with 1D. The mesh
mapper must now handle non-axis-aligned Voronoi cells of
non-uniform size against a uniform Cartesian SHDOM grid, and the
vertical discretization becomes terrain-following rather than
flat. Both extensions are strict additions to the Phase 1
algorithm: the area-weighted formula
{eq}`eq-v2c-c2v` carries over unchanged, with a general
polygon-intersection routine replacing the hex-rectangle
specialization and a vertical remapping stage composed after the
horizontal stage.

**Phase 5 — adaptive cell-set selection.** The Phase 3 cell mask
is static; Phase 5 makes it dynamic, reselecting SHDOM-active
cells at each radiation step based on the instantaneous cloud
field (for instance, cells with any cloud condensate above a
threshold within their column). This requires only reconstruction
of the mapper weights when the active cell set changes and is
transparent to the remaining modules.

**Later — feedback, full-band coverage, longwave.** Beyond the
scoped Phase 1–5 roadmap, three extensions await physical
motivation and scoping. Feedback of SHDOM-derived heating rates
into the dynamical core — inverting the Phase 1 diagnostic
invariant for selected (cell, band) pairs — is unblocked once
Phase 2 and Phase 3 have established confidence in the numerics.
Extending SHDOM's band coverage to the full RRTMG band set retires
the 1D approximation entirely within the 3D region and brings the
framework to parity with column-wise operational RT. A
corresponding longwave path — a second add-only RRTMG entry
point, a longwave dispatcher, a longwave diagnostic schema —
parallels the shortwave path through the same pipeline with no
new research questions beyond spectral coverage.

## Limitations

Four limitations of the Phase 1 configuration are scoped
intentionally and deserve explicit statement. First, the nominal
1 km cell spacing of the MPAS-A supercell test case is coarse
relative to the sub-cloud-scale structure on which 3D radiative
effects develop most strongly; SHDOM's own guidance
{cite:p}`evans1998shdom` favours finer grids on cloud-field
inputs. The Phase 1 outputs therefore characterize the 3D effect
as it manifests *at the resolution the host model provides*, and
should not be interpreted as an upper bound on the full 3D signal
that a higher-resolution host would produce. Second, the visible
band is the only spectrum exercised; the longwave path, the full
shortwave band set, and the near-infrared cloud regime are all
scoped to later phases. Third, SHDOM is diagnostic only: the
feedback to MPAS-A dynamics that would make SHDOM the
authoritative solver within the 3D region is scoped out of Phase 1
and requires the per-band confidence that Phase 2 and the
sub-domain confidence that Phase 3 will establish. Fourth,
bypassing MCICA for the SHDOM optical-property pathway trades some
sub-column overlap information against the determinism required
for like-with-like SHDOM-vs-RRTMG comparison; reinstating a
stochastic sub-column treatment for SHDOM, if needed at higher
resolution, is a tractable Phase 2+ option.

## Relationship to other 3D-RT-in-ESM efforts

The recent literature on 3D radiative transfer in atmospheric
models partitions along two axes: the method family the 3D solver
belongs to, and the host model's horizontal discretization. Along
the first axis, the present framework uses the spherical-harmonic
discrete-ordinate formulation of {cite:t}`evans1998shdom`, which
returns the full directional radiance field; the SPARTACUS
approach of {cite:t}`hogan2016spartacus` retains a two-stream
angular discretization but encodes lateral transport through a
matrix correction, and is operational in the ECMWF IFS; the
TenStream approach of
{cite:t}`jakub2015tenstream,jakub2016tenstream_les` discretizes
ten streams and has been coupled to UCLA-LES and ICON LES
configurations; MYSTIC {cite:p}`mayer2009mystic` uses Monte Carlo
integration and is the reference against which approximate methods
are most often calibrated. The four method families differ in
computational cost (Monte Carlo most expensive;
two-stream-plus-correction least) and in what they return
(radiance fields for SHDOM and MYSTIC; fluxes and heating rates
for SPARTACUS; ten-stream intensities for TenStream).

Along the second axis — the host model's horizontal mesh — every
coupling of a 3D solver to an atmospheric host model known to the
author operates on a structured or logically-rectangular grid:
SPARTACUS in IFS on a reduced Gaussian grid, TenStream in UCLA-LES
and ICON on Cartesian LES domains and on the structured ICON
native grid. The unstructured-host case that MPAS-A presents has
no published precedent. This narrative's contribution is not a
new 3D solver but rather the mesh-coupling infrastructure that
connects an existing 3D solver to an unstructured host: the
conservative bound-preserving Voronoi-to-Cartesian optical-
property map of [](methods/mesh-mapping.md), the add-only RRTMG
patching strategy of [](methods/optical-properties.md), and the
heterogeneous MPI redistribution of [](methods/mpi.md). These
components are, with minor specializations, reusable for any of
the alternative 3D solver choices; the 3D engine is
interchangeable inside the framework as described.
