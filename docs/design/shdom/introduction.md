(sec-introduction)=
# Introduction

## Motivation

Plane-parallel (1D) radiative transfer underpins essentially every
operational weather and climate model. The computational economics
are inescapable: a correlated-$k$ solver such as RRTMG
{cite:p}`mlawer1997rrtm,iacono2008rrtmg`, combined with a Monte Carlo
independent-column approximation (MCICA) for sub-grid cloud overlap
{cite:p}`pincus2003mcica,barker2003monte`, resolves each model column
independently along the local vertical and thereby reduces a
fundamentally 3D transport problem to $N_{\mathrm{cell}}$ decoupled
1D solves. This independence across columns is what makes radiation
affordable to call at every physics timestep; a full 3D solver on
the same domain is several orders of magnitude more expensive per
call.

The physical price of the 1D approximation is by now well
characterized. Horizontal photon transport — side-illumination of
convective towers, cloud-edge brightening, and inter-cloud
scattering — is suppressed by construction in the plane-parallel
framework and therefore omitted entirely from column-by-column
solvers
{cite:p}`marshak2005,varnai1999effects,barker1996gamma,petty2006`.
The resulting errors are small when individual clouds are sub-grid,
which is the regime for which RRTMG and its relatives were designed;
they grow rapidly once horizontal cloud structure is resolved,
because in that regime the neglected 3D transport operates precisely
on the scales the dynamical core is advecting.

Model resolutions are now passing through the cloud-resolving
threshold globally. Kilometre-scale atmospheric simulations are
feasible within a single integration, and convection-permitting
configurations of Earth system models are entering production use.
The Deep Convective Clouds and Chemistry (DC3) campaign
{cite:p}`barth2015dc3` both motivated and made available
observational constraints on the structure of exactly the scenes in
which 1D error is no longer negligible. Full 3D radiative transfer
at every cell, at every timestep, remains prohibitive for production
simulations; what is feasible is a *hybrid 1D/3D* arrangement in
which a 3D solver operates only where resolved cloud structure
demands it, and a 1D solver continues to service the remainder of
the domain. The central unresolved problem for such an arrangement
on a modern unstructured-mesh atmospheric model — and the problem
this narrative addresses — is the coupling of a 3D solver's natural
structured-grid representation to the host's unstructured Voronoi
cells.

## Prior work

### MPAS-A as an unstructured-mesh host

The Model for Prediction Across Scales atmospheric core
{cite:p}`skamarock2012mpas` discretizes the governing equations on a
centroidal Voronoi tessellation of the sphere with C-grid
staggering. The representation of geostrophic modes
{cite:p}`thuburn2009cgrid` and the unified treatment of energy and
potential-vorticity conservation {cite:p}`ringler2010unified` derive
from the TRiSK family of operators, which extend classical
finite-difference C-grid schemes to arbitrary polygonal meshes. The
multiresolution capability {cite:p}`ringler2011multiresolution` is
the direct antecedent of MPAS-A's variable-resolution atmospheric
configurations in current production, and will become the principal
setting in which a 3D radiation patch embedded on a refined
sub-domain is physically motivated.

### SHDOM and contemporaneous 3D solvers

The Spherical Harmonic Discrete Ordinate Method of
{cite:t}`evans1998shdom` resolves the full radiance field in three
dimensions by expanding the source function in spherical harmonics
while sweeping radiances along discrete ordinates, using each
representation where it is cheapest; an adaptive Cartesian cell
structure concentrates degrees of freedom on regions of strong
radiative gradient. Contemporaneous and more recent 3D solvers in
atmospheric science span a range of method families: Monte Carlo
integration of the radiative transfer equation (the MYSTIC family,
{cite:alp}`mayer2009mystic`), two-stream schemes augmented with a
tractable matrix representation of lateral transport (SPARTACUS,
{cite:alp}`hogan2016spartacus`), and discrete-ordinate solvers
aimed at cloud-resolving LES (TenStream,
{cite:alp}`jakub2015tenstream,jakub2016tenstream_les`). Of these,
SHDOM provides the full directional radiance field at moderate cost
and with a well-documented open implementation, and is the solver
embedded in the present framework; we return to its algorithmic
structure in [](methods/shdom-algorithm.md).

### 3D radiative transfer in atmospheric host models

Every published coupling of a 3D solver to an atmospheric host model
known to the author operates on a structured or
logically-rectangular horizontal grid. TenStream has been integrated
with UCLA-LES {cite:p}`jakub2016tenstream_les`; SPARTACUS is used
operationally in the ECMWF IFS {cite:p}`hogan2016spartacus`; MYSTIC
is coupled to libRadtran and to various LES hosts
{cite:p}`mayer2009mystic`. None of these integrations face the
structured/unstructured boundary that a Voronoi host model presents.
The unstructured-host case is the open problem addressed here.

### Conservative mesh remapping

The numerical machinery for mapping fields between dissimilar meshes
was established for climate modelling by {cite:t}`jones1999scrip` in
the form of first- and second-order conservative remapping on
spherical grids, and has been generalized in the Earth System
Modeling Framework {cite:p}`hill2004esmf` and the coupler families
YAC {cite:p}`hanke2016yac` and OASIS {cite:p}`valcke2013oasis`. These
frameworks establish the pattern of area-weighted conservative
overlap that this narrative adapts; they do not, however, address
the additional constraints specific to optical-property transfer,
where single-scattering albedo and asymmetry parameter must remain
inside their physical bounds after remapping. That specialization is
developed in [](methods/mesh-mapping.md).

## Contribution

This narrative describes:

- A hybrid 1D/3D radiative-transfer architecture in which an
  unstructured-mesh host (MPAS-A) dispatches different cells to
  different solvers (1D RRTMG, 3D SHDOM) within a single radiation
  step, with a clean treatment of the interior interface where the
  3D and 1D regions meet.
- A conservative Voronoi↔Cartesian mapping for optical properties
  and fluxes that preserves single-scattering albedo and
  asymmetry-parameter bounds by construction.
- An add-only set of patches to the bundled RRTMG shortwave code
  that expose per-band optical properties and per-band fluxes
  without modifying any existing upstream routine.
- A Phase 1 realization on the idealized supercell test case, in
  which the complete coupling pipeline — from optical-property
  handoff through heterogeneous MPI redistribution to diagnostic
  output — is exercised, and the mesh-mapping error floor is
  established as the baseline for the multi-phase roadmap.
