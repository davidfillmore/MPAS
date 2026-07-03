(sec-mesh)=
# Mesh sensitivity

The order-of-accuracy estimate of [](discretization.md) holds on the
regular hexagonal meshes of the verification test suite, and the
spherical MMS shows that the twelve pentagonal cells of a quasi-uniform
icosahedral mesh already degrade the global spherical rate to
$L^2 = 1.38$. The mesh-sensitivity study asks the operational question:
what happens on the *variable-resolution* meshes that are a principal
motivation for MPAS-A, where the bulk of the mesh is no longer a regular
hexagonal tessellation? The answer, established here by direct numerical
experiment in the production solver, is that the second-order accuracy
is a property of the near-uniform hexagonal mesh and is lost on the
irregular tessellations that variable resolution requires.

## Mesh topology under variable resolution

The quasi-uniform SCVT meshes of the spherical MMS are icosahedral:
every cell is a hexagon except for exactly twelve pentagons ($0.5\%$ of
cells). A variable-resolution mesh cannot have this structure.
Generating a refined mesh with the standard MPAS-Tools spherical
generator (JIGSAW), using a circular $4\times$ refinement patch (coarse
$480$ km background, fine $120$ km core, smooth transition annulus),
produces a *general* centroidal-Voronoi tessellation in which
$5$–$10\%$ of cells are non-hexagonal (five-, seven-, and eight-sided),
distributed across the whole sphere and densest in the transition band.
The non-hexagonal fraction is an intrinsic consequence of the spatially
varying cell density: a control mesh of *uniform* density built on the
same fine generator grid is also a general centroidal-Voronoi
tessellation, with $231$ non-hexagonal cells rather than twelve. The
twelve-pentagon icosahedral mesh is the special case obtained only for a
globally uniform density on a coarse generator grid.

## Matched-control design

To separate the effect of mesh *irregularity* from the effect of
variable *resolution*, each variable-resolution mesh is compared not to
the icosahedral mesh but to a uniform-density control built on the same
generator grid, which is a general centroidal-Voronoi tessellation of
the same family. The spherical-harmonic MMS (horizontally isolated
source) is run in the production solver on a nested $4\times$ sequence
($480$–$120$, $240$–$60$, $120$–$30$ km) and on the matched uniform
control sequence ($480$, $240$, $120$ km), to PCG residuals below
$10^{-12}$ throughout.

## Convergence collapse

The table below reports the finest-pair solved-potential $L^2$ slopes.
Both general centroidal-Voronoi sequences — variable-resolution and
uniform-control — converge at a near-stalled rate of $\sim\!0.35$,
against the icosahedral mesh's $1.38$. The step slopes decrease toward
zero under refinement (the uniform-control interior is essentially
flat), the signature of an error floor: the operator carries an
$\mathcal{O}(1)$ truncation component that does not vanish under
refinement on an irregular mesh. The collapse is not caused by
degenerate cells. Every mesh in both sequences is well-centered — the
minimum two-point weight ratio $\ell_e/d_e$ is $\ge 0.20$ (median
$\sim\!0.5$), with no negative or vanishing weights — so the Voronoi
geometry is sound and the loss of order is a property of the
finite-volume operator on an irregular mesh, not of a broken
tessellation.

| Mesh family | $L^2$ slope |
|---|---:|
| Icosahedral SCVT (twelve pentagons) | $1.38$ |
| Uniform general-CVT (matched control) | $0.35$ |
| Variable-resolution general-CVT ($4\times$) | $0.38$ |

Finest-pair solved-potential $L^2$ convergence slope of the spherical
MMS for the icosahedral mesh and the two general centroidal-Voronoi
sequences. All solves use the horizontally isolated source with PCG
residuals below $10^{-12}$.

## Where the error lives

A per-cell decomposition of the solved-potential error on the matched
pair localizes the mechanism. On the uniform control, non-hexagonal
cells carry $3.4\times$ the area-normalized squared error of hexagonal
cells: the twelve-pentagon defect of the spherical MMS is the minimal
instance of a general penalty on every non-hexagonal cell. On the
variable-resolution mesh the refinement-transition band, where the
non-hexagonal fraction is highest ($14\%$ versus $9\%$ in the far
field), is a local error hotspot carrying $\sim\!1.4\times$ the
far-field relative error. Within fixed local-resolution bands — which
control for cell size — the per-cell error correlates with the local
cell-size gradient (Spearman rank correlation up to $+0.75$) and with
loss of well-centeredness ($+0.55$); the non-hexagonal flag alone is a
weaker, binary proxy ($+0.25$). The error tracks precisely the geometric
irregularity that variable resolution introduces.

## Verification and scope

The analysis pipeline reproduces the icosahedral spherical-MMS ground
truth to all reported digits (step $L^2$ slopes $1.73 / 1.48 / 1.38$
across the $480$–$60$ km sequence), confirming that the general-CVT
result is a property of the meshes and operator, not of the diagnostic.
The general centroidal-Voronoi meshes here are the direct output of the
standard mesh generator; they are well-centered but only moderately
optimized (median $\ell_e/d_e \approx 0.5$ against the regular-hexagon
value $1/\sqrt{3} \approx 0.58$). Whether a more aggressively optimized
variable-resolution SCVT recovers part of the lost order is an open
question, but the meshes that current variable-resolution MPAS-A
configurations actually use are general centroidal-Voronoi
tessellations of this class, so the result is operationally relevant.
For the present simple finite-volume operator the conclusion is
unambiguous: its second-order accuracy is confined to near-uniform
hexagonal meshes, and restoring order on the irregular and
variable-resolution meshes that motivate MPAS-A requires a higher-order
or finite-element reconstruction ([](discussion.md)), not a different
diagonal finite-volume weight.
