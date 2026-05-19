(tier-c)=
# Tier C — Mesh sensitivity (planned follow-on study)

Tier C is the planned follow-on numerical-methods study: a
measurement programme for how the TRiSK div-of-grad Laplacian's
convergence rate responds to the irregular Voronoi geometries that
arise in production MPAS-A configurations. The formal
order-of-accuracy result from Tier A holds on the regular
hexagonal meshes of its test suite; the Voronoi meshes used in
practice are quasi-uniform centroidal tessellations in the bulk
and contain refinement-transition bands on variable-resolution
configurations, and the convergence rate across such bands is not
known a priori.

Three sub-studies address this.

## Tier C.1 — Perturbed quasi-uniform mesh

Starting from the ideal $7.5$ km regular-hex Voronoi tessellation
used in Tier A.1, cell centres are displaced by an isotropic
Gaussian perturbation with standard deviation $\delta \cdot h$,
for $\delta \in \{0.10, 0.25, 0.50\}$ of the mean cell spacing,
and the Voronoi tessellation is regenerated about the perturbed
generators. The Tier A.1 Cartesian method of manufactured solutions
(MMS) is then re-run on each perturbed mesh in the same $h$ sequence,
and the log–log convergence slope is computed and plotted against
$\delta$. The dependence of measured slope on $\delta$ is the
Tier C.1 result.

## Tier C.2 — Variable-resolution global mesh

The $60$–$15$ km MPAS-A X-mesh, which is the canonical published
variable-resolution configuration, is used as a second test bed.
The Tier A.2 spherical-harmonic MMS is run on this mesh in two
placements: (i) with the refined region located diametrically away
from the harmonic's peak (control); and (ii) with the refined region
centred on the peak (sensitivity). Global and regional $L^2$ error
norms are reported for both placements, and the error contribution
localized to the refinement-transition band is extracted by difference.

## Tier C.3 — Voronoi-offset correlation

On the $\delta = 0.50$ perturbed mesh of Tier C.1, the local
truncation error of the operator is computed cell-by-cell from the
known analytic solution, and its magnitude is regressed against
the per-cell offset
$|\boldsymbol{x}^{\mathrm{Vor}}_i - \boldsymbol{x}^{\mathrm{area}}_i|$
between the Voronoi generator and the true area centroid. The
hypothesis — supported by the standard second-order finite-volume
analysis on dual meshes — is that the truncation error at a cell
is first-order in this offset; the regression reports the
measured functional dependence.

No pass/fail criterion is applied to Tier C. The reported data,
including the measured slope-vs-$\delta$ curve from C.1, the
error attribution from C.2, and the offset regression from C.3,
will constitute the numerical finding of this planned follow-on
mesh-sensitivity study.
