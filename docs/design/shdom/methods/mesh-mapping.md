(sec-mesh-mapping)=
# Structured–unstructured mesh mapping

This is the core numerical contribution.

## Conservative area-weighted remapping

Let $\{c_i\}$ denote the Voronoi cells of MPAS-A and $\{q_j\}$ the
Cartesian cells of SHDOM. For each overlapping pair $(c_i, q_j)$,
the horizontal polygon intersection area

$$
A(i,j) \;=\; \mathrm{area}\bigl(c_i \cap q_j\bigr)
$$

defines the weight with which the value of an intensive field on
cell $c_i$ contributes to the remapped value on cell $q_j$, and
conversely for the reverse map. The area-weighted mean of a field
$\phi$ from Voronoi to Cartesian cells is

$$
\phi^{\,q}_j \;=\;
    \frac{\displaystyle\sum_{i} A(i,j)\,\phi^{\,c}_i}
         {\displaystyle\sum_{i} A(i,j)},
\qquad
\phi^{\,c}_i \;=\;
    \frac{\displaystyle\sum_{j} A(i,j)\,\phi^{\,q}_j}
         {\displaystyle\sum_{j} A(i,j)},
$$ (eq-v2c-c2v)

and extensive fields (fluxes integrated over area, flux divergences
integrated over volume) are remapped by the analogous area- or
volume-weighted sums. Mass consistency is exact: because $\sum_j
A(i,j) = \mathrm{area}(c_i)$ and $\sum_i A(i,j) =
\mathrm{area}(q_j)$ by the definition of $A$, the integral of any
intensive field over the overlapping region is preserved to
rounding under both directions of the map.

The central approximation is that horizontal overlap is the sole
determinant of the remapping weight. The Phase 1 configuration is
flat (the supercell test case has no terrain), and the Phase 1
vertical discretization is identical between MPAS-A and SHDOM (the
same layer thicknesses, the same number of levels); no vertical
remapping is required. Phase 4 will extend the mapper to
terrain-following vertical coordinates on MPAS-A's
variable-resolution meshes, a strict extension of the present
algorithm that adds a second remapping stage without altering the
horizontal weights.

Two specializations keep the Phase 1 implementation compact. The
MPAS-A supercell mesh is a regular hexagonal tiling, oriented
axis-aligned with respect to the SHDOM Cartesian grid, so each
$A(i,j)$ is a hexagon-rectangle polygon intersection with a
closed-form area in terms of the relative offset of cell centres
and the two cell sizes. The weights depend only on the mesh and
box configuration — not on the field being mapped — and are
computed once at framework initialization. They are stored as a
sparse CSR-like structure with one row per Cartesian cell and one
non-zero entry per overlapping Voronoi cell (typically two to
seven), at a cost negligible compared to the MPAS-A mesh itself.

Future phases will replace the hex-rectangle specialization with a
general polygon-intersection routine: MPAS-A variable-resolution
meshes do not in general admit an axis-aligned regular-hex
structure. The area-weighted formula {eq}`eq-v2c-c2v` is
independent of which polygon-intersection algorithm supplies the
weights, so this is a drop-in replacement that does not touch the
rest of the framework. Complete algebraic details of the
hex-rectangle intersection, including the edge cases near the
periodic boundary of the supercell domain, are given in
[](../appendices/mesh-mapper.md).

## Bounded, physically meaningful transfer of scattering quantities

Naively area-averaging single-scattering albedo $\omega$ and
asymmetry parameter $g$ violates radiative energy balance. We map
extinction $\beta$, scattering coefficient $\beta_s = \beta
\omega$, and $\beta_s g$ as three independent intensive
quantities, recovering $\omega$ and $g$ on the target mesh via

$$
\omega_{\mathrm{target}}
  = \frac{\beta_{s,\mathrm{target}}}{\beta_{\mathrm{target}}},
\qquad
g_{\mathrm{target}}
  = \frac{(\beta_s g)_{\mathrm{target}}}{\beta_{s,\mathrm{target}}}.
$$ (eq-omega-g-recover)

This preserves physical bounds $\beta \geq 0$, $\omega \in [0,1]$,
$g \in [-1,1]$ by construction.

## Bidirectional mapping and round-trip error

The mapper is not an isometry. The composition $\mathcal{M}_{c\to
v} \circ \mathcal{M}_{v\to c}$ applied to a field defined on the
Voronoi mesh returns a field on the same Voronoi mesh, but one
that has been smoothed by the successive area-weighted averages —
the round-trip error. Characterizing this error as a function of
the ratio of Cartesian to Voronoi cell sizes, of the smoothness of
the mapped field, and of the location within the domain (interior
vs. near the periodic edge) is the numerical price of the hybrid
1D/3D arrangement, and a principal numerical finding of the
eventual validation work.

The characterization is established by exercising a standalone
mapper test driver on a battery of synthetic fields of known
analytic form: spatially constant fields (for which the round-trip
error must be at rounding); linear fields (for which the
round-trip error must be first-order accurate in the cell size);
and prescribed sinusoidal fields at resolved and
marginally-resolved wavelengths (for which the error attenuation
with wavenumber is measured). Alongside the round-trip test, a
conservation identity — the integral $\int \beta \, dA$ of the
mapped extinction over the full domain — is checked after each of
the forward and reverse maps and must be preserved to rounding.

The reporting form for the round-trip characterization is two
scalar summaries per field — the $L^2$ and $L^\infty$ norms of the
Voronoi-mesh difference between the original and round-tripped
field, normalized by the original's $L^2$ or $L^\infty$ norm. On
the regular hexagonal supercell mesh with its axis-aligned
Cartesian target, the round-trip error is expected to follow the
standard area-weighted remapping behaviour: at rounding for
constants, at first order in cell size for linear fields, and
decaying as a power of the ratio of Cartesian to Voronoi cell size
for resolved oscillatory fields. Later phases that map between
dissimilar mesh types (Phase 4's variable-resolution MPAS-A mesh
to uniform Cartesian, Phase 3's 3D sub-domain to 1D exterior) will
reuse this framework with modified weights; the round-trip test
procedure carries over without change.
