(app-mesh-mapper)=
# Mesh-mapper algorithm details

This appendix records the three implementation details of the
Phase 1 mesh mapper that are referenced, but not developed in
full, in the main text: the closed-form area of the
hexagon-rectangle overlap that provides the Phase 1 weights; the
sparse storage layout in which those weights are held; and the
treatment of the two periodic boundaries of the supercell domain.

## Hexagon-rectangle overlap area

For a regular hexagon $c_i$ of circumradius $r$ centred at
$(x_i, y_i)$, oriented with one edge parallel to the $x$-axis, and
a rectangle $q_j$ of dimensions $(\Delta x, \Delta y)$ centred at
$(x_j, y_j)$, the signed intersection area $A(i,j)$ is a piecewise
polynomial function of the offset $(x_j - x_i, y_j - y_i)$.
Because both polygons are convex and axis-aligned (up to the
hexagon's $60^\circ$ rotational symmetry), their intersection is
itself a convex polygon whose vertices are a subset of the hexagon
vertices contained in the rectangle, the rectangle vertices
contained in the hexagon, and the pairwise edge intersections.
For the parameter regimes of interest on the supercell mesh —
rectangle size comparable to hexagon size, both small relative to
the domain — the intersection polygon has between three and eight
vertices, and its area is summed contribution by contribution in
closed form. The analytic formula is implemented in
`mpas_atmphys_rt_mesh_mapper` as a pure function of the offset
vector and the two shape parameters $(r, \Delta x, \Delta y)$.

## Sparse storage

The weight table is stored as a sparse row-major structure with
one row per Cartesian cell $q_j$, holding, for each overlapping
Voronoi cell $c_i$, the pair $(i, A(i,j))$. For the regular
hex/axis-aligned rectangle configuration each row has between two
and seven non-zero entries. The same table is read in reverse for
the $c\!\to\!v$ direction without recomputation. The total
storage is $O(N_q)$ where $N_q$ is the number of Cartesian cells,
with the constant depending on the relative cell-size ratio; for
Phase 1 the table is a few megabytes and is allocated once at
framework initialization.

## Periodic boundary treatment

The supercell domain is doubly periodic on its lateral edges. For
each hexagon-rectangle overlap computation that straddles a
periodic edge, the rectangle centre is translated by the
appropriate domain-length offset and the overlap is computed in
the translated frame before being added to the sparse table. The
procedure is a straightforward extension of the interior
computation and produces no special cases beyond bookkeeping; the
alternative — using a copy of the hexagon translated into the
domain interior — produces identical weights at double the
per-overlap work.

Further details on the numerical subtleties of bound-preserving
optical-property remapping are developed in
[](../methods/mesh-mapping.md) of the main text.
