(app-exterior-calculus)=
# Exterior-calculus formulation

This appendix recasts the continuous and discrete problems in the
language of differential forms and discrete exterior calculus
(DEC). The main text and the earlier appendices treat $\varphi$
as a scalar function on a Euclidean domain,
$\bE = -\nabla\varphi$ as a vector field, and the discrete
operator as a matrix assembled from mesh metrics. The
exterior-calculus perspective is equivalent, coordinate-free, and
makes several features of the TRiSK construction transparent at
the level of *structure* rather than *coefficients*:

- The continuous Poisson problem is a single equation on the
  de Rham complex:
  $\delta(\varepsilon\, \mathsf{d}\varphi) = \rho\,\mu$, where
  $\mu$ is the volume form.
- The TRiSK div-of-grad discretization is the composition
  $\mathsf{d}^{\top}\mathsf{H}\,\mathsf{d}$ with $\mathsf{H}$ a
  *diagonal* Hodge-star matrix.
- $\mathsf{H}$ is diagonal precisely because the MPAS-A mesh is
  a *circumcentric Voronoi–Delaunay pair*: primal edges are
  orthogonal to their corresponding dual edges. This structural
  property is what makes the TRiSK operators mimetic.
- Symmetry and positive-definiteness of the assembled operator
  follow from the factorization
  $\mathsf{d}^{\top}\mathsf{H}\,\mathsf{d}$ with $\mathsf{H}$
  SPD, without any coefficient-level calculation.
- Discrete Stokes' theorem, $\mathsf{d}^{2} = 0$ at the chain
  level, and integration-by-parts all hold as *identities*, not
  approximations.

The material below is self-contained but assumes comfort with
the vector-calculus / divergence-theorem derivations of
[](horizontal-operator.md)–[](vertical-operator.md).

(app-ec-continuous)=
## Continuous setup: differential forms on the atmosphere

Let $M \subset \mathbb{R}^{3}$ denote the atmospheric domain, a
3D oriented Riemannian manifold with boundary
$\partial M = \Gamma_{\mathrm{g}} \cup \Gamma_{\mathrm{top}} \cup \Gamma_{\mathrm{lat}}$.
Write $\Omega^{k}(M)$ for the space of smooth $k$-forms on $M$.
The de Rham complex is the sequence

$$
  \Omega^{0}(M) \xrightarrow{\mathsf{d}} \Omega^{1}(M)
    \xrightarrow{\mathsf{d}} \Omega^{2}(M)
    \xrightarrow{\mathsf{d}} \Omega^{3}(M),
$$ (eq-ec-derham)

where $\mathsf{d}$ is the exterior derivative, satisfying
$\mathsf{d}^{2} = 0$. Under the Euclidean metric, $\Omega^{0}$ is
identified with scalar fields, $\Omega^{1}$ with vector fields
(via the musical isomorphism $\flat$), $\Omega^{2}$ with vector
fields (via $\flat$ after Hodge-duality), and $\Omega^{3}$ with
scalar fields (as volume densities). Under these identifications:

$$
\begin{aligned}
  \mathsf{d}: \Omega^{0} \to \Omega^{1} & \;\leftrightarrow\; \nabla, \\
  \mathsf{d}: \Omega^{1} \to \Omega^{2} & \;\leftrightarrow\; \nabla\times, \\
  \mathsf{d}: \Omega^{2} \to \Omega^{3} & \;\leftrightarrow\; \nabla\cdot\,.
\end{aligned}
$$

The Hodge star $\star : \Omega^{k} \to \Omega^{n-k}$ (with
$n = 3$) uses the metric and the volume form $\mu \in \Omega^{n}$
to identify forms with their duals. The codifferential
$\delta : \Omega^{k} \to \Omega^{k-1}$ is defined by

$$
  \delta = (-1)^{n(k+1) + 1}\, \star \mathsf{d} \star.
$$ (eq-ec-codiff-def)

On $n = 3$ for $k = 1$, $\delta = -\star \mathsf{d}\star$, which
under the form-to-vector identification corresponds to
$-\nabla\cdot$ on vector fields.

(app-ec-poisson-continuous)=
## Continuous Poisson as a de Rham equation

Poisson's equation
$\nabla \cdot (\varepsilon \nabla \varphi) = -\rho$ becomes,
under the identifications above:

$$
  \delta\bigl(\varepsilon \, \mathsf{d}\varphi\bigr) \;=\; \rho,
  \qquad \varphi \in \Omega^{0}(M).
$$ (eq-ec-poisson-forms)

Equivalently, multiplying by $\mu$ to view $\rho$ as a 3-form
$\rho\,\mu \in \Omega^{3}$,

$$
  \mathsf{d}\star(\varepsilon\,\mathsf{d}\varphi) \;=\; -\rho\,\mu.
$$ (eq-ec-poisson-forms-3)

This is the form in which the equation will be discretized
below: $\mathsf{d}\varphi$ is a 1-form,
$\star(\varepsilon\,\mathsf{d}\varphi)$ is a 2-form (representing
a flux), and $\mathsf{d}\star(\varepsilon\,\mathsf{d}\varphi)$ is
a 3-form (a volume density) that must equal $-\rho\,\mu$.

## Boundary conditions in form language

The Phase 1 boundary conditions translate as:

$$
\begin{aligned}
  \varphi|_{\Gamma_{\mathrm{g}}} &= 0
    \qquad \text{(0-form Dirichlet trace)}, \\
  \iota^{*}_{\Gamma_{\mathrm{top}}}\!\star(\varepsilon\,\mathsf{d}\varphi) &= 0
    \qquad \text{(2-form flux trace; Neumann)},
\end{aligned}
$$

where $\iota^{*}_{\Gamma}$ denotes pullback of a form to the
boundary surface $\Gamma$. The top-face Neumann condition is a
flux-trace vanishing condition on the 2-form
$\star(\varepsilon\,\mathsf{d}\varphi)$, which is the natural
(Neumann-type) boundary datum for the second-order elliptic
operator in {eq}`eq-ec-poisson-forms`.

(app-ec-discrete-setup)=
## Discrete setup: the cochain complex on MPAS-A's mesh

The MPAS-A domain is tiled by a 3D cell complex $\mathcal{K}$
built as the extrusion of a planar Voronoi tessellation in the
horizontal by a stack of $K$ vertical layers. The primal cell
complex $\mathcal{K}$ has:

| Cell | Description |
|---|---|
| Primal 0-cells | vertices of the extruded mesh |
| Primal 1-cells | edges (vertical edges + horizontal Voronoi edges) |
| Primal 2-cells | faces (lateral faces of the prism + horizontal caps) |
| Primal 3-cells | prisms (Voronoi cell $\times$ layer) |

The *dual* complex $\star\mathcal{K}$ is defined by placing a
dual 0-cell at the circumcentre of each primal 3-cell (the
cell-column centre) and taking the dual 1-cells to be segments
between adjacent dual 0-cells, passing through the centre of the
corresponding primal 2-face. Dual 2-cells are polygons orthogonal
to primal 1-edges, and dual 3-cells are dual volumes around
primal 0-vertices.

The scalar potential $\varphi$ is carried at dual 0-cells (cell
centres): a *dual 0-cochain*. The integer dimension

$$
  N \;=\; n_{\mathrm{cells}} \cdot K
$$

equals the number of dual 0-cells, i.e. the number of unknowns.

**Cochain spaces and pairings.** Let $C^{k}$ denote the space of
real-valued primal $k$-cochains (functions on primal $k$-cells)
and $\widetilde{C}^{k}$ the space of dual $k$-cochains. Pairings
$\langle\cdot,\cdot\rangle_k$ are integration-like sums over
cells.

**Coboundary operators.** The discrete exterior derivative at
level $k$ is the coboundary operator

$$
  \mathsf{d}_{k} : C^{k} \to C^{k+1},
  \qquad
  (\mathsf{d}_k \alpha)(\sigma)
  \;=\; \sum_{\tau \in \partial\sigma} [\sigma : \tau]\, \alpha(\tau),
$$ (eq-ec-d-discrete)

where $\sigma$ is a $(k+1)$-cell, $\partial\sigma$ is its
boundary as a formal chain of $k$-cells, and
$[\sigma : \tau] \in \{-1, +1\}$ is the orientation coefficient.
The matrix representation of $\mathsf{d}_k$ is the signed
incidence matrix of $(k+1)$-cells on $k$-cells. The fundamental
identity

$$
  \mathsf{d}_{k+1} \circ \mathsf{d}_k = 0
$$ (eq-ec-d-squared)

holds *exactly* at the chain level (no mesh-refinement error),
by the same combinatorial argument as $\partial^{2} = 0$ on
simplicial chains.

(app-ec-d-grad)=
## Discrete gradient on dual 0-cochains

The unknown $\varphi$ is a dual 0-cochain, i.e. a function
$\widetilde{\varphi} : \{\text{cell columns}\} \to \mathbb{R}$.
By duality the dual coboundary operator
$\widetilde{\mathsf{d}}_0 : \widetilde{C}^{0} \to \widetilde{C}^{1}$
acts on cell-centre values to produce values on dual 1-cells:

$$
  (\widetilde{\mathsf{d}}_0 \widetilde{\varphi})(\tilde{e})
  \;=\; \varphi_j - \varphi_i,
$$ (eq-ec-dual-d0)

where $\tilde{e}$ is the dual 1-cell connecting cell centres $i$
and $j$, oriented from $i$ to $j$.

This is the discrete analogue of $\mathsf{d}\varphi$ on a 0-form.
It carries only the *combinatorial* content of the gradient; the
geometric content (how the gradient is measured on the ground)
enters only through the Hodge star.

A primal 1-edge $e$ pierces exactly one dual 2-face, which is in
bijection with one dual 1-edge $\tilde{e}$. Under this bijection
the primal edge length $\ell_e$ and the dual edge length
$\tilde{\ell}_e = d_e$ are the two relevant geometric quantities
associated with edge $e$.

(app-ec-hodge)=
## Discrete Hodge star on MPAS-A's circumcentric mesh

The discrete Hodge star
$\mathsf{H}_k : C^{k} \to \widetilde{C}^{n-k}$ maps primal
$k$-cochains to dual $(n-k)$-cochains via local integration over
dual cells. The central structural fact of DEC on orthogonal
primal–dual meshes is:

:::{admonition} Diagonal Hodge star theorem
:class: note

On a mesh where every primal $k$-cell $\sigma$ is *orthogonal*
to its dual $(n-k)$-cell $\star\sigma$ (in the sense that every
primal edge is orthogonal to the corresponding dual face, and
every primal face is orthogonal to the corresponding dual edge),
the discrete Hodge star $\mathsf{H}_k$ is diagonal with respect
to the natural cochain bases, with entries
$(\mathsf{H}_k)_{\sigma\sigma} = |\star\sigma|/|\sigma|$.
:::

MPAS-A's extruded Voronoi mesh satisfies this orthogonality
*exactly*:

- In the horizontal, the primal 1-edge $e$ is a Voronoi face
  edge, perpendicular by construction to the line segment
  between cell centres (the dual 1-edge through $e$).
- In the vertical, primal edges are vertical by assumption (flat
  terrain), and layer faces are horizontal; orthogonality is
  automatic.

We need specifically $\mathsf{H}_1$, which maps dual 1-cochains
to primal 2-cochains. Its diagonal entries at edge $\tilde{e}$
(dual) / $e$ (primal, 2-face dual to $\tilde{e}$) are:

- Horizontal lateral edge $e$ (Voronoi face, length $\ell_e$,
  layer thickness $\Delta z_k$): dual 2-face area is
  $\ell_e \cdot \Delta z_k$; dual 1-edge length is $d_e$.

  $$
    (\mathsf{H}_1)_{ee}^{(\text{horiz, layer } k)}
    \;=\; \frac{\ell_e \, \Delta z_k}{d_e}.
  $$ (eq-ec-H1-horiz)

- Vertical edge at the upper face of layer $k$ (horizontal
  section, area $A_i$): dual 2-face area is $A_i$, dual 1-edge
  length is $\Delta z_{k+\hlf}$.

  $$
    (\mathsf{H}_1)_{ee}^{(\text{vert, upper})}
    \;=\; \frac{A_i}{\Delta z_{k+\hlf}}.
  $$ (eq-ec-H1-vert-up)

Analogously for the lower vertical face. These diagonal entries
are *exactly* the mesh-dependent weights `h_weight(k, e)/ε` and
`v_weight_upper(k, i)/ε` of the implementation (see
[](implementation-note.md)). The correspondence is not a
coincidence — the TRiSK construction is exactly the DEC
discretization of $\delta\!\circ\!(\varepsilon\,\mathsf{d}\cdot)$
on the circumcentric Voronoi–Delaunay pair.

(app-ec-d-star-d)=
## TRiSK div-of-grad as $\mathsf{d}^{\top}\mathsf{H}\,\mathsf{d}$

Applying the discrete analogues of {eq}`eq-ec-poisson-forms` to
the cell-centred potential $\widetilde{\varphi}$:

$$
  (-\mathbb{L})\,\widetilde{\varphi}
  \;=\; \varepsilon \cdot \mathsf{d}_0^{\top} \, \mathsf{H}_1 \,
                   \mathsf{d}_0 \, \widetilde{\varphi},
$$ (eq-ec-L-as-dTHd)

where $\mathsf{d}_0$ is {eq}`eq-ec-dual-d0` (the discrete dual
gradient), $\mathsf{H}_1$ is the diagonal Hodge star
{eq}`eq-ec-H1-horiz`–{eq}`eq-ec-H1-vert-up`, and
$\mathsf{d}_0^{\top}$ is its matrix transpose (which represents
the adjoint coboundary, i.e. the discrete divergence). The
permittivity $\varepsilon$ factors out as an overall scalar
because it is constant in Phase 1; for variable $\varepsilon$ it
is absorbed into $\mathsf{H}_1$ as an edge-dependent factor.

**Unwinding the composition row by row.** For a cell $i$, the
row $\widetilde{\varphi} \mapsto
(\mathsf{d}_0^{\top} \mathsf{H}_1 \mathsf{d}_0 \widetilde{\varphi})_i$
expands as

$$
\begin{aligned}
  (\mathsf{d}_0^{\top} \mathsf{H}_1 \mathsf{d}_0
   \widetilde{\varphi})_i
  &\;=\; \sum_{\tilde{e} \in \partial^{\top}(i)}
         [\,i : \tilde{e}\,] \,
         (\mathsf{H}_1)_{\tilde{e}\tilde{e}} \,
         (\mathsf{d}_0 \widetilde{\varphi})(\tilde{e}) \\
  &\;=\; \sum_{\tilde{e} : e \in \partial i}
         [\,i : \tilde{e}\,] \,
         \frac{\ell_e \Delta z_k}{d_e} \,
         (\varphi_j - \varphi_i) \\
  &\;=\; \sum_{e \in \partial i}
         \frac{\ell_e \Delta z_k}{d_e} \,
         (\varphi_i - \varphi_j)
         \qquad \text{(sign convention absorbed)}
\end{aligned}
$$ (eq-ec-row-expand)

plus the analogous vertical-face terms. Comparing to the
volume-integrated stencil derived in
[](horizontal-operator.md), Eq. {eq}`eq-A-horiz-sum`, the two
expressions agree identically.

Thus:

$$
  A \;=\; -\mathbb{L} \;=\; \varepsilon \cdot \mathsf{d}_0^{\top}
                             \mathsf{H}_1 \mathsf{d}_0.
$$ (eq-ec-A-equals)

The ground-Dirichlet boundary condition is imposed by removing
the ground face from the primal 2-chain in the definition of
$\mathsf{d}_0$ and adding the ground-face Hodge-star contribution
to the diagonal of $A$, equivalent to the half-distance
ground-face entry of `v_weight_lower(1, i)`.

(app-ec-spd)=
## Structural symmetry and positive-definiteness

From the factorization
$A = \varepsilon\, \mathsf{d}^{\top} \mathsf{H}\,\mathsf{d}$:

**Symmetry.** Transposition gives
$A^{\top} = \varepsilon\, \mathsf{d}^{\top} \mathsf{H}^{\top} \mathsf{d}$.
Since $\mathsf{H}$ is diagonal, $\mathsf{H}^{\top} = \mathsf{H}$,
so $A^{\top} = A$.

**Positive semi-definiteness.** For any $\vu$,

$$
  \vu^{\top} A \vu
  \;=\; \varepsilon \, \vu^{\top} \mathsf{d}^{\top} \mathsf{H}\,\mathsf{d}\, \vu
  \;=\; \varepsilon \, (\mathsf{d}\vu)^{\top} \mathsf{H} (\mathsf{d}\vu)
  \;\ge\; 0,
$$ (eq-ec-spsd)

since $\mathsf{H}$ is diagonal with positive entries. This is
the DEC form of Dirichlet's principle.

**Strict positive-definiteness under ground Dirichlet.** Equality
in {eq}`eq-ec-spsd` holds iff $\mathsf{d}\vu = 0$, i.e. iff
$\vu$ is a discrete constant. The ground-Dirichlet condition
removes the constant mode from the admissible space (any vector
with $v_{i, 1} \ne 0$ has a nonzero ground-face contribution to
$\vu^{\top} A \vu$), so $\vu^{\top} A \vu > 0$ for nonzero
admissible $\vu$. Hence $A$ is SPD.

The pair-by-pair entry comparison of [](symmetry.md) is therefore
subsumed by a two-line structural argument in DEC language.

(app-ec-mimetic)=
## Mimetic identities

Several identities that are approximations in coefficient-level
calculations hold exactly at the cochain level:

**Curl–grad vanishing.** $\mathsf{d}_1 \mathsf{d}_0 = 0$
combinatorially. Consequence for MPAS-A: the discrete curl of a
gradient is exactly zero, not $\mathcal{O}(h^{2})$.

**Div–curl vanishing.** $\mathsf{d}_2 \mathsf{d}_1 = 0$
combinatorially. Consequence: the discrete divergence of a curl
vanishes exactly.

**Discrete Stokes' theorem.** For any primal $k$-chain $c$ and
primal $k$-cochain $\alpha$,

$$
  \langle \mathsf{d}\alpha, c \rangle
  \;=\; \langle \alpha, \partial c \rangle
$$ (eq-ec-stokes)

exactly, by the definition of the coboundary $\mathsf{d}$ as the
dual of the boundary $\partial$. Applied to the divergence
stencil, this is discrete integration by parts: the sum of
fluxes out of a cell equals minus the sum of fluxes into
neighbouring cells, with the shared edge contributing equally to
both sides. This was verified coefficient-level in
{ref}`app-horizontal-consistency`; the DEC statement makes it
structural.

**Integration-by-parts for the Laplacian.** For admissible
$\vu$, $\vv$ (satisfying the Dirichlet BC),

$$
  \vu^{\top} A \vv
  \;=\; \varepsilon\, (\mathsf{d}\vu)^{\top} \mathsf{H} (\mathsf{d}\vv)
  \;=\; \varepsilon\, (\mathsf{d}\vv)^{\top} \mathsf{H} (\mathsf{d}\vu)
  \;=\; \vv^{\top} A \vu,
$$ (eq-ec-ibp)

which is the discrete form of Green's first identity.

(app-ec-convergence)=
## Convergence from the DEC perspective

The structural advantage of the DEC discretization is that
*combinatorial* identities like $\mathsf{d}^{2} = 0$ and
adjointness hold exactly, so the discretization error is
concentrated entirely in the *Hodge star*. Intuitively: the
mesh's combinatorial content is always right, and the geometric
content is right to the extent that $\mathsf{H}$ approximates the
continuous Hodge.

On a circumcentric Voronoi–Delaunay pair with quasi-uniform cell
sizes, the diagonal Hodge star {eq}`eq-ec-H1-horiz` is
second-order accurate as an approximation to the continuous
$\star$. Combined with first-order accuracy of the discrete
gradient {eq}`eq-ec-dual-d0`, the composition
$\mathsf{d}^{\top}\mathsf{H}\mathsf{d}$ yields an
$\mathcal{O}(h^{2})$ approximation of the continuous Laplacian in
$L^{2}$ under mesh-regularity assumptions
{cite:p}`hirani2003discrete`.

On distorted or variable-resolution meshes, the quality of the
diagonal Hodge-star approximation degrades. Specifically, the
DEC framework suggests that convergence-rate losses on distorted
meshes reflect the deviation of the discrete Hodge from the
continuous one, which is dominated by the deviation from perfect
orthogonality ("pathology" of non-centroidal Voronoi cells).
This prediction motivates the Voronoi-centroid-offset correlation
analysis in
[](../../verification-plan/tier-c-mesh.md) (Tier C.3).

(app-ec-feec)=
## Connection to finite element exterior calculus

Finite element exterior calculus (FEEC) is the
Galerkin/variational counterpart of DEC {cite:p}`arnold2006finite`.
Using Whitney forms (piecewise polynomial differential forms
compatible with the de Rham complex) on a simplicial
triangulation, FEEC produces a Galerkin discretization of
{eq}`eq-ec-poisson-forms` that is structurally similar to the
DEC formulation above.

On the Delaunay triangulation associated with MPAS-A's Voronoi
mesh, the lowest-order Whitney 0-forms are piecewise-linear nodal
elements (P1), Whitney 1-forms are lowest-order Nédélec edge
elements {cite:p}`bossavit1998computational`, and so on. The
*mass-lumped* Galerkin discretization of
$\delta(\varepsilon\,\mathsf{d}\varphi) = \rho$ using these
Whitney spaces coincides with the DEC formulation when the mass
matrices are replaced by diagonal lumped approximations, and
these lumped mass matrices coincide with the diagonal Hodge
stars of {eq}`eq-ec-H1-horiz`–{eq}`eq-ec-H1-vert-up` on
circumcentric Voronoi–Delaunay pairs.

This observation provides two complementary characterizations of
TRiSK's elliptic operator:

1. As a *DEC discretization* on the Voronoi–Delaunay complex,
   with closed-form diagonal Hodge stars derived from mesh
   metrics.
2. As a *lumped-mass Galerkin FEEC discretization* on the
   Delaunay triangulation using lowest-order Whitney elements,
   with mass lumping justified by the orthogonality of the
   circumcentric pair.

The two views are equivalent on circumcentric meshes, but
diverge when the orthogonality assumption fails. This is the
deep reason why MPAS-A meshes are required to be "well-centred"
(cell circumcentre lies inside the cell) for the operators to
behave well.

(app-ec-summary)=
## Summary

Reading the implementation
([](horizontal-operator.md), [](implementation-note.md)) in DEC
language:

| Code quantity | DEC interpretation |
|---|---|
| $\varphi$ at cell centres | dual 0-cochain $\widetilde{\varphi}$ |
| $\varphi_j - \varphi_i$ at edge | $(\mathsf{d}_0 \widetilde{\varphi})(\tilde{e})$ |
| `h_weight(k, e) / ε` | $(\mathsf{H}_1)_{\tilde{e}\tilde{e}}$ (horizontal) |
| `v_weight_upper(k, i) / ε` | $(\mathsf{H}_1)_{\tilde{e}\tilde{e}}$ (vertical) |
| Matvec $A \widetilde{\varphi}$ | $\varepsilon\,\mathsf{d}_0^{\top}\mathsf{H}_1\mathsf{d}_0\widetilde{\varphi}$ |
| Right-hand side $b = \rho \cdot V$ | $\star\rho\mu$ evaluated on dual 3-cochains |
| Ground-Dirichlet step in $\mathsf{H}_1$ | ground-face contribution in the lower primal 2-chain |
| $-\mathbb{L}$ SPD | factorization $\mathsf{d}^{\top}\mathsf{H}\mathsf{d}$ with $\mathsf{H}$ SPD |

The exterior-calculus formulation adds no computational content
that is not already in the coefficient-level implementation, but
it gives clean structural reasons for the properties proven
combinatorially in the other appendices: mesh-centric *design
choices* (e.g. insisting on circumcentric meshes, embracing the
Voronoi–Delaunay pair as the primal–dual complex) become
*consequences* of demanding a diagonal Hodge star that makes the
DEC operators computationally efficient.
