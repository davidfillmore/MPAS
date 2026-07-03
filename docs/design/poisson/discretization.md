(sec-discretization)=
# Discretization

## Mesh geometry and notation

MPAS-A discretizes the sphere as a centroidal Voronoi tessellation
{cite:p}`du1999centroidal`. Quantities carried at cell centres are
indexed by $i$, quantities at cell-face edges by $e$. The mesh
provides three primary horizontal metrics at each edge $e$:

- $\ell_e =$ `dvEdge`: the length of the Voronoi face (primal
  edge);
- $d_e =$ `dcEdge`: the distance between the two cell centres on
  either side of edge $e$;
- and, at each cell $i$, $A_i =$ `areaCell`: the cell area.

The base discretization assumes the layered-vertical coordinate is
flat ($z$-levels rather than terrain-following) so that the 3D
cell volume factors as $V_{i,k} = A_i \cdot \Delta z_k$ with
$\Delta z_k$ the thickness of layer $k$. Terrain-following meshes
can instead use the uniform-altitude `zgrid` remap mode described
in [](discussion.md), where the operator is assembled from
precomputed active-band finite-volume weights rather than the
native per-level layer thicknesses $\Delta z_k$ used below.

Layer $k$ indices run from $k = 1$ (the ground-adjacent layer) to
$k = K$ (the model-top layer). Layer-midpoint altitude is
$z_{i,k}$; layer-face altitudes are
$z^{\mathrm{f}}_{i,k+1/2}$ with thickness
$\Delta z_k = z^{\mathrm{f}}_{i,k+1/2}
           - z^{\mathrm{f}}_{i,k-1/2}$.

Each cell has a set of edge-neighbours $\mathcal{N}(i)$; for an
edge $e \in \partial i$, let $j = \mathrm{nbr}(i,e)$ denote the
neighbour cell across edge $e$.

## Continuous-to-discrete via the divergence theorem

For each cell-column $(i, k)$ we integrate
Eq. {eq}`eq-poisson` over the 3D prism $\Omega_{i,k}$ and apply
the divergence theorem:

$$
\oint_{\partial \Omega_{i,k}}
  \varepsilon \, (\nabla \varphi) \cdot \boldsymbol{n} \, dA
\;=\; -\int_{\Omega_{i,k}} \rho \, dV.
$$ (eq-int-form)

The boundary $\partial \Omega_{i,k}$ consists of the lateral
Voronoi faces (one per edge $e \in \partial i$) plus the upper and
lower horizontal faces at $z^{\mathrm{f}}_{i,k+1/2}$ and
$z^{\mathrm{f}}_{i,k-1/2}$.

## Horizontal stencil (TRiSK div-of-grad)

For each lateral face $e \in \partial i$, the flux through the
face is approximated using the cell-centre-to-cell-centre
finite-difference gradient:

$$
(\nabla \varphi) \cdot \boldsymbol{n}_e
  \;\approx\; \frac{\varphi_j - \varphi_i}{d_e},
\qquad j = \mathrm{nbr}(i, e).
$$ (eq-edge-gradient)

The face area is $\ell_e \cdot \Delta z_k$, so the contribution of
edge $e$ to the volume-integrated horizontal divergence is
$\varepsilon \, (\ell_e \Delta z_k / d_e) \,
(\varphi_j - \varphi_i)$. Summing over all edges of cell $i$:

$$
\sum_{e \in \partial i}
  \varepsilon \, \frac{\ell_e \Delta z_k}{d_e}
  \, (\varphi_j - \varphi_i).
$$ (eq-horiz-stencil)

## Vertical stencil

The upper-face flux is approximated with a second-order centred
difference:

$$
\varepsilon \, A_i \,
  \frac{\varphi_{i,k+1} - \varphi_{i,k}}{\Delta z_{k+1/2}},
$$ (eq-vert-upper)

where $\Delta z_{k+1/2} = z_{i,k+1} - z_{i,k}$ is the distance
between layer midpoints. The lower-face flux is analogous with
$\Delta z_{k-1/2}$. Boundary treatments are:

| Boundary | Condition | Discrete effect |
|---|---|---|
| $k = 1$ lower face | $\varphi = 0$ (Dirichlet) | Flux $= \varepsilon A_i (0 - \varphi_{i,1}) / (\tfrac{1}{2} \Delta z_1)$, split into diagonal and right-hand side |
| $k = K$ upper face | $\partial_z \varphi = 0$ (Neumann) | Flux set to zero |

## Assembled discrete operator

Combining the horizontal and vertical flux sums equates to the
volume-integrated form of the discrete operator at cell-column
$(i, k)$:

$$
(\mathbb{L} \, \varphi)_{i,k}
  \;=\; \sum_{e \in \partial i}
           \varepsilon \frac{\ell_e \Delta z_k}{d_e}
           (\varphi_j - \varphi_i)
        \,+\, \varepsilon A_i \,
           \frac{\varphi_{i,k+1} - \varphi_{i,k}}{\Delta z_{k+1/2}}
        \,+\, \varepsilon A_i \,
           \frac{\varphi_{i,k-1} - \varphi_{i,k}}{\Delta z_{k-1/2}},
$$ (eq-L-full)

with $j = \mathrm{nbr}(i, e)$ and with boundary faces treated per
the table above. The right-hand side of the discrete Poisson
system is the volume-integrated charge,

$$
b_{i,k} \;=\; -\rho_{i,k} \, V_{i,k}
              + \text{(ground-Dirichlet contribution)}.
$$ (eq-rhs)

The discrete system is then

$$
\mathbb{L} \, \boldsymbol{\varphi} \;=\; \boldsymbol{b}.
$$ (eq-linear-system)

## Properties of the discrete operator

The assembled matrix $\mathbb{L}$ has two essential properties.

**(P1) Sparsity.** Each row has at most `maxEdges + 2 + 1`
non-zero entries: one horizontal neighbour contribution for each
entry of `nEdgesOnCell(i)`, two vertical neighbours, and the
diagonal self-entry. On the regular hexagonal meshes used in the
Cartesian benchmark this reduces to $6 + 2 + 1 = 9$ entries.

**(P2) Symmetry and positive definiteness.** $\mathbb{L}$ is
symmetric in the standard $\ell^2$ inner product and, under the
ground-Dirichlet boundary condition, strictly negative definite.
Equivalently, $-\mathbb{L}$ is symmetric positive definite (SPD).

Property (P1) bounds matrix storage and dictates matvec cost; (P2)
permits the use of the conjugate-gradient method as the Krylov
iterative solver of choice.

## Order of accuracy

On a regular-hexagonal Voronoi tessellation with uniform vertical
spacing, the stencil {eq}`eq-L-full` has truncation error
$\mathcal{O}(h^2)$ in the bulk and $\mathcal{O}(h)$ pointwise near
the ground face, with $\mathcal{O}(h^2)$ global $L^2$ convergence
once the boundary-layer contribution is volume-weighted. On
realistic MPAS-A Voronoi meshes with quasi-uniform distortion,
supraconvergence delivers $\mathcal{O}(h^2)$ in $L^2$ under modest
regularity assumptions; on variable-resolution meshes the
transition-zone error behaviour is an empirical question that is
addressed in [](verification-plan/tier-c-mesh.md).
