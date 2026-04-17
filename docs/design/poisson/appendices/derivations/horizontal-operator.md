(app-horizontal)=
# Horizontal TRiSK div-of-grad stencil

This appendix derives the horizontal stencil
{eq}`eq-horiz-stencil` step by step starting from the continuous
Poisson equation.

(app-horizontal-cv)=
## The control volume

For a cell-column labeled by horizontal index $i$ and vertical
level $k$, the 3D control volume
$\Omega_{i, k} \subset \Rthree$ is a prism: a polygon in the
horizontal (the Voronoi cell $\omega_i \subset \mathbb{R}^{2}$)
extruded vertically from $z = z^{\mathrm{f}}_{i, k - \hlf}$ to
$z = z^{\mathrm{f}}_{i, k + \hlf}$. Under the flat-terrain
assumption, the face heights
$z^{\mathrm{f}}_{i, k \pm \hlf}$ do not depend on $i$, and the
vertical extent of the prism is $\dz{k}$.

The horizontal boundary of the prism consists of lateral faces,
one per edge $e \in \partial i$. Each lateral face is a rectangle
of horizontal length $\elen$ and vertical extent $\dz{k}$, hence
of area $\Aedge = \elen \cdot \dz{k}$.

The top and bottom of the prism are two horizontal polygons, each
of area $A_i$, normal vectors $\pm\hat{\bm{z}}$.

(app-horizontal-integrate)=
## Integrating Poisson over the prism

Integrate both sides of equation {eq}`eq-poisson` over
$\Omega_{i, k}$:

$$
  \int_{\Omega_{i, k}}
     \nabla \cdot \big[\varepsilon \nabla \varphi \big] \, dV
  \;=\;
  -\int_{\Omega_{i, k}} \rho \, dV.
$$ (eq-A-poisson-int)

Apply the divergence theorem to the left side:

$$
  \oint_{\partial \Omega_{i, k}}
    \varepsilon \, (\nabla \varphi) \cdot \bn \; dA
  \;=\; -\int_{\Omega_{i, k}} \rho \, dV.
$$ (eq-A-divthm)

This is the starting point. The left side splits into three
contributions corresponding to the three kinds of prism face:

$$
  \underbrace{
    \sum_{e \in \partial i} \int_{\mathrm{face}\,e}
        \varepsilon (\nabla \varphi) \cdot \bn_e \; dA
  }_{\mathcal{F}^{h}_{i, k}}
  \;+\;
  \underbrace{
    \int_{\mathrm{top\,face}}
        \varepsilon (\nabla \varphi) \cdot \hat{\bm{z}} \; dA
  }_{\mathcal{F}^{+}_{i, k}}
  \;-\;
  \underbrace{
    \int_{\mathrm{bot\,face}}
        \varepsilon (\nabla \varphi) \cdot \hat{\bm{z}} \; dA
  }_{\mathcal{F}^{-}_{i, k}}
  \;=\; -\int_{\Omega_{i, k}} \rho \, dV.
$$ (eq-A-split)

The minus sign on the bottom-face term arises because the outward
normal points in the $-\hat{\bm{z}}$ direction: the outward flux
is $-(\nabla \varphi) \cdot \hat{\bm{z}}$.

(app-horizontal-face)=
## Approximating one horizontal face flux

Consider a single lateral face at edge $e$. The face is a
rectangle with outward normal $\bn_e$ pointing from cell $i$ to
cell $j = \mathrm{nbr}(i, e)$. Its area is
$\elen \cdot \dz{k}$. On this face,

$$
  \int_{\mathrm{face}\,e}
    \varepsilon (\nabla \varphi) \cdot \bn_e \; dA
  \;\approx\;
  \elen \cdot \dz{k} \cdot \varepsilon_e \cdot
  \big[ (\nabla \varphi) \cdot \bn_e \big]_{\bar{\bx}_e, z_k},
$$ (eq-A-face-flux-raw)

where the bracket on the right is the edge-normal gradient
evaluated at the face midpoint $(\bar{\bx}_e, z_k)$ with
$\bar{\bx}_e$ the horizontal edge midpoint and $z_k$ the
layer-midpoint altitude. This is a one-point midpoint quadrature
over the face: second-order accurate in both horizontal and
vertical direction because we sample at the face centroid
{cite:p}`leveque2007finite`.

The edge-normal gradient at the face midpoint is approximated by
the finite difference between the cell-center values,

$$
  \big[(\nabla \varphi) \cdot \bn_e \big]_{\bar{\bx}_e, z_k}
  \;\approx\;
  \frac{\varphi(\bm{x}_j, z_k) - \varphi(\bm{x}_i, z_k)}{\dcen}
  \;\approx\;
  \frac{\varphi_{j, k} - \varphi_{i, k}}{\dcen},
$$ (eq-A-face-grad)

where $\bm{x}_i, \bm{x}_j$ are the cell centres on either side.
The final approximation replaces the continuous $\varphi$ at the
cell centres by the discrete unknowns
$\varphi_{i, k}, \varphi_{j, k}$. The truncation error of
{eq}`eq-A-face-grad` is discussed in [](taylor.md).

Substituting {eq}`eq-A-face-grad` into {eq}`eq-A-face-flux-raw`,

$$
  \int_{\mathrm{face}\,e}
    \varepsilon (\nabla \varphi) \cdot \bn_e \; dA
  \;\approx\;
  \varepsilon_e \frac{\elen \, \dz{k}}{\dcen}
  \, (\varphi_{j, k} - \varphi_{i, k}).
$$ (eq-A-face-flux-final)

This is the single-face discrete flux. Summing over all edges of
cell $i$ produces the horizontal-flux term in {eq}`eq-A-split`:

$$
  \mathcal{F}^{h}_{i, k}
  \;\approx\;
  \sum_{e \in \partial i}
    \varepsilon_e \frac{\elen \, \dz{k}}{\dcen}
    \, (\varphi_{j(e), k} - \varphi_{i, k}),
  \qquad j(e) = \mathrm{nbr}(i, e).
$$ (eq-A-horiz-sum)

(app-horizontal-consistency)=
## Consistency of adjacent-cell face fluxes

We verify that a face shared between cells $i$ and $j$
contributes with equal magnitude and opposite sign to the two
cells' volume integrals. This is the discrete analogue of
Newton's third law; it is what makes the discrete operator
symmetric (see [](symmetry.md)).

Let edge $e$ separate cells $i$ and $j$. The horizontal edge
length $\elen$, the face thickness $\dz{k}$, the cell-centre
distance $\dcen$, and the edge permittivity $\varepsilon_e$ are
all properties of the face itself and do not depend on whether we
view the face from cell $i$ or cell $j$.

Viewing from cell $i$ with outward normal $\bn_e$ pointing toward
$j$: the flux contribution into $i$'s right-hand side is (from
{eq}`eq-A-face-flux-final`):

$$
  +\,\varepsilon_e \frac{\elen \, \dz{k}}{\dcen} \, (\varphi_{j, k} - \varphi_{i, k}).
$$

Viewing from cell $j$ with outward normal $-\bn_e$ pointing
toward $i$: the edge-normal gradient changes sign,

$$
  (\nabla \varphi) \cdot (-\bn_e)
  = -\frac{\varphi_{j, k} - \varphi_{i, k}}{\dcen}
  = \frac{\varphi_{i, k} - \varphi_{j, k}}{\dcen},
$$

so the flux into $j$'s right-hand side is

$$
  \varepsilon_e \frac{\elen \, \dz{k}}{\dcen} \, (\varphi_{i, k} - \varphi_{j, k})
  \;=\;
  -\, \varepsilon_e \frac{\elen \, \dz{k}}{\dcen}\, (\varphi_{j, k} - \varphi_{i, k}).
$$

The two contributions are equal and opposite. This antisymmetry
at the face level is the discrete statement of flux conservation
and is the source of the $\ell^{2}$ symmetry of $\mathbb{L}$.

(app-horizontal-row)=
## Constructing the horizontal row of $\mathbb{L}$

Equation {eq}`eq-A-horiz-sum` expands to

$$
\begin{aligned}
  \mathcal{F}^{h}_{i, k}
  &\;\approx\; \sum_{e \in \partial i}
      \varepsilon_e \frac{\elen \, \dz{k}}{\dcen} \, \varphi_{j(e), k}
      \;-\; \sum_{e \in \partial i}
      \varepsilon_e \frac{\elen \, \dz{k}}{\dcen} \, \varphi_{i, k} \\
  &\;=\; \sum_{e \in \partial i}
      \varepsilon_e \frac{\elen \, \dz{k}}{\dcen} \, \varphi_{j(e), k}
      \;-\; \varphi_{i, k} \sum_{e \in \partial i}
      \varepsilon_e \frac{\elen \, \dz{k}}{\dcen}.
\end{aligned}
$$

So the horizontal contribution to row $(i, k)$ of the matrix
$\mathbb{L}$ is:

- off-diagonal: for each $e \in \partial i$, the entry at column
  $(j(e), k)$ is $+\,\varepsilon_e \, \elen \dz{k} / \dcen$;
- diagonal: the entry at column $(i, k)$ picks up a contribution
  of $-\sum_{e \in \partial i} \varepsilon_e \elen \dz{k}/\dcen$.

(app-horizontal-matrix)=
## From flux to matrix equation

Combining {eq}`eq-A-split` and {eq}`eq-A-horiz-sum` with the
vertical fluxes (derived in [](vertical-operator.md)), and
applying the midpoint quadrature to the right-hand side volume
integral,

$$
  \int_{\Omega_{i, k}} \rho \, dV
  \;\approx\;
  \rho_{i, k} \cdot \vol,
  \qquad \vol = A_i \, \dz{k},
$$ (eq-A-rhs-int)

we obtain the discrete equation for cell-column $(i, k)$:

$$
  \mathcal{F}^{h}_{i, k} + \mathcal{F}^{+}_{i, k} - \mathcal{F}^{-}_{i, k}
  \;=\; -\rho_{i, k} \, \vol.
$$ (eq-A-discrete-eq)

The linear-system form is $\mathbb{L} \bphi = \bb$ with
right-hand side $b_{i, k} = -\rho_{i, k} \vol$ (plus BC
contributions to be added from [](vertical-operator.md)). The
horizontal contribution to the row is exactly equation
{eq}`eq-horiz-stencil`.
