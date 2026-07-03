(app-taylor)=
# Taylor analysis and order of accuracy

This appendix demonstrates second-order accuracy on a
regular-hexagonal Voronoi tessellation with uniform vertical
spacing. We focus on the horizontal stencil; the vertical is the
standard 1D centred second-derivative stencil with well-known
$\mathcal{O}(\Delta z^{2})$ truncation.

## Setup: regular hex mesh

On a regular hexagonal mesh with edge length $h$, a cell $i$ has
six edge-neighbours located at cell-centre offsets

$$
  \bm{r}_{\theta} = h\sqrt{3} \, (\cos\theta, \sin\theta),
  \qquad \theta = 0, \pi/3, 2\pi/3, \pi, 4\pi/3, 5\pi/3.
$$

The mesh metrics are constant: $\elen = h$, $\dcen = h\sqrt{3}$,
$A_i = (3\sqrt{3}/2) h^{2}$. The vertical thickness
$\dz{k} = \Delta z$ is constant.

## Per-unit-volume form of the horizontal operator

Divide the horizontal part of equation {eq}`eq-L-full` by
$\vol = A_i \Delta z$:

$$
  \frac{1}{\vol}\, \mathcal{F}^{h}_{i, k}
  \;=\; \frac{1}{A_i \Delta z}
        \sum_{e \in \partial i} \varepsilon \frac{h \cdot \Delta z}{h\sqrt{3}}
        (\varphi_{j, k} - \varphi_{i, k})
  \;=\; \frac{\varepsilon}{A_i \sqrt{3}}
        \sum_{\theta} (\varphi_{\theta, k} - \varphi_{i, k}),
$$

where $\varphi_{\theta, k}$ is the cell-centre value at offset
$\bm{r}_\theta$. With $A_i = (3\sqrt{3}/2) h^{2}$ this simplifies
to a prefactor we will recompute carefully below:

$$
  \frac{1}{\vol}\, \mathcal{F}^{h}_{i, k}
  \;=\; \frac{\varepsilon}{A_i} \cdot \frac{1}{\sqrt{3}}
        \sum_{\theta} (\varphi_{\theta, k} - \varphi_{i, k}).
$$ (eq-D-unit-vol)

## Taylor expansion of $\varphi_\theta$

Let $\varphi(\bm{x})$ be a smooth function. Expand at
$\bm{x}_i$:

$$
  \varphi(\bm{x}_i + \bm{r}_{\theta})
  \;=\; \varphi(\bm{x}_i)
        + \bm{r}_\theta \cdot \nabla \varphi
        + \tfrac{1}{2} (\bm{r}_\theta \cdot \nabla)^{2} \varphi
        + \tfrac{1}{6} (\bm{r}_\theta \cdot \nabla)^{3} \varphi
        + \mathcal{O}(h^{4}).
$$ (eq-D-taylor)

Sum over the six neighbours:

$$
  \sum_{\theta} \varphi_{\theta}
  \;=\; 6\varphi
         + \left(\sum_{\theta} \bm{r}_\theta\right) \cdot \nabla \varphi
         + \tfrac{1}{2} \sum_{\theta} (\bm{r}_\theta \cdot \nabla)^{2} \varphi
         + \tfrac{1}{6} \sum_{\theta} (\bm{r}_\theta \cdot \nabla)^{3} \varphi
         + \mathcal{O}(h^{4}).
$$

By symmetry of the hex (rotation by $\pi/3$),
$\sum_\theta \bm{r}_\theta = 0$ and
$\sum_\theta (\bm{r}_\theta \cdot \nabla)^{3} \varphi = 0$.

**Second-order term.** With $|\bm{r}_\theta| = h\sqrt{3}$,

$$
\begin{aligned}
  \sum_{\theta} (\bm{r}_\theta \cdot \nabla)^{2} \varphi
  &\;=\; \sum_{\theta} \left(
           r_{\theta, x}^{2} \partial_x^{2}
           + 2 r_{\theta, x} r_{\theta, y} \partial_x \partial_y
           + r_{\theta, y}^{2} \partial_y^{2}
         \right) \varphi \\
  &\;=\;
    \left( \sum_{\theta} r_{\theta, x}^{2}\right) \partial_x^{2} \varphi
    + 2 \left( \sum_{\theta} r_{\theta, x} r_{\theta, y}\right) \partial_x \partial_y \varphi
    + \left( \sum_{\theta} r_{\theta, y}^{2}\right) \partial_y^{2} \varphi.
\end{aligned}
$$

By symmetry, $\sum_\theta r_{\theta, x} r_{\theta, y} = 0$, and

$$
\begin{aligned}
  \sum_{\theta} r_{\theta, x}^{2}
  &\;=\; (h\sqrt{3})^{2} \sum_{\theta} \cos^{2}\theta
  \;=\; 3 h^{2} \cdot 3 \;=\; 9 h^{2}, \\
  \sum_{\theta} r_{\theta, y}^{2} &\;=\; 9 h^{2}
  \qquad \text{(same by symmetry)}.
\end{aligned}
$$

Therefore

$$
  \sum_{\theta} (\bm{r}_\theta \cdot \nabla)^{2} \varphi
  \;=\; 9 h^{2} \left( \partial_x^{2} + \partial_y^{2} \right) \varphi
  \;=\; 9 h^{2} \, \nabla_{h}^{2} \varphi.
$$ (eq-D-second-order)

## Putting it together

Substituting into {eq}`eq-D-unit-vol`:

$$
\begin{aligned}
  \frac{1}{\vol} \mathcal{F}^{h}_{i, k}
  &\;=\; \frac{\varepsilon}{A_i} \cdot \frac{1}{\sqrt{3}}
         \left[ \sum_\theta \varphi_{\theta, k} - 6 \varphi \right] \\
  &\;=\; \frac{\varepsilon}{A_i} \cdot \frac{1}{\sqrt{3}}
         \left[ \tfrac{1}{2} \cdot 9 h^{2} \nabla_{h}^{2} \varphi
                + \mathcal{O}(h^{4}) \right] \\
  &\;=\; \frac{\varepsilon}{(3\sqrt{3}/2) h^{2}}
         \cdot \frac{9 h^{2}}{2\sqrt{3}}\,
         \nabla_{h}^{2}\varphi + \mathcal{O}(h^{2}) \\
  &\;=\; \varepsilon \, \nabla_{h}^{2} \varphi + \mathcal{O}(h^{2}).
\end{aligned}
$$ (eq-D-final)

The discrete horizontal operator, per unit volume, converges to
$\varepsilon \nabla_{h}^{2} \varphi$ with second-order accuracy.

## Vertical stencil truncation

The vertical stencil in {eq}`eq-B-vert-combined`, written per
unit volume:

$$
  \frac{1}{\vol} \left(\mathcal{F}^{+}_{i, k} - \mathcal{F}^{-}_{i, k}\right)
  \;=\; \frac{\varepsilon}{\dz{k}} \left[
         \frac{\varphi_{i, k+1} - \varphi_{i, k}}{\dz{k+\hlf}}
         - \frac{\varphi_{i, k} - \varphi_{i, k-1}}{\dz{k-\hlf}}
       \right],
$$

is the classical 1D centred second derivative on a non-uniform
grid. With uniform
$\dz{k} = \dz{k \pm \hlf} = \Delta z$, this reduces to

$$
  \frac{\varepsilon}{\Delta z^{2}}
  (\varphi_{i, k+1} - 2 \varphi_{i, k} + \varphi_{i, k-1})
  \;=\; \varepsilon \partial_z^{2} \varphi + \mathcal{O}(\Delta z^{2}).
$$

On a non-uniform vertical, the same stencil is first-order at
each grid point but second-order in $L^{2}$ after summation by
trapezoidal quadrature {cite:p}`leveque2007finite`, a form of
supraconvergence.

## Combined 3D convergence

Summing horizontal and vertical truncation,

$$
  \frac{1}{\vol} (\mathbb{L} \bphi)_{i, k}
  \;=\; \varepsilon \nabla^{2} \varphi_{i, k}
        + \mathcal{O}(h^{2}) + \mathcal{O}(\Delta z^{2}).
$$

The right-hand side is similarly of second-order local
truncation. The global $L^{2}$ error
$\|\bphi - \bphi_{\mathrm{exact}}\|_{L^{2}}$ then converges at the
expected rate $\mathcal{O}(\max\{h, \Delta z\}^{2})$ for smooth
$\varphi$.

On distorted Voronoi meshes, the per-row truncation picks up
terms that do not cancel via the hex symmetry argument above.
Global $L^{2}$ convergence persists via supraconvergence under
modest regularity assumptions, but the rate may be reduced. The
empirical behaviour on MPAS-A meshes is quantified in the
mesh-sensitivity verification plan
([](../../mesh-sensitivity.md)).
