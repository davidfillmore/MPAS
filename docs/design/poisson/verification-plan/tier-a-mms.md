(tier-a)=
# Tier A — Verification by manufactured solutions

Tier A establishes the formal order of accuracy of the discrete
operator against two manufactured-solution problems, one in planar
Cartesian geometry and one in global spherical geometry. In each
case an analytic $\varphi_{\mathrm{exact}}(\boldsymbol{x})$ is
chosen that satisfies the boundary conditions of
[](../governing-equations.md) exactly, and the source
$\rho = -\nabla\!\cdot\!(\varepsilon\nabla\varphi_{\mathrm{exact}})$
is computed analytically and supplied to the solver. The discrete
$\varphi$ returned by the solver is then compared cell-by-cell
against $\varphi_{\mathrm{exact}}$ on a sequence of progressively
refined meshes, and the $L^2$ and $L^\infty$ error norms are
reported as functions of nominal cell size $h$.

## Tier A.1 — Cartesian sinusoid

On a doubly-periodic Voronoi nh-mesh with $K = 50$ uniform
vertical layers, the manufactured potential is

$$
\varphi_{\mathrm{exact}}(x, y, z)
  = \sin(k_x x)\,\sin(k_y y)\,\sin(k_z z),
$$ (eq-mms-cart)

with $k_x = k_y = 2\pi / L$ and $k_z = \pi / (2H)$ chosen so that
the ground Dirichlet and model-top Neumann conditions
{eq}`eq-bc-ground`–{eq}`eq-bc-top` are satisfied exactly.
Substitution into $-\nabla\!\cdot\!(\varepsilon_0 \nabla\varphi)$
yields
$\rho = \varepsilon_0 (k_x^2 + k_y^2 + k_z^2)\,\varphi_{\mathrm{exact}}$.
The mesh sequence is $h = 15, 7.5, 3.75, 1.875$ km horizontal cell
spacing, holding the vertical discretization fixed. The pass
criterion is that the log–log slope of the $L^2$ error against
$h$, computed from the two finest levels, satisfies
$\text{slope} \geq 1.9$.

## Tier A.2 — Spherical-harmonic manufactured solution

On a global quasi-uniform MPAS-A mesh, the manufactured potential
is the tensor product of a spherical harmonic in the horizontal
and a radial function in the vertical,

$$
\varphi_{\mathrm{exact}}(\theta, \lambda, r)
  = Y_{l}^{m}(\theta, \lambda) \, R_{n}(r),
$$ (eq-mms-sphere)

with $l = 4$, $m = 2$, and $R_n(r)$ chosen so that the ground and
top boundary conditions are exactly satisfied. Substitution into
the continuous operator in spherical coordinates yields the
analytic $\rho$. The mesh sequence is $h = 120, 60, 30$ km global
quasi-uniform spacing. The pass criterion is $\text{slope} \geq
1.9$ on the two finest levels; some degradation on the coarsest
mesh is anticipated from harmonic resolvability and is not
disqualifying.

### Current diagnostic status

The implemented Tier A.2 runner separates the coupled three-dimensional
manufactured-solution smoke case from the horizontal-isolated formal accuracy
case, matching the Tier A.1 split. The official SCVT mesh sequence has not yet
met the original global second-order acceptance criterion. Operator-only and
solution-error diagnostics localize the convergence loss to graph rings around
the twelve pentagonal cells in the spherical mesh; excluding a fixed number of
defect rings restores near-second-order $L^2$ behavior in the regular
hexagonal region.

A follow-on operator-only prototype tested two local correction ideas. Bounded
positive shared edge-factor retuning preserves the current symmetric two-point
operator structure, but did not restore the A.2 global convergence rate. A
local tangent-plane quadratic least-squares replacement did restore global
$L^2$ convergence in the prototype, but it leaves the current two-point
SPD/mimetic operator class and is therefore diagnostic only. The A.2 gate
remains deferred until the project either documents the current
defect-limited global behavior or adopts a deliberately designed
mimetic/multi-point defect correction.

The manufactured-solution sources {eq}`eq-mms-cart` and
{eq}`eq-mms-sphere` are obtained by direct substitution into the
continuous operator; the standard algebraic expansion is omitted.
