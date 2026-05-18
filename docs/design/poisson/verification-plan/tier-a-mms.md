(tier-a)=
# Tier A — Verification by manufactured solutions

Tier A establishes the manufactured-solution behavior of the
discrete operator in planar Cartesian geometry and on a global
spherical MPAS-A mesh. In each case an analytic
$\varphi_{\mathrm{exact}}(\boldsymbol{x})$ is chosen that satisfies
the boundary conditions of
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
analytic $\rho$. The accepted Phase 1 mesh sequence uses the official
MPAS-A SCVT bundles at approximately $h = 480, 240, 120,$ and
$60$ km. The A.2 run is split into a coupled three-dimensional smoke
case and a horizontal-isolated diagnostic whose right-hand side applies the
discrete vertical operator to the exact solution, matching the Tier A.1 split.

### Accepted Phase 1 interpretation

Tier A.2 is accepted for Phase 1 as a characterization gate for the current
volume-integrated, symmetric two-point spherical operator. The global
horizontal-isolated solution errors on the official SCVT sequence decrease but
do not meet the original second-order global criterion: the finest-pair slopes
are $L^2 = 1.380$ and $L^\infty = 0.253$. This is not treated as an algebraic
solver failure; residuals are below $10^{-12}$ in the accepted diagnostic run.

Operator-only and solution-error diagnostics localize the convergence loss to
graph rings around the twelve pentagonal cells in the spherical mesh. Excluding
through six to eight defect rings restores near-second-order $L^2$ behavior in
the regular hexagonal region, while the global $L^\infty$ norm remains
dominated by the first pentagon-adjacent ring. The accepted A.2 result therefore
reports both the global defect-limited norms and the defect-ring-exclusion
diagnostic. The latter is a decomposition of the observed error mechanism, not
a replacement for the global error.

The production Phase 1 operator remains the conservative SPD two-point
operator. A follow-on operator-only prototype showed that bounded positive
shared edge-factor retuning does not restore the global A.2 rate. A local
tangent-plane quadratic least-squares replacement can recover global $L^2$
convergence, but it leaves the current two-point SPD/mimetic operator class and
is diagnostic only. The planned numerical-methods upgrade is therefore a
defect-aware mimetic or multi-point correction, introduced only after a
red/green operator-residual test demonstrates restored spherical harmonic
behavior without sacrificing the conservation and symmetry properties required
by the solver.

The manufactured-solution sources {eq}`eq-mms-cart` and
{eq}`eq-mms-sphere` are obtained by direct substitution into the
continuous operator; the standard algebraic expansion is omitted.
