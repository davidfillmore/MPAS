(sec-verification)=
# Verification

The verification evaluates the discrete operator against two
method-of-manufactured-solutions (MMS) problems, one in planar
Cartesian geometry and one in global spherical geometry. In each case
an analytic $\varphi_{\mathrm{exact}}(\boldsymbol{x})$ is chosen that
satisfies the boundary conditions of [](governing-equations.md)
exactly, and the source
$\rho = -\nabla\!\cdot\!(\varepsilon\nabla\varphi_{\mathrm{exact}})$ is
computed analytically — or by a controlled discrete split that isolates
the horizontal truncation error (described below) — and supplied to the
solver. The discrete $\varphi$ returned by the solver is then compared
cell-by-cell against $\varphi_{\mathrm{exact}}$ on a sequence of
progressively refined meshes, and the $L^2$ and $L^\infty$ error norms
are reported as functions of nominal cell size $h$.

## Cartesian sinusoid

On a doubly periodic planar hexagonal Voronoi mesh with $K = 50$
uniform vertical layers, the manufactured potential is

$$
\varphi_{\mathrm{exact}}(x, y, z)
  = \sin(k_x x)\,\sin(k_y y)\,\sin(k_z z),
$$ (eq-mms-cart)

with $k_x = 2\pi/L_x$ and $k_y = 2\pi/L_y$ matched to the periodic
domain lengths and $k_z = \pi/(2H)$ chosen so that the ground Dirichlet
and model-top Neumann conditions {eq}`eq-bc-ground`–{eq}`eq-bc-top` are
satisfied exactly. Substitution into
$-\nabla\!\cdot\!(\varepsilon_0\nabla\varphi)$ yields

$$
\rho = \varepsilon_0 (k_x^2 + k_y^2 + k_z^2)\,\varphi_{\mathrm{exact}}.
$$

The mesh sequence is $h = 15,\,7.5,\,3.75,\,1.875$ km horizontal cell
spacing, holding the vertical discretization fixed. The coupled
three-dimensional mode, retained as an integration ("smoke") test,
converges but is limited by the fixed $K = 50$ vertical grid, giving a
finest two-mesh $L^2$ slope of $1.25$. A horizontally isolated source
mode applies the discrete vertical operator to
$\varphi_{\mathrm{exact}}$ inside the right-hand side, removing the
fixed vertical truncation error from the horizontal refinement
sequence. In that formal horizontal-accuracy mode, the finest two-mesh
slopes are $L^2 = 2.000$ and $L^\infty = 2.000$, with final PCG
residuals below $7 \times 10^{-9}$ on all four meshes. This confirms
second-order horizontal convergence on the regular hexagonal mesh.

## Spherical-harmonic manufactured solution

On a global quasi-uniform MPAS-A mesh, the manufactured potential is
the tensor product of a spherical harmonic in the horizontal and a
vertical profile,

$$
\varphi_{\mathrm{exact}}(\theta, \lambda, z)
  = Y_l^m(\theta, \lambda)\, R(z),
$$ (eq-mms-sphere)

with $l = 4$, $m = 2$, and $R(z)$ chosen so that the ground and top
boundary conditions are exactly satisfied. Substitution into the
thin-shell (shallow-atmosphere) form of the operator — the
Laplace–Beltrami operator evaluated at the mean sphere radius $a$, whose
eigenvalue on $Y_l^m$ is $-l(l+1)/a^2$, plus a plain vertical second
derivative — yields the analytic $\rho$. The thin-shell form is the
continuous counterpart of the discrete operator {eq}`eq-L-full`, whose
horizontal metrics are likewise evaluated at the mean sphere radius
independent of height. The mesh sequence uses the official MPAS-A SCVT
bundles at approximately $h = 480,\,240,\,120,$ and $60$ km. The coupled
three-dimensional run is retained as a smoke test; the convergence
diagnostic uses a horizontally isolated source mode analogous to the
Cartesian case.

The global horizontally isolated errors decrease, but the current
two-point spherical operator does not satisfy the original second-order
global criterion: the finest two-mesh slopes are $L^2 = 1.380$ and
$L^\infty = 0.253$, with PCG residuals below $10^{-12}$. Re-solving the
coarsest mesh at $K = 4$, $16$, and $32$ vertical levels changes the
error by under $1.3\%$, with the error identical across vertical levels
to machine precision; the horizontally isolated source therefore
removes the vertical truncation cleanly, and the rate loss is the
horizontal pentagon defect rather than a fixed-vertical-grid limitation
(in contrast to the coupled Cartesian mode at slope $1.25$).
Operator-only and solution-error diagnostics localize the rate loss to
rings of cells (in graph distance) around the twelve pentagonal cells.
Excluding the first six to eight defect rings restores near-second-order
$L^2$ behavior in the regular hexagonal bulk, while the global
$L^\infty$ norm remains dominated by the first pentagon-adjacent ring.
The spherical MMS therefore characterizes the present conservative SPD
two-point operator rather than certifying global second-order spherical
convergence. The twelve pentagons are the minimal instance of a more
general mesh sensitivity: [](mesh-sensitivity.md) shows that when the
entire tessellation is irregular, as on the variable-resolution meshes
that motivate MPAS-A, the same mechanism degrades the convergence rate
globally. A defect-aware mimetic or multi-point correction is a future
numerical-methods problem.

The manufactured-solution sources {eq}`eq-mms-cart` and
{eq}`eq-mms-sphere` are obtained by direct substitution into the
continuous operator; the standard algebraic expansion is omitted.
