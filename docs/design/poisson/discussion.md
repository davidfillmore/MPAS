(sec-discussion)=
# Discussion and outlook

The Poisson solver presented here provides a general-purpose
elliptic capability on MPAS-A's unstructured Voronoi mesh. The
immediate follow-on work falls into numerical-methods extensions
and atmospheric applications.

## Defect-aware spherical operator

Tier A.2 identifies the twelve pentagonal cells and their
neighboring graph rings as the dominant error source for the global
spherical method of manufactured solutions (MMS). The Phase 1
operator intentionally remains the conservative symmetric positive
definite (SPD) two-point stencil used by preconditioned conjugate
gradients (PCG). A Phase 2 correction should be defect-aware and
should preserve the conservation and symmetry properties needed by
the solver, likely through a mimetic multi-point or otherwise
locally enriched spherical operator rather than scalar edge-factor
retuning.

## Lightning parameterization

Dielectric breakdown occurs when $|\boldsymbol{E}|$ exceeds a
threshold $E_{\mathrm{brk}}$ of roughly $150$–$300$ kV m$^{-1}$ at
storm altitudes, with a pressure-dependent correction. A lightning
scheme would detect breakdown cells from the Poisson solver
output, construct a flash channel by stochastic or deterministic
propagation, neutralize charge along the channel, and re-solve
Poisson with the updated $\rho$. Schemes of this type exist for
regular-grid models {cite:p}`mansell2005charge`; the MPAS-A
implementation would be the first on an unstructured spherical
mesh.

## Lightning NO$_x$ (LNO$_x$) emission

Following the lightning scheme, NO$_x$ emission per flash is
estimated from flash channel geometry
{cite:p}`barthe2012lightning` and coupled into chemistry. In the
DAVINCI-MPAS and MPAS-Model-ACOM-dev ecosystem forks, this routes
LNO$_x$ directly into the in-tree or MUSICA/MICM chemistry
solvers, respectively.

## Global atmospheric electric circuit

Replacing constant $\varepsilon$ with spatially variable
conductivity $\sigma(\boldsymbol{x})$ in {eq}`eq-poisson` and
imposing an effective ionospheric-top boundary condition turns
the solver into a quasi-static global electric-circuit solver.
The discrete operator {eq}`eq-L-full` remains structurally
identical; only the edge coefficient is changed. The solver
scales to this use case without structural modification.

## Terrain-following vertical coordinate

MPAS-A's operational vertical coordinate is terrain-following.
The flat-terrain assumption of [](discretization.md) restricts
Phase 1 benchmarks to meshes with
$z_{\mathrm{g}}(\boldsymbol{x}_h) = 0$. The generalization to
terrain-following requires the addition of cross-derivative metric
terms to the Laplacian; this is a straightforward but nontrivial
extension planned for Phase 2.

## Larger problems: algebraic multigrid

For very high-resolution applications ($\gtrsim 10^{8}$
unknowns), the iteration-count growth of diagonal-preconditioned
conjugate gradient (CG) motivates replacing the in-tree solver with
a multigrid method. Algebraic multigrid via hypre
{cite:p}`falgout2006hypre`
is a natural candidate; hypre's BoomerAMG is known to deliver
near-optimal scaling on elliptic operators of this class.
Integrating hypre preserves the operator developed here and
replaces only the iterative-solver layer.
