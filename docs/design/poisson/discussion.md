(sec-discussion)=
# Discussion and outlook

The Poisson solver presented here provides a general-purpose elliptic
capability on MPAS-A's unstructured Voronoi mesh. The immediate
follow-on work falls into numerical-methods extensions and atmospheric
applications.

## Defect-aware spherical operator

The spherical MMS identifies the twelve pentagonal cells and their
surrounding rings of cells as the dominant error source. The operator
intentionally remains the conservative SPD two-point stencil used by
PCG. Exploratory prototyping narrows the correction space. A *locally*
enriched Hodge — an un-lumped Whitney edge-mass applied only in the
defect neighborhoods — creates a discretization interface between the
enriched and bulk regions and loses consistency there
($\mathcal{O}(1/h)$ residual in the first defect ring), so any
consistent correction must be applied *globally*. A globally consistent
cotangent (P1-equivalent) Hodge, however, coincides with the baseline
two-point weight $\ell_e/d_e$ to within $\sim\!6 \times 10^{-4}$ on the
well-centered SCVT mesh — the identity
$\tfrac{1}{2}(\cot\alpha + \cot\beta) = \ell_e/d_e$ — so it cannot raise
the order above the baseline. Closing the spherical defect therefore
appears to require a genuinely higher-order reconstruction — for example
a mass-consistent (Galerkin) right-hand side or a higher-order gradient
reconstruction — rather than a different diagonal Hodge weight, while
still preserving the conservation and symmetry the solver relies on. The
mesh-sensitivity study ([](mesh-sensitivity.md)) shows that this
conclusion is not specific to the twelve pentagons: the same
higher-order reconstruction is what restoring second-order accuracy on
the irregular, variable-resolution meshes requires, because the diagonal
two-point weight is order-limited on any non-hexagonal tessellation. A
finite-element (Galerkin) Poisson discretization on the same mesh, which
is second order independent of cell topology, is the natural route. This
remains a future numerical-methods problem.

## Lightning parameterization

Dielectric breakdown occurs when $|\boldsymbol{E}|$ exceeds a threshold
$E_{\mathrm{brk}}$ of roughly $150$–$300$ kV m$^{-1}$ at storm
altitudes, with a pressure-dependent correction. A lightning scheme
would detect breakdown cells from the Poisson solver output, construct a
flash channel by stochastic or deterministic propagation, neutralize
charge along the channel, and re-solve Poisson with the updated $\rho$.
Schemes of this type exist for regular-grid models
{cite:p}`mansell2005charge`; the MPAS-A implementation would be the
first on an unstructured spherical mesh.

## Lightning NO$_x$ (LNO$_x$) emission

Following the lightning scheme, NO$_x$ emission per flash is estimated
from flash channel geometry {cite:p}`barthe2012lightning` and coupled
into chemistry. In the DAVINCI-MPAS and MPAS-Model-ACOM-dev ecosystem
forks, this routes LNO$_x$ directly into the in-tree or MUSICA/MICM
chemistry solvers, respectively.

## Global atmospheric electric circuit

Replacing constant $\varepsilon$ with spatially variable conductivity
$\sigma(\boldsymbol{x})$ in {eq}`eq-poisson` and imposing an effective
ionospheric-top boundary condition turns the solver into a quasi-static
global electric-circuit solver. The discrete operator {eq}`eq-L-full`
remains structurally identical; only the edge coefficient is changed.
The solver scales to this use case without structural modification.

## Terrain-following vertical coordinate

MPAS-A's operational vertical coordinate is terrain-following. Rather
than add cross-derivative metric terms to the Laplacian, the solver
lifts the flat-terrain assumption of [](discretization.md) by the
altitude-grid remap of [](terrain-following.md): charge is remapped
conservatively onto a uniform vertical grid with shaved cut cells and
grounded walls, the same SPD operator is solved, and the fields are
returned to the terrain-following coordinate. The gentle-hill
manufactured solution recovers the second-order interior rate; the one
remaining approximation is that the near-surface field over steep slopes
is first order.

## Larger problems: algebraic multigrid

For very high-resolution applications ($\gtrsim 10^8$ unknowns), the
iteration-count growth of diagonal-preconditioned CG motivates replacing
the in-tree solver with a multigrid method. Algebraic multigrid via
hypre {cite:p}`falgout2006hypre` is a natural candidate; hypre's
BoomerAMG is known to deliver near-optimal scaling on elliptic operators
of this class. Integrating hypre preserves the operator developed here
and replaces only the iterative-solver layer.
