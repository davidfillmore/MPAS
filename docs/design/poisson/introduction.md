(sec-introduction)=
# Introduction

Thunderstorm electrification, lightning, and their associated
production of nitrogen oxides (LNO$_x$) couple convective dynamics,
cloud microphysics, and atmospheric chemistry through the electric
field $\boldsymbol{E}$. Resolving $\boldsymbol{E}$ inside a
convection-resolving model requires solving an elliptic
boundary-value problem at each physics step: the Poisson equation
for the electrostatic potential $\varphi$ arising from the
space-charge distribution $\rho$. Storm-interior electric-field
strengths routinely exceed $\sim 150$–$300$ kV m$^{-1}$,
triggering dielectric breakdown and lightning channels that
rapidly redistribute charge and emit NO$_x$
{cite:p}`macgorman1998electrical`. On larger scales, the same
elliptic operator, with permittivity $\varepsilon$ replaced by
conductivity $\sigma(\boldsymbol{x})$, describes the global
atmospheric electric circuit
{cite:p}`williams2009global,rycroft2008investigation`.

The Model for Prediction Across Scales Atmosphere (MPAS-A) is a
nonhydrostatic atmospheric model discretized on centroidal Voronoi
tessellations of the sphere with C-grid staggering
{cite:p}`skamarock2012multiscale,ringler2010unified,thuburn2009numerical`.
Its horizontal discrete operators — divergence, gradient, curl —
are built from a small set of mesh metrics (cell areas, edge
lengths, cell-centre separations) and have been extensively
characterized for hyperbolic dynamics. To the author's knowledge,
no *elliptic* operator has previously been constructed on MPAS-A
meshes, and no production atmospheric model based on unstructured
Voronoi tessellations currently supports the Poisson equation as a
first-class capability.

This narrative establishes such a capability. A discrete Laplacian
is defined on MPAS-A's unstructured-horizontal, layered-vertical
mesh by applying the TRiSK div-of-grad construction in the
horizontal and a centered-difference stencil in the vertical.
After multiplication by cell volume, the resulting operator is
sparse and symmetric positive definite (SPD), and is amenable to
standard Krylov iteration. The linear system is solved by
preconditioned conjugate gradients (PCG). Phase 1 evaluates the
method through a five-tier framework, with the first, second, and
fifth tiers reported here as completed results and the
mesh-sensitivity and scaling tiers retained as follow-on
measurements:

1. **Tier A** ([](verification-plan/tier-a-mms.md)) — formal
   verification via the method of manufactured solutions (MMS) in
   Cartesian and spherical geometries.
2. **Tier B** ([](verification-plan/tier-b-idealized.md)) —
   qualitative tests on canonical storm-like charge distributions.
3. **Tier C** ([](verification-plan/tier-c-mesh.md)) —
   characterization of convergence-rate sensitivity to mesh
   distortion and to variable-resolution mesh-refinement
   transitions.
4. **Tier D** ([](verification-plan/tier-d-scaling.md)) — parallel
   strong- and weak-scaling measurements.
5. **Tier E** ([](verification-plan/tier-e-showcase.md)) — an
   end-to-end illustration on a simulated supercell thunderstorm
   using a lightweight electrification stub.

The Phase 1 results reported here establish the in-tree operator
and solver, verify second-order horizontal convergence on the
Cartesian MMS sequence, characterize the spherical pentagon-defect
limitation, and demonstrate physically interpretable diagnostic
fields for idealized and dynamically evolving charge sources.
Tier C remains the principal future numerical-methods study. On
variable-resolution MPAS-A meshes — which are a principal
motivation for the TRiSK framework — the convergence rate of the
elliptic operator across refinement transitions is not known
a priori; Tier C specifies the measurement programme that will
characterize it directly by numerical experiment. Tier D will
measure parallel scaling and preconditioner performance.

This narrative proceeds as follows.
[](governing-equations.md) states the governing equations and
boundary conditions. [](discretization.md) defines the discrete
operator. [](solver.md) describes the solver.
[](verification-plan/index.md) specifies the five-tier
verification programme. [](discussion.md) discusses applications
and [](conclusions.md) concludes. A summary table of stencil
coefficients appears in
[](appendices/stencil-summary.md), and the extended derivation
appendices expand the mathematical details of the operator,
solver, and field reconstruction.
