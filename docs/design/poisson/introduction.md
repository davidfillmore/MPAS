(sec-introduction)=
# Introduction

Thunderstorm electrification, lightning, and their associated production
of nitrogen oxides (LNO$_x$) couple convective dynamics, cloud
microphysics, and atmospheric chemistry through the electric field
$\boldsymbol{E}$. Resolving $\boldsymbol{E}$ inside a
convection-resolving model requires solving an elliptic boundary-value
problem at each physics step: the Poisson equation for the electrostatic
potential $\varphi$ arising from the space-charge distribution $\rho$.
Storm-interior electric field strengths grow toward the
$\sim 150$–$300$ kV m$^{-1}$ dielectric-breakdown range at storm
altitudes, triggering breakdown and lightning channels that rapidly
redistribute charge and emit NO$_x$ {cite:p}`macgorman1998electrical`.
On larger scales, the same elliptic operator, with permittivity
$\varepsilon$ replaced by conductivity $\sigma(\boldsymbol{x})$,
describes the global atmospheric electric circuit
{cite:p}`williams2009global,rycroft2008investigation`.

The Model for Prediction Across Scales Atmosphere (MPAS-A) is a
nonhydrostatic atmospheric model discretized on centroidal Voronoi
tessellations of the sphere with C-grid staggering
{cite:p}`skamarock2012multiscale,ringler2010unified,thuburn2009numerical`.
Its horizontal discrete operators — divergence, gradient, curl — are
built from a small set of mesh metrics (cell areas, edge lengths,
cell-center separations) and have been extensively characterized for
hyperbolic dynamics. To our knowledge, no *elliptic* Poisson operator
has previously been constructed for MPAS-A's atmospheric dynamics, and
no production atmospheric model based on unstructured Voronoi
tessellations currently supports the Poisson equation as a first-class
capability.

This narrative establishes such a capability. A discrete Laplacian is
defined on MPAS-A's unstructured-horizontal, layered-vertical mesh by
applying the TRiSK div-of-grad construction in the horizontal and a
centered-difference stencil in the vertical. After multiplication by
cell volume, the resulting operator is sparse, symmetric, and negative
definite — its negation is symmetric positive definite (SPD) — and is
amenable to standard Krylov iteration. The linear system is solved by
preconditioned conjugate gradients (PCG). The method is evaluated
through five complementary studies. Four — verification, idealized
applications, mesh sensitivity, and an end-to-end demonstration — are
reported here as completed results; the fifth, parallel scaling, is
specified but retained as a follow-on measurement:

1. **Verification** ([](verification.md)) — via the method of
   manufactured solutions (MMS) in Cartesian and spherical geometries.
2. **Idealized applications** ([](idealized-applications.md)) —
   qualitative demonstration on canonical storm-like charge
   distributions.
3. **Mesh sensitivity** ([](mesh-sensitivity.md)) — characterization of
   convergence-rate sensitivity to mesh distortion and to
   variable-resolution mesh-refinement transitions.
4. **Parallel scaling** ([](parallel-performance.md)) — strong- and
   weak-scaling behavior, retained as a follow-on measurement.
5. **End-to-end example** ([](end-to-end-supercell.md)) — an
   illustration on a simulated supercell thunderstorm using a
   lightweight electrification stub.

The results reported here establish the in-tree operator and solver,
verify second-order horizontal convergence on the Cartesian MMS
sequence, characterize the spherical pentagon-defect limitation,
quantify the operator's mesh-sensitivity, and demonstrate physically
interpretable diagnostic fields for idealized and dynamically evolving
charge sources. The central numerical finding is that the operator's
second-order accuracy is a property of near-uniform hexagonal meshes: on
the irregular general centroidal-Voronoi tessellations used for
variable-resolution MPAS-A configurations — a principal motivation for
the TRiSK framework — the solved-potential convergence collapses to
$\sim\!0.35$, and restoring second order requires a higher-order or
finite-element reconstruction rather than a different finite-volume
weight. Parallel scaling and preconditioner performance remain follow-on
measurements.

This narrative proceeds as follows. [](governing-equations.md) states
the governing equations and boundary conditions. [](discretization.md)
defines the discrete operator. [](solver.md) describes the solver.
[](verification.md) through [](end-to-end-supercell.md) present the
studies, and [](terrain-following.md) extends the operator to
terrain-following meshes. [](discussion.md) discusses applications and
[](conclusions.md) concludes. Detailed step-by-step mathematical
derivations are provided in the [appendices](sec-appendices).
