(sec-equations)=
# Governing equations

## Continuous problem

In the quasi-static limit of Maxwell's equations, valid when the
timescales and length scales of interest satisfy
$\omega L / c \ll 1$, the electric field is irrotational,
$\nabla \times \boldsymbol{E} = 0$, and therefore derives from a
scalar potential,

$$
\boldsymbol{E}(\boldsymbol{x}, t)
  = -\nabla \varphi(\boldsymbol{x}, t).
$$ (eq-E-from-phi)

Gauss's law,
$\nabla \cdot (\varepsilon \boldsymbol{E}) = \rho$, then yields the
Poisson equation for $\varphi$ in divergence form,

$$
\nabla \cdot \left[ \varepsilon(\boldsymbol{x}) \,
                    \nabla \varphi(\boldsymbol{x}, t) \right]
  = -\rho(\boldsymbol{x}, t),
$$ (eq-poisson)

where $\varepsilon$ is the permittivity of the medium and $\rho$
is the free space-charge density. In the lower atmosphere, the
relative permittivity of air departs from unity by less than
$10^{-3}$, so we take $\varepsilon = \varepsilon_0 = 8.854 \times
10^{-12}$ F m$^{-1}$ constant throughout this work. The
divergence form is retained nonetheless, so that the same code
supports future work in which $\varepsilon$ is replaced by a
spatially variable conductivity $\sigma(\boldsymbol{x})$ to solve
the current-continuity equation $\nabla \cdot (\sigma \nabla
\varphi) = S$ for the global atmospheric electric circuit.

Throughout, $\varphi$ is in volts (V), $\rho$ is in C m$^{-3}$,
and $\boldsymbol{E}$ is in V m$^{-1}$.

## Boundary conditions

Let $\Omega$ denote the 3D model domain, with horizontal extent
given by the MPAS-A mesh and vertical extent spanning $z \in
[z_{\mathrm{g}}(\boldsymbol{x}_h), z_{\mathrm{top}}]$ where
$z_{\mathrm{g}}$ is ground altitude as a function of horizontal
position $\boldsymbol{x}_h$ and $z_{\mathrm{top}}$ is the model
top. We impose

$$
\varphi = 0
  \quad \text{on } z = z_{\mathrm{g}}
  \quad \text{(perfectly conducting Earth, Dirichlet)},
$$ (eq-bc-ground)

$$
\partial_z \varphi = 0
  \quad \text{on } z = z_{\mathrm{top}}
  \quad \text{(Neumann)},
$$ (eq-bc-top)

plus periodic conditions in the horizontal on Cartesian
doubly-periodic meshes or no lateral boundary conditions on global
spherical meshes (the sphere is closed). Regional limited-area
configurations with lateral Dirichlet, Neumann, or Robin
conditions are deferred to future work.

The Dirichlet condition {eq}`eq-bc-ground` eliminates the constant
null mode of $\nabla \cdot (\varepsilon \nabla)$, ensuring a
unique solution $\varphi$ for any admissible $\rho$.

## Source term and coupling

The space-charge density $\rho$ is supplied externally to the
Poisson solver. In the idealized-physics configurations used for
the benchmarks of [](verification.md) and
[](idealized-applications.md), $\rho$ is an analytic
function of position. In the end-to-end supercell configuration
([](end-to-end-supercell.md)), $\rho$ is
produced by a lightweight electrification stub that parameterizes
charge separation as proportional to vertical velocity and to
hydrometeor mass mixing ratios. When graupel and cloud ice are
active, the positive and negative proxies are graupel and ice; in
the Kessler supercell used here, the fallback proxies are rain and
cloud water. This is an order-of-magnitude placeholder for future
full electrification schemes of Saunders–Peck or Takahashi type
{cite:p}`saunders1991effect,takahashi1978riming`, not a physical
charging model.

The coupling of the solver into the MPAS-A time loop is one-way in
this work: $\varphi$ and $\boldsymbol{E}$ are diagnostic outputs
and do not feed back into dynamics or microphysics. Follow-on
applications of the solver — lightning parameterizations, LNO$_x$
emission, and two-way feedback experiments — are discussed in
[](discussion.md).
