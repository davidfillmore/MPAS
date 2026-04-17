(tier-b)=
# Tier B — Idealized electrostatic applications

Tier B exercises the solver on two idealized charge distributions
of atmospheric-electrostatics significance, for which analytic or
near-analytic reference solutions exist and against which the
qualitative plausibility of the solver output can be assessed.
Both tests use a doubly-periodic $200 \times 200 \times 20$ km
Cartesian domain at nominal $1$ km horizontal spacing.

## Tier B.1 — Regularized point charge

A three-dimensional isotropic Gaussian charge distribution,

$$
\rho_{\mathrm{pt}}(\boldsymbol{x}) \;=\;
  \frac{Q}{(2\pi\sigma^2)^{3/2}}
  \exp\!\left(-\frac{|\boldsymbol{x} - \boldsymbol{x}_0|^2}{2\sigma^2}\right),
\qquad
\int \rho_{\mathrm{pt}} \, dV = Q,
$$ (eq-rho-point)

with $\sigma = 2$ km and total charge $Q = 1$ C centred at the
middle of the domain, regularizes the bare Coulomb singularity
enough to be representable on the mesh. Outside the Gaussian core
(i.e. at $r = |\boldsymbol{x} - \boldsymbol{x}_0| \gg \sigma$) the
computed field magnitude $|\boldsymbol{E}|$ is expected to
converge to the Coulomb far-field
$|\boldsymbol{E}| \simeq Q / (4\pi\varepsilon_0 r^2)$. The Tier B.1
measurement is a radial cross-section of $|\boldsymbol{E}|$ from
$\boldsymbol{x}_0$ outward, compared against the analytic
far-field profile on a log–log axis; agreement to within plotting
accuracy over $r \gtrsim 3\sigma$ is the qualitative pass
criterion.

## Tier B.2 — Thundercloud tripole

A three-Gaussian vertical stack represents the classical
mature-storm tripole. The parameters reproduce observational
charge magnitudes and vertical offsets typical of deep continental
convection:

| Layer | $z$ (km) | $Q$ (C) | $(\sigma_h, \sigma_v)$ (km) |
|---|---|---|---|
| Upper positive            | 10 | $+40$ | $(5.0, 1.5)$ |
| Main negative             |  6 | $-40$ | $(5.0, 1.5)$ |
| Lower positive (screening)|  3 | $+5$  | $(5.0, 1.5)$ |

Each Gaussian is centred at horizontal domain centre and the
indicated altitude, with horizontal and vertical standard
deviations $\sigma_h$, $\sigma_v$. The Tier B.2 measurement
reports the peak $|\boldsymbol{E}|$ magnitude and its vertical
location, a vertical cross-section of $\varphi$ and
$|\boldsymbol{E}|$ through the storm axis, and the 3D isosurface
of $|\boldsymbol{E}|$ at a fixed threshold. The qualitative pass
criterion is that the peak $|\boldsymbol{E}|$ occurs near the main
negative layer and attains a magnitude in the
$100$–$300$ kV m$^{-1}$ range consistent with storm-interior
observations.
