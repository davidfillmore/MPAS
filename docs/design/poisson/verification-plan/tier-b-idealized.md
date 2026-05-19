(tier-b)=
# Tier B — Idealized electrostatic applications

Tier B exercises the solver on two idealized charge distributions of
atmospheric-electrostatics significance. The tests use the local
idealized supercell mesh, a doubly-periodic Cartesian domain spanning
approximately $83.5 \times 83.8 \times 20$ km with nominal 1 km
horizontal spacing. These cases are not formal convergence tests; they
check whether the computed potential and electric field have the
expected physical structure for interpretable charge distributions.

## Tier B.1 — Regularized point charge

A three-dimensional isotropic Gaussian charge distribution,

$$
\rho_{\mathrm{pt}}(\boldsymbol{x}) \;=\;
  \frac{Q}{(2\pi\sigma^2)^{3/2}}
  \exp\!\left(-\frac{|\boldsymbol{x} - \boldsymbol{x}_0|^2}{2\sigma^2}\right),
\qquad
\int \rho_{\mathrm{pt}} \, dV = Q,
$$ (eq-rho-point)

with $\sigma = 2$ km and total charge $Q = 20$ C centred at the
middle of the domain, regularizes the bare Coulomb singularity enough
to be representable on the mesh. Outside the Gaussian core
($r = |\boldsymbol{x} - \boldsymbol{x}_0| \gg \sigma$), the computed
field magnitude $|\boldsymbol{E}|$ should approach the Coulomb
far-field $|\boldsymbol{E}| \simeq Q / (4\pi\varepsilon_0 r^2)$ in
the interior region where the conducting ground, zero-flux model top,
and horizontal periodic images are still weak perturbations.

The completed Phase 1 diagnostic uses median radial bins over
$6 \le r \le 25$ km. Over this interior window, the fitted log-log
slope is $-1.94$, close to the inverse-square value $-2$, and the
median ratio to the Coulomb profile is $1.04$. At larger radii the
finite periodic box and vertical boundary conditions become visible, so
those distances are excluded from the fit. The preconditioned
conjugate-gradient (PCG) residual for the solve was
$9.93 \times 10^{-11}$.

## Tier B.2 — Thundercloud tripole

A three-Gaussian vertical stack represents the classical mature-storm
tripole. The current diagnostic source is a neutral stack with modest
charge magnitudes; it tests field topology rather than breakdown-scale
field strength.

| Layer | $z$ (km) | $Q$ (C) | $(\sigma_h, \sigma_v)$ (km) |
|---|---:|---:|---:|
| Upper positive | 10 | $+20$ | $(8.0, 1.5)$ |
| Main negative | 7 | $-40$ | $(8.0, 1.5)$ |
| Lower positive | 4 | $+20$ | $(8.0, 1.5)$ |

Each Gaussian is centred at the horizontal domain centre and the
indicated altitude, with horizontal and vertical standard deviations
$\sigma_h$ and $\sigma_v$. The completed Phase 1 diagnostic reports
the vertical cross-section of $\varphi$ and $|\boldsymbol{E}|$ through
the source axis. The maximum field in the current configuration is
2.86 kV m$^{-1}$ at $z = 8.75$ km, with PCG residual
$9.83 \times 10^{-11}$. The potential has the expected negative well
around the main negative layer and positive lobes above and below it.

The tripole result should not be read as a storm-breakdown simulation:
the field strength is well below the common $100$-$300$ kV m$^{-1}$
storm-interior threshold range. Increasing the source amplitude or
narrowing the Gaussian widths is the appropriate next step if a later
showcase requires breakdown-scale fields from an idealized analytic
source.
