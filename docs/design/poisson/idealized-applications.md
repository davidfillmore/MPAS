(sec-idealized)=
# Idealized electrostatic applications

The idealized applications exercise the solver on two canonical charge
distributions of atmospheric-electrostatics significance. The tests use
the idealized-supercell mesh, a doubly periodic Cartesian domain
spanning approximately $83.5 \times 83.8 \times 20$ km with nominal
$1$ km horizontal spacing. These cases are not formal convergence
tests; they check whether the computed potential and electric field
have the expected physical structure for interpretable charge
distributions.

## Regularized point charge

A three-dimensional isotropic Gaussian charge distribution,

$$
\rho_{\mathrm{pt}}(\boldsymbol{x}) \;=\;
  \frac{Q}{(2\pi\sigma^2)^{3/2}}
  \exp\!\left(-\frac{|\boldsymbol{x} - \boldsymbol{x}_0|^2}{2\sigma^2}\right),
\qquad
\int \rho_{\mathrm{pt}} \, dV = Q,
$$ (eq-rho-point)

with $\sigma = 2$ km and total charge $Q = 20$ C centered at the middle
of the domain, regularizes the bare Coulomb singularity enough to be
representable on the mesh. Outside the Gaussian core (i.e. at
$r = |\boldsymbol{x} - \boldsymbol{x}_0| \gg \sigma$) the computed field
magnitude $|\boldsymbol{E}|$ is expected to converge to the Coulomb
far-field $|\boldsymbol{E}| \simeq Q / (4\pi\varepsilon_0 r^2)$.
{numref}`fig-point-charge` shows per-bin medians of $|\boldsymbol{E}|$
in radial bins over $6 \le r \le 25$ km. Over this interior window the
fitted log-log slope is $-1.94$, close to the inverse-square value
$-2$, and the median ratio to the Coulomb profile is $1.04$. At larger
radii the finite periodic box and vertical boundary conditions become
visible, so those distances are excluded from the fit. The PCG residual
for the solve was $9.93 \times 10^{-11}$.

```{figure} figures/tier_B1_point_charge.png
:name: fig-point-charge
:alt: Regularized point-charge diagnostic.
:width: 78%

Regularized point-charge diagnostic on the idealized supercell mesh.
Blue markers show per-bin medians of the MPAS field magnitude, while
the solid and dashed curves show the analytic Gaussian field and
Coulomb far-field reference, respectively.
```

## Thundercloud tripole

A three-Gaussian vertical stack represents the classical mature-storm
tripole. The parameters are given below; the amplitudes are
intentionally modest in the present diagnostic source, so this case
tests field topology rather than breakdown-scale field strength.

| Layer | $z$ (km) | $Q$ (C) | $(\sigma_h, \sigma_v)$ (km) |
|---|---:|---:|---:|
| Upper positive | 10 | $+20$ | $(8.0,\,1.5)$ |
| Main negative | 7 | $-40$ | $(8.0,\,1.5)$ |
| Lower positive | 4 | $+20$ | $(8.0,\,1.5)$ |

Each Gaussian is centered at the horizontal domain center and the
indicated altitude, with horizontal and vertical standard deviations
$\sigma_h$, $\sigma_v$. The tripole diagnostic reports the vertical
cross-section of $\varphi$ and $|\boldsymbol{E}|$ through the source
axis ({numref}`fig-tripole`). The maximum field in the current
configuration is $2.86$ kV m$^{-1}$ at $z = 8.75$ km, with PCG residual
$9.83 \times 10^{-11}$. The potential has the expected negative well
around the main negative layer and positive lobes above and below it;
increasing the source amplitude or narrowing the Gaussian widths is the
appropriate next step if a later showcase requires storm-interior
fields of order $100$–$300$ kV m$^{-1}$.

```{figure} figures/tier_B2_tripole.png
:name: fig-tripole
:alt: Synthetic tripole diagnostic.
:width: 100%

Synthetic tripole diagnostic. The left panel shows the electrostatic
potential with electric-field streamlines on the vertical
source-axis section; the right panel shows $|\boldsymbol{E}|$ on the
same section. Charge-density contours indicate the positive, negative,
and positive Gaussian stack.
```
