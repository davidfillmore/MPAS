(app-shdom-sources)=
# Thermal and solar pseudo-source terms

The source function {eq}`eq-J-physical` contains three physical
contributions: in-scattering, thermal emission, and the
singly-scattered direct solar beam. The first was treated in
[](scattering.md) and absorbed into the multiplicative kernel
$\omega \chi_l / (2l+1)$ of Eq. (A5). The second and third
contribute the additive term $S_{lm}$ in that equation, and are
derived here, reproducing Eqs. (A6) and (A7) of Evans (1998).

(app-shdom-sources-thermal)=
## Thermal emission: Eq. (A6)

Within the Kirchhoff (local thermodynamic equilibrium)
approximation and neglecting polarization, thermal emission is
isotropic and has radiance equal to $(1 - \omega)\,B(T)$, where
$B$ is the Planck function at the local kinetic temperature $T$
and $1 - \omega$ is the single-scattering absorption coefficient
fraction. Writing the thermal piece of the source function as

$$
  J^{\mathrm{thermal}}(\hat{\Omega}) \;=\; (1 - \omega)\,B(T)
  \qquad \text{(isotropic)},
$$

its spherical-harmonic coefficients are obtained by direct
projection against {eq}`eq-I-lm-continuous`,

$$
  S_{lm}^{\mathrm{thermal}} \;=\; (1 - \omega)\,B(T)
    \int_{0}^{2\pi}\!\int_{-1}^{1} Y_{lm}(\mu, \phi)\,d\mu\,d\phi.
$$ (eq-S-thermal-integral)

Evaluating the integral of $Y_{lm}$ over the full sphere: for
$(l, m) \ne (0, 0)$, $Y_{lm}$ is orthogonal to the constant
function $1$ (which is proportional to $Y_{00}$), so the integral
vanishes. For $(l, m) = (0, 0)$,
$Y_{00}(\mu, \phi) = \Lambda_{00}\cdot u(0) = \Lambda_{00}$ is a
constant; orthonormality
$\int\!\int Y_{00}^{2}\,d\mu\,d\phi = 1$ combined with
$\int\!\int d\mu\,d\phi = 4\pi$ yields
$\Lambda_{00} = 1/\sqrt{4\pi}$ and therefore

$$
  \int_{0}^{2\pi}\!\int_{-1}^{1} Y_{00}(\mu, \phi)\,d\mu\,d\phi
  \;=\; \frac{1}{\sqrt{4\pi}}\cdot 4\pi
  \;=\; (4\pi)^{1/2}.
$$ (eq-int-Y00)

Substituting back,

$$
  \boxed{\;
    S_{lm}^{\mathrm{thermal}}
    \;=\; (1 - \omega)\,B(T)\,(4\pi)^{1/2}\,\delta_{l0}\,\delta_{m0},
  \;}
$$ (eq-A6)

which is Eq. (A6) of Evans (1998). The thermal source deposits all
its strength in the $(0, 0)$ mode — which is exactly
$(1/(4\pi)^{1/2})\times S_{00}^{\mathrm{thermal}}\,Y_{00}(\hat{\Omega})
= (1 - \omega) B(T)$, recovering the isotropic emission as
expected.

(app-shdom-sources-solar-motivation)=
## Solar direct beam and pseudo-source: motivation

Atmospheric radiation problems driven by solar insolation have a
collimated boundary condition: at the top of the atmosphere the
downwelling radiance is confined to a single direction
$(\mu_0, \phi_0)$ with $\mu_0 = \cos\theta_0$ the cosine of the
solar zenith angle, carrying an irradiance $F_0$ on a horizontal
surface. Were one to carry this direct beam inside the discrete
ordinate representation, the collimated radiance would span a
single angular grid point surrounded by cells in which $I = 0$,
leading to catastrophic angular aliasing. SHDOM sidesteps this by
the standard *direct–diffuse decomposition*: split the total
radiance into a collimated direct component (integrated
analytically along the solar ray by Beer's law) plus a diffuse
remainder (solved by SHDOM),

$$
  I(\mathbf{r}, \hat{\Omega})
  \;=\; I^{\mathrm{direct}}(\mathbf{r}, \hat{\Omega})
    \;+\; I^{\mathrm{diff}}(\mathbf{r}, \hat{\Omega}),
$$

where

$$
  I^{\mathrm{direct}}(\mathbf{r}, \hat{\Omega})
  \;=\; \frac{F_0}{\mu_0}\,e^{-\tau_s(\mathbf{r})}\,
        \delta\bigl(\hat{\Omega} - \hat{\Omega}_0\bigr),
$$ (eq-I-direct)

with $\tau_s(\mathbf{r}) \equiv \int_{\mathbf{r}}^{\mathrm{sun}}
k\,ds'$ the optical depth from $\mathbf{r}$ to space along the
solar ray. The prefactor $F_0 / \mu_0$ converts the horizontal
irradiance $F_0$ to slanted-path radiance. The diffuse component
$I^{\mathrm{diff}}$ satisfies the same radiative transfer
equation {eq}`eq-A11` but with the direct-beam scattering
contribution appearing as a *pseudo-source* in the source
function. Only the diffuse component is represented in the
spherical-harmonic / discrete-ordinate machinery.

(app-shdom-sources-solar)=
## Solar pseudo-source: Eq. (A7)

The solar pseudo-source is the scattering integral of the direct
beam {eq}`eq-I-direct`:

$$
  S^{\odot}(\hat{\Omega})
  \;=\; \frac{\omega}{4\pi}
    \int_{4\pi} P(\hat{\Omega}, \hat{\Omega}')\,
      I^{\mathrm{direct}}(\hat{\Omega}')\,d\hat{\Omega}'
  \;=\; \frac{\omega}{4\pi}\,\frac{F_0}{\mu_0}\,e^{-\tau_s}\,
        P(\hat{\Omega}, \hat{\Omega}_0),
$$ (eq-S-solar-angular)

where the $\delta$-function has collapsed the scattering integral
onto $\hat{\Omega}' = \hat{\Omega}_0$. Projecting
{eq}`eq-S-solar-angular` onto the spherical-harmonic basis and
substituting the Legendre expansion {eq}`eq-A4` and addition
theorem {eq}`eq-addition-real`,

$$
\begin{aligned}
  S_{lm}^{\odot}
  &\;=\; \frac{\omega}{4\pi}\,\frac{F_0}{\mu_0}\,e^{-\tau_s}
    \int_{4\pi} P(\hat{\Omega}, \hat{\Omega}_0)\,Y_{lm}(\hat{\Omega})\,d\hat{\Omega} \\
  &\;=\; \frac{\omega}{4\pi}\,\frac{F_0}{\mu_0}\,e^{-\tau_s}
    \sum_{l'} \chi_{l'}\,\frac{4\pi}{2l'+1}
    \sum_{m'} Y_{l'm'}(\hat{\Omega}_0)
    \int_{4\pi} Y_{l'm'}(\hat{\Omega})\,Y_{lm}(\hat{\Omega})\,d\hat{\Omega} \\
  &\;=\; \frac{\omega}{4\pi}\,\frac{F_0}{\mu_0}\,e^{-\tau_s}\,
    \chi_l\,\frac{4\pi}{2l+1}\,Y_{lm}(\hat{\Omega}_0),
\end{aligned}
$$

where orthonormality {eq}`eq-Ylm-ortho` collapses the $l', m'$
sums. Simplifying,

$$
  \boxed{\;
    S_{lm}^{\odot}
    \;=\; \frac{F_0}{\mu_0}\,e^{-\tau_s}\,
          Y_{lm}(\mu_0, \phi_0)\,
          \frac{\omega\,\chi_l}{2l+1},
  \;}
$$ (eq-A7)

which is Eq. (A7) of Evans (1998). The solar pseudo-source
distributes its strength across all $(l, m)$ modes according to
$Y_{lm}$ evaluated at the solar direction, with the same
$\omega \chi_l / (2l+1)$ multiplier that appears in the
in-scattering term of {eq}`eq-A5`. Combining {eq}`eq-A5`,
{eq}`eq-A6`, and {eq}`eq-A7` yields the full source function
update used by SHDOM in step 4 of the Picard iteration
(see {ref}`app-shdom-rte-picard`).
