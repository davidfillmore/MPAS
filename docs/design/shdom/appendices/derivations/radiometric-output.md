(app-shdom-radout)=
# Radiometric output quantities

Once the Picard iteration has converged, SHDOM computes
radiometric output quantities — hemispheric fluxes, mean
radiance, net flux vector, net flux convergence (heating rate),
and radiance at user-specified viewing directions — from the
converged radiance and source function fields. This section
derives Eqs. (A21)–(A27) of Evans (1998). Throughout, $\bar{I}$
denotes the diffuse mean radiance (the radiance integrated over
all directions and divided by $4\pi$) and subscripts $(lm)$
refer to the spherical-harmonic coefficients of the diffuse
field; direct-beam contributions are added back explicitly below.

(app-shdom-radout-hemflux)=
## Hemispheric fluxes from the discrete ordinates: Eq. (A21)

The upward ($+$) and downward ($-$) hemispheric flux on a
horizontal surface is, in continuous form,

$$
  F^{\pm}(\mathbf{x}) \;=\; \int_{\pm\mu > 0}
                            I(\mathbf{x}, \hat{\Omega})\,|\mu|\,d\hat{\Omega},
$$ (eq-Fpm-continuous)

where $F^{+}$ integrates over the upward hemisphere ($\mu > 0$)
and $F^{-}$ over the downward hemisphere ($\mu < 0$).
Approximating {eq}`eq-Fpm-continuous` on the reduced Gaussian
discrete-ordinate grid of
{ref}`app-shdom-sphharm-grid`,

$$
  \boxed{\;
    F^{\pm}(\mathbf{x}) \;=\;
    \sum_{j \in \mathcal{H}_{\pm}}\! w_j
    \sum_{k=1}^{N_{\phi, j}} \hat{w}_{jk}\,|\mu_j|\,I_{jk}(\mathbf{x}),
  \;}
$$ (eq-A21)

which is Eq. (A21) of Evans (1998). Here
$\mathcal{H}_{+} = \{j : \mu_j > 0\}$ selects the $N_\mu / 2$
zenith abscissae in the upper hemisphere (and $\mathcal{H}_{-}$
those in the lower), $w_j$ and $\hat{w}_{jk}$ are the zenith
Gauss–Legendre and azimuthal uniform-grid weights defined in
{ref}`app-shdom-sphharm-grid`, and
$I_{jk}(\mathbf{x})$ is the discrete-ordinate radiance at
spatial grid point $\mathbf{x}$. Because Gauss–Legendre
abscissae are symmetric about $\mu = 0$ and $|\mu|$ is an even
polynomial in $\mu$, the upward and downward integrations use
identical weight sets differing only by which abscissae are
included.

(app-shdom-radout-Imean)=
## Mean radiance from the $(0, 0)$ spherical-harmonic coefficient: Eq. (A22)

The diffuse mean radiance is

$$
  \bar{I}(\mathbf{x}) \;\equiv\; \frac{1}{4\pi}\int_{4\pi}
         I(\mathbf{x}, \hat{\Omega})\,d\hat{\Omega}.
$$

Expanding $I$ in the orthonormal basis and using
$\int Y_{lm}\,d\hat{\Omega} = (4\pi)^{1/2}\,\delta_{l0}\delta_{m0}$
(derived in
{ref}`app-shdom-sources-thermal`,
Eq. {eq}`eq-int-Y00`),

$$
  \bar{I}(\mathbf{x}) \;=\; \frac{1}{4\pi}
    \sum_{l, m} I_{lm}(\mathbf{x})\,\int Y_{lm}\,d\hat{\Omega}
  \;=\; \frac{I_{00}(\mathbf{x})}{(4\pi)^{1/2}}.
$$

Rearranging to Evans's form,

$$
  \boxed{\;
    \bar{I}(\mathbf{x})
    \;=\; \frac{1}{(4\pi)^{1/2}}\,I_{l=0, m=0}(\mathbf{x}),
  \;}
$$ (eq-A22)

which is Eq. (A22) of Evans (1998). The mean radiance requires
only the $(0, 0)$ spherical-harmonic coefficient, which is always
stored and is the cheapest of the radiometric outputs.

(app-shdom-radout-Fnet)=
## Net flux vector from the $(l = 1)$ coefficients: Eq. (A23)

The diffuse net flux vector
$\mathbf{F} = \int I\,\hat{\Omega}\,d\hat{\Omega}$ has Cartesian
components

$$
  F_x \;=\; \int I\,\sqrt{1-\mu^{2}}\,\cos\phi\,d\hat{\Omega},
  \quad
  F_y \;=\; \int I\,\sqrt{1-\mu^{2}}\,\sin\phi\,d\hat{\Omega},
  \quad
  F_z \;=\; \int I\,\mu\,d\hat{\Omega}.
$$ (eq-F-components-def)

The three direction cosines
$(\sqrt{1-\mu^{2}}\cos\phi, \sqrt{1-\mu^{2}}\sin\phi, \mu)$ are
linear combinations of the three $l = 1$ spherical harmonics
$(Y_{1, 1}, Y_{1, -1}, Y_{1, 0})$. Specifically, projecting each
direction cosine onto $Y_{lm}$:

$$
\begin{aligned}
  \int Y_{lm}(\hat{\Omega})\,\mu\,d\hat{\Omega}
    &\;=\; (4\pi/3)^{1/2}\,\delta_{l1}\,\delta_{m0}, \\
  \int Y_{lm}(\hat{\Omega})\,\sqrt{1-\mu^{2}}\,\cos\phi\,d\hat{\Omega}
    &\;=\; \pm(4\pi/3)^{1/2}\,\delta_{l1}\,\delta_{m, +1}, \\
  \int Y_{lm}(\hat{\Omega})\,\sqrt{1-\mu^{2}}\,\sin\phi\,d\hat{\Omega}
    &\;=\; \pm(4\pi/3)^{1/2}\,\delta_{l1}\,\delta_{m, -1},
\end{aligned}
$$

where the signs on the last two lines depend on whether
$Y_{1, \pm 1}$ is defined with or without the Condon–Shortley
phase factor $(-1)^{m}$. (The first line is sign-unambiguous:
$Y_{1, 0} \propto P_{1}(\mu) = \mu$ has no Condon–Shortley
factor.) The magnitude $(4\pi/3)^{1/2}$ is derived by direct
computation: for the $F_z$ case,

$$
  \int Y_{1, 0}\,\mu\,d\hat{\Omega}
  \;=\; \sqrt{\tfrac{3}{4\pi}}\!\int_{-1}^{1}\!\!\mu^{2}d\mu\!\int_{0}^{2\pi}\!\!d\phi
  \;=\; \sqrt{\tfrac{3}{4\pi}}\cdot\tfrac{2}{3}\cdot 2\pi
  \;=\; \bigl(\tfrac{4\pi}{3}\bigr)^{1/2},
$$

using the orthonormal convention
$Y_{1, 0} = \sqrt{3/(4\pi)}\,\mu$. Substituting into
{eq}`eq-F-components-def` via the expansion
$I = \sum I_{lm} Y_{lm}$,

$$
  \boxed{\;
    F_x \;=\; \mp\bigl(\tfrac{4\pi}{3}\bigr)^{1/2} I_{1, 1},
    \qquad
    F_y \;=\; \mp\bigl(\tfrac{4\pi}{3}\bigr)^{1/2} I_{1, -1},
    \qquad
    F_z \;=\; \bigl(\tfrac{4\pi}{3}\bigr)^{1/2} I_{1, 0},
  \;}
$$ (eq-A23)

reproducing Eq. (A23) of Evans (1998) up to the $\pm$ choice on
$F_x, F_y$, which is fixed by the sign convention on the
$m = \pm 1$ real spherical harmonics. Evans's published formulas
use the negative signs for $F_x$ and $F_y$ (i.e. Condon–Shortley
phase in the transverse components); any particular SHDOM
implementation must pick one convention and stick with it. The
$z$-component sign is unambiguous:
$F_z = +(4\pi/3)^{1/2} I_{1, 0}$.

(app-shdom-radout-direct)=
## Direct-beam additions: Eqs. (A24)–(A25)

The radiometric quantities above were derived for the diffuse
radiance. Under the direct–diffuse decomposition
(see {ref}`app-shdom-sources-solar-motivation`), the
direct solar beam contributes separately and must be added back.
Using the collimated direct radiance
$I^{\mathrm{direct}} = (F_0/\mu_0)\,e^{-\tau_s}\,
\delta(\hat{\Omega} - \hat{\Omega}_0)$:

**Mean-radiance contribution.**
$\bar{I}^{\odot} = (1/4\pi)\int I^{\mathrm{direct}}\,d\hat{\Omega}
= (1/4\pi)\,(F_0/\mu_0)\,e^{-\tau_s}$,

$$
  \boxed{\;
    \bar{I}^{\odot}(\mathbf{x})
    \;=\; \frac{F_0\,e^{-\tau_s}}{4\pi\,\mu_0},
  \;}
$$ (eq-A24)

reproducing Eq. (A24) of Evans (1998).

**Flux-vector contribution.**
$\mathbf{F}^{\odot} = \int I^{\mathrm{direct}}\,\hat{\Omega}\,d\hat{\Omega}
= (F_0/\mu_0)\,e^{-\tau_s}\,\hat{\Omega}_0$, with components

$$
  \boxed{\;
    F^{\odot}_x \;=\; \frac{F_0}{\mu_0}\,e^{-\tau_s}\,\sin\theta_0\cos\phi_0,
    \quad
    F^{\odot}_y \;=\; \frac{F_0}{\mu_0}\,e^{-\tau_s}\,\sin\theta_0\sin\phi_0,
    \quad
    F^{\odot}_z \;=\; -F_0\,e^{-\tau_s},
  \;}
$$ (eq-A25)

reproducing Eq. (A25) of Evans (1998). The $F^{\odot}_z$
component reduces to $-F_0\,e^{-\tau_s}$ because the
$\hat{\Omega}_0$ direction has
$\hat{\Omega}_0\cdot\hat{z} = -\mu_0$ (the direct beam propagates
downward, opposite to the $\mu > 0$ upward convention), and
$(F_0/\mu_0)\,e^{-\tau_s}\,(-\mu_0) = -F_0\,e^{-\tau_s}$.

(app-shdom-radout-convergence)=
## Net flux convergence (heating rate): Eq. (A26)

The *net flux convergence*
$-\nabla\cdot\mathbf{F}_{\mathrm{net}}$ is the rate at which
radiative energy is deposited per unit volume and is proportional
(via the heat capacity) to the radiative heating rate of the
atmosphere.

**Diffuse contribution.** Start from the angular integral of the
RTE {eq}`eq-A11`:

$$
  \nabla\cdot\mathbf{F}^{\mathrm{diff}}
  \;=\; \int_{4\pi}\hat{\Omega}\cdot\nabla I\,d\hat{\Omega}
  \;=\; \int_{4\pi}\frac{dI}{ds}\Big|_{\hat{\Omega}}\,d\hat{\Omega}
  \;=\; -k\!\int_{4\pi}\bigl[I - J\bigr]\,d\hat{\Omega}.
$$

Using the full source function expression
{eq}`eq-J-physical` (suppressing the direct-beam direct-radiance
piece, which is already outside the diffuse solver),

$$
\begin{aligned}
  \int J\,d\hat{\Omega}
  &\;=\; \frac{\omega}{4\pi}\!\int\!\!\int P(\hat{\Omega}, \hat{\Omega}')\,
         I(\hat{\Omega}')\,d\hat{\Omega}'\,d\hat{\Omega}
       + (1 - \omega) B(T)\!\int d\hat{\Omega}
       + \int S^{\odot}\,d\hat{\Omega} \\
  &\;=\; \omega\!\int I\,d\hat{\Omega}'
       + 4\pi (1 - \omega) B
       + \omega\,\frac{F_0}{\mu_0}\,e^{-\tau_s},
\end{aligned}
$$

where $\int P(\hat{\Omega}, \hat{\Omega}')\,d\hat{\Omega} = 4\pi$
(phase function normalization) and
$\int S^{\odot}\,d\hat{\Omega} = \omega\,(F_0/\mu_0)\,e^{-\tau_s}$
(computed directly from {eq}`eq-S-solar-angular` using the same
normalization identity). Combining and recognizing
$\int I\,d\hat{\Omega} = 4\pi\bar{I}$:

$$
  \nabla\cdot\mathbf{F}^{\mathrm{diff}}
  \;=\; -\,4\pi k(1 - \omega)\,\bar{I}
  \;+\; 4\pi k(1 - \omega)\,B(T)
  \;+\; k\omega\,\frac{F_0}{\mu_0}\,e^{-\tau_s}.
$$ (eq-divF-diff)

**Direct-beam contribution.** The direct-beam flux vector is
$\mathbf{F}^{\odot} = (F_0/\mu_0)\,e^{-\tau_s}\,\hat{\Omega}_0$
with $\hat{\Omega}_0$ constant in space. Taking the divergence,

$$
  \nabla\cdot\mathbf{F}^{\odot}
  \;=\; \frac{F_0}{\mu_0}\,\hat{\Omega}_0\cdot\nabla e^{-\tau_s}
  \;=\; -\,\frac{F_0}{\mu_0}\,e^{-\tau_s}\,\bigl(\hat{\Omega}_0\cdot\nabla\tau_s\bigr)
  \;=\; -\,k\,\frac{F_0}{\mu_0}\,e^{-\tau_s},
$$ (eq-divF-direct)

using $\hat{\Omega}_0\cdot\nabla\tau_s = +k$ (moving along the
solar ray *away* from the sun increases $\tau_s$ at rate $k$).

**Combined.** Adding {eq}`eq-divF-diff` and {eq}`eq-divF-direct`,

$$
\begin{aligned}
  \nabla\cdot\mathbf{F}_{\mathrm{net}}
  &\;=\; -\,4\pi k(1 - \omega)\,\bar{I}
     \;+\; 4\pi k(1 - \omega)\,B(T)
     \;+\; k\omega\,\frac{F_0}{\mu_0}\,e^{-\tau_s}
     \;-\; k\,\frac{F_0}{\mu_0}\,e^{-\tau_s} \\[4pt]
  &\;=\; -\,k(1 - \omega)\,\Bigl[\,4\pi\bar{I} + \frac{F_0}{\mu_0}\,e^{-\tau_s} - 4\pi B\,\Bigr].
\end{aligned}
$$

Negating,

$$
  \boxed{\;
    -\,\nabla\cdot\mathbf{F}_{\mathrm{net}}
    \;=\; k\,(1 - \omega)\,\Bigl[\,4\pi\,\bar{I}
          \;+\; \frac{F_0}{\mu_0}\,e^{-\tau_s}
          \;-\; 4\pi\,B(T)\,\Bigr],
  \;}
$$ (eq-A26)

reproducing Eq. (A26) of Evans (1998). The physical
interpretation is a local energy balance: the bracketed quantity
is the *total radiation field integrated over direction* (the
diffuse $4\pi\bar{I}$ plus the direct-beam solid-angle-integrated
intensity $(F_0/\mu_0)e^{-\tau_s}$) minus the *locally emitted*
radiation $4\pi B$, and the prefactor $k(1 - \omega)$ is the
absorption coefficient. Net energy deposited per unit volume is
thus absorption $-$ emission.

(app-shdom-radout-TMS)=
## Radiance at a specified viewing direction (TMS): Eq. (A27)

Radiance at a specified viewing direction $(\mu, \phi)$ is
computed by integrating the source function along the outgoing
ray from the point of observation backward into the medium. For
solar problems with delta-M scaling enabled, the source function
must be corrected by the TMS (Nakajima–Tanaka) procedure described
in {ref}`app-shdom-deltaM-tms`: the $\delta$-M-scaled
solver supplies the multiply-scattered source accurately, but the
single-scattered solar piece must be replaced by its un-truncated
form evaluated at the exact scattering angle
$\cos\Theta = \hat{\Omega}\cdot\hat{\Omega}_0$.

The composite formula is {eq}`eq-A27` derived in
{ref}`app-shdom-deltaM-tms`; we reproduce it here for
completeness:

$$
\begin{aligned}
  J(\mu, \phi) \;=\;\,
  &\sum_{l, m}\,J'_{lm}\,Y_{lm}(\mu, \phi) \\
  &\;-\; \sum_{l, m}\,
    \frac{F_0}{\mu_0}\,e^{-\tau_s}\,
    \frac{\omega'\,\chi'_l}{2l+1}\,
    Y_{lm}(\mu_0, \phi_0)\,Y_{lm}(\mu, \phi) \\
  &\;+\; \frac{F_0}{\mu_0}\,\frac{e^{-\tau_s}}{1 - f\omega}\,\omega\,
    \sum_l \chi_l\,P_l(\cos\Theta),
\end{aligned}
$$

where $J'_{lm}$ are the $\delta$-M-scaled converged source-
function coefficients and the three terms have the interpretation
(a) multi-scatter + $\delta$-M-single-scatter together,
(b) subtract off the $\delta$-M single-scatter, (c) add back the
un-truncated single-scatter evaluated at the exact scattering
angle. The radiance at the specified viewing angle is then
obtained by integrating $J(\mu, \phi)$ along the outgoing ray
through the medium via the cell-integration formulas of
[](cell-integration.md), with the multi-cell path split into
segments of optical depth at most $0.1$ to control integration
error (Evans, 1998, appendix §f).

---

This concludes the SHDOM-method derivation appendices.
Radiometric outputs delivered by SHDOM on each gridpoint of the
base grid comprise $\bar{I}$, the net flux vector, the
hemispheric fluxes, the net flux convergence, and — for
user-specified viewing directions and locations — the
TMS-corrected radiance at that direction. All quantities are
assembled from the converged spherical-harmonic source function
$J_{lm}$ and the associated discrete-ordinate radiance $I_{jk}$
via the formulas of this section.
