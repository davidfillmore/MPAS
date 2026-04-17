(app-shdom-scat)=
# Scattering integral in spherical harmonic space

The efficiency of SHDOM rests on the observation that the
scattering integral — the most expensive angular operation in the
radiative transfer equation — reduces to a *pointwise
multiplication* in spherical-harmonic space when the phase
function depends only on the scattering angle. This section
derives that reduction step by step, reproducing Eq. (1) on
page 430 and Eqs. (A4)–(A5) of Evans (1998).

(app-shdom-scat-P-expand)=
## The scattering integral and the Legendre expansion of the phase function

The scattering contribution to the source function,

$$
  J^{\mathrm{scat}}(\hat{\Omega})
  \;=\; \frac{\omega}{4\pi}
    \int_{4\pi} P(\hat{\Omega}, \hat{\Omega}')\,
      I(\hat{\Omega}')\,d\hat{\Omega}',
$$ (eq-Jscat-def)

couples scattering to every direction $\hat{\Omega}'$ through the
phase function $P$. For unpolarized radiative transfer in a medium
whose scatterers are randomly oriented, the phase function depends
only on the scattering angle $\Theta$ between incident and
scattered directions,

$$
  \cos\Theta \;=\; \hat{\Omega}\cdot\hat{\Omega}'
           \;=\; \mu\mu' + \sqrt{1-\mu^{2}}\sqrt{1-\mu'^{2}}\,\cos(\phi - \phi'),
$$ (eq-cosTheta)

so $P$ admits an expansion in Legendre polynomials of $\cos\Theta$
alone:

$$
  \boxed{\;
    P(\cos\Theta) \;=\; \sum_{l=0}^{\infty} \chi_l\,P_l(\cos\Theta),
  \;}
$$ (eq-A4)

reproducing Eq. (A4) of Evans (1998). Here $P_l$ is the standard
(unnormalized) Legendre polynomial of degree $l$, with
$\int_{-1}^{1} P_l(x) P_{l'}(x)\,dx = \tfrac{2}{2l+1}\delta_{ll'}$,
and $\chi_l$ are the dimensionless Legendre expansion coefficients
of the phase function. Normalization of $P$ to
$\int_{4\pi} P\,d\hat{\Omega}/(4\pi) = 1$ fixes $\chi_0 = 1$,
since $P_l$ integrates to zero over $\mu \in [-1, 1]$ for
$l \ge 1$ and to 2 for $l = 0$.

(app-shdom-scat-addition)=
## Addition theorem for real orthonormal spherical harmonics

The Legendre polynomial of $\cos\Theta$ decouples into a sum of
products of spherical harmonics evaluated separately at
$\hat{\Omega}$ and $\hat{\Omega}'$. For *complex* orthonormal
spherical harmonics $Y_{lm}^{\mathbb{C}}(\mu, \phi)$ this is the
classical addition theorem

$$
  P_l(\hat{\Omega}\cdot\hat{\Omega}')
  \;=\; \frac{4\pi}{2l+1}\,
    \sum_{m=-l}^{l} Y_{lm}^{\mathbb{C}\,*}(\hat{\Omega})\,
                    Y_{lm}^{\mathbb{C}}(\hat{\Omega}').
$$ (eq-addition-complex)

The real orthonormal basis of
{ref}`app-shdom-sphharm-Ylm` is related to the
complex basis by a unitary change of basis within each fixed $l$:
for $m > 0$ one takes
$\bigl(Y_{l,m}^{\mathbb{C}} + (-1)^{m}Y_{l,-m}^{\mathbb{C}}\bigr)/\sqrt{2}$,
for $m < 0$ one takes the corresponding sine combination, and for
$m = 0$ one leaves $Y_{l0}^{\mathbb{C}}$ unchanged. Because the
sum $\sum_m Y^{\mathbb{C}\,*}(\hat{\Omega}) Y^{\mathbb{C}}(\hat{\Omega}')$
over a complete orthonormal set within a fixed $l$ block is
unitary-invariant, substituting the real combinations into
{eq}`eq-addition-complex` gives the real-basis version

$$
  \boxed{\;
    P_l(\hat{\Omega}\cdot\hat{\Omega}')
    \;=\; \frac{4\pi}{2l+1}\,
      \sum_{m=-l}^{l} Y_{lm}(\hat{\Omega})\,Y_{lm}(\hat{\Omega}'),
  \;}
$$ (eq-addition-real)

where $Y_{lm}$ are Evans's real orthonormal spherical harmonics
{eq}`eq-Ylm-def`. (The complex conjugate has disappeared because
the real basis equals its own conjugate.)

(app-shdom-scat-derivation)=
## Scattering integral in spherical harmonic space

Expanding the radiance in the real basis,
$I(\hat{\Omega}') = \sum_{l', m'} I_{l'm'}\,Y_{l'm'}(\hat{\Omega}')$,
and combining {eq}`eq-Jscat-def`–{eq}`eq-addition-real`,

$$
\begin{aligned}
  J^{\mathrm{scat}}(\hat{\Omega})
  &\;=\; \frac{\omega}{4\pi}
    \sum_{l=0}^{\infty} \chi_l
    \int_{4\pi} P_l(\hat{\Omega}\cdot\hat{\Omega}')\,I(\hat{\Omega}')\,d\hat{\Omega}'
  \\[4pt]
  &\;=\; \frac{\omega}{4\pi}
    \sum_{l=0}^{\infty} \chi_l\,\frac{4\pi}{2l+1}
    \sum_{m=-l}^{l} Y_{lm}(\hat{\Omega})
    \int_{4\pi} Y_{lm}(\hat{\Omega}')\,I(\hat{\Omega}')\,d\hat{\Omega}'
  \\[4pt]
  &\;=\; \omega\,\sum_{l, m}\frac{\chi_l}{2l+1}\,I_{lm}\,Y_{lm}(\hat{\Omega}),
\end{aligned}
$$

where the last line uses the projection formula
{eq}`eq-I-lm-continuous` to identify the remaining angular
integral as $I_{lm}$. Reading off the $(l, m)$ spherical-harmonic
coefficient,

$$
  \boxed{\;
    J^{\mathrm{scat}}_{lm}
    \;=\; \frac{\omega\,\chi_l}{2l+1}\,I_{lm},
  \;}
$$ (eq-Jscat-lm)

which is Eq. (1) on page 430 of Evans (1998). Scattering has
become pointwise in $(l, m)$-space: each spherical-harmonic
coefficient of the radiance field is multiplied by a scalar that
depends on $l$ but neither on $m$ nor on any other coefficient.
Computing the scattering integral therefore costs
$O(N_{\mathrm{lm}}) = O(N)$ operations, in contrast to the
$O(N^{2})$ cost of direct summation in the angular domain. This
is the central reason why SHDOM's per-iteration cost is
$O(N^{3/2})$ rather than $O(N^{2})$: the angular transforms of
{ref}`app-shdom-sphharm-forward`–{ref}`app-shdom-sphharm-inverse`,
not the scattering integral, dominate the operation count.

(app-shdom-scat-full)=
## Full source function in spherical harmonic space

Adding the thermal and solar pieces of the source function
{eq}`eq-J-physical` (each derived in [](sources.md)),

$$
  \boxed{\;
    J_{lm}
    \;=\; \frac{\omega\,\chi_l}{2l+1}\,I_{lm}
    \;+\; S_{lm},
  \;}
$$ (eq-A5)

which is Eq. (A5) of Evans (1998), with

$$
  S_{lm} \;\equiv\; S_{lm}^{\mathrm{thermal}} + S_{lm}^{\odot}.
$$ (eq-S-split)

The two pieces of $S_{lm}$ are derived separately in the sections
that follow.
