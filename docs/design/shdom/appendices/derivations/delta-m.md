(app-shdom-deltaM)=
# Delta-M scaling and the TMS correction

Cloud and aerosol phase functions are typically sharply
forward-peaked: for water clouds the asymmetry parameter $g$ is
$0.85$–$0.87$ in the solar, with $\chi_l$ decaying slowly as $l$
increases. At the truncation level $L = N_\mu - 1 \sim 16$ typical
of SHDOM production runs, the Legendre expansion {eq}`eq-A4`
simply cannot resolve the forward peak, and the single-scattered
radiance computed from the truncated
$P(\cos\Theta) = \sum_{l=0}^{L} \chi_l P_l$ suffers large errors
in the forward hemisphere. Two standard remedies are used together
in SHDOM:

- *Delta-M scaling* {cite:p}`wiscombe1977deltam`: subtract an
  idealized forward-delta peak of strength $f$ from $P$ and
  rescale the optical properties accordingly; the smoothed
  remainder $P'$ has a Legendre series that truncates cleanly at
  $L$.
- *TMS correction* {cite:p}`nakajima1988tms`: for the
  singly-scattered *solar* contribution to the output radiance at
  a specified viewing direction, replace the $\delta$-M-scaled
  phase function by the un-truncated original phase function
  evaluated at the exact scattering angle, keeping $\delta$-M
  scaling only for the multiply-scattered piece.

We derive the $\delta$-M scaling in full and state the TMS
composite formula Eq. (A27) with an interpretation of each term;
the full derivation of (A27) follows {cite:t}`nakajima1988tms` and
is not reproduced here.

(app-shdom-deltaM-scaling)=
## Delta-M scaling: decomposition of the phase function

Write the phase function as a weighted sum of a forward-delta peak
and a smoother residual,

$$
  P(\cos\Theta)
  \;=\; 2f\,\delta(1 - \cos\Theta)
  \;+\; (1 - f)\,P'(\cos\Theta),
$$ (eq-P-deltaM-decomp)

where $f \in [0, 1)$ is the fraction of the phase function
contained in the forward peak, the factor $2$ is chosen so that
the decomposition preserves the normalization
$\int P\,d\mu/2 = \int P'\,d\mu/2 = 1$ (since
$\int_{-1}^{1}\delta(1-\mu)\,d\mu = 1$, the delta piece
contributes $2f \cdot 1 = 2f$ which halves to $f$ when divided by
2), and $P'$ is the normalized residual phase function.
Legendre-expanding each side of {eq}`eq-P-deltaM-decomp` using the
identity

$$
  \delta(1 - \mu) \;=\; \sum_{l=0}^{\infty} \frac{2l+1}{2}\,P_l(\mu),
$$ (eq-delta-Legendre)

(which follows from the standard Legendre-polynomial completeness
relation $\sum_l (2l+1) P_l(\mu) P_l(\mu')/2 = \delta(\mu - \mu')$
at $\mu' = 1$ using $P_l(1) = 1$), and matching $P_l$
coefficients,

$$
  \chi_l \;=\; f\,(2l+1) \;+\; (1 - f)\,\chi'_l
  \qquad (\text{all } l \ge 0),
$$ (eq-chi-match)

which inverts to

$$
  \chi'_l \;=\; \frac{\chi_l - f\,(2l+1)}{1 - f}.
$$ (eq-chi-prime)

Truncation at $L$ requires $\chi'_l = 0$ for all $l > L$; the
standard choice {cite:p}`wiscombe1977deltam` enforces this at
$l = L + 1$ by selecting

$$
  f \;=\; \frac{\chi_{L+1}}{2L + 3}.
$$ (eq-f-choice)

The scaled Legendre coefficients retained in the SHDOM solver are
then

$$
  \chi'_l
  \;=\; \frac{\chi_l - \chi_{L+1}\,\dfrac{2l+1}{2L+3}}
             {1 - \dfrac{\chi_{L+1}}{2L+3}},
  \qquad l = 0, 1, \ldots, L.
$$ (eq-chi-prime-final)

By construction $\chi'_0 = 1$ (the $l = 0$ normalization is
preserved because $\chi_0 = 1$ and $f(2\cdot 0 + 1) = f$, so
$\chi'_0 = (1 - f)/(1 - f) = 1$).

(app-shdom-deltaM-optprops)=
## Scaled single-scattering albedo and optical depth

The delta-function piece of the decomposition represents
forward-scattering events that leave photon direction unchanged
and therefore are physically indistinguishable from no scattering
at all. The observable radiative properties — extinction,
absorption, effective scattering — must be rescaled so that the
$\delta$-function events are removed from the scattering
bookkeeping:

$$
\begin{aligned}
  \text{absorption coefficient:} \quad
  \beta_{\mathrm{abs}} &\;=\; \beta\,(1 - \omega)
  \qquad (\text{unchanged}), \\[4pt]
  \text{effective scattering:} \quad
  \beta_{\mathrm{scat}}' &\;=\; \beta\,\omega\,(1 - f)
  \qquad (\text{only the }(1 - f)\text{ portion is real scattering}), \\[4pt]
  \text{scaled extinction:} \quad
  \beta' &\;\equiv\; \beta_{\mathrm{abs}} + \beta_{\mathrm{scat}}'
          \;=\; \beta\,(1 - \omega f).
\end{aligned}
$$

Scaled single-scattering albedo and optical depth follow
immediately,

$$
  \boxed{\;
    \omega' \;=\; \frac{\omega\,(1 - f)}{1 - \omega f},
    \qquad
    \tau' \;=\; \tau\,(1 - \omega f).
  \;}
$$ (eq-omega-tau-prime)

Equations {eq}`eq-chi-prime-final` and {eq}`eq-omega-tau-prime`
together are Wiscombe's delta-M scaling. The SHDOM solver carries
$\chi'_l$, $\omega'$, $\tau'$ internally in its optical-property
fields; the un-scaled $\chi_l$, $\omega$, $\tau$ are retained only
to construct the TMS correction to the output radiance.

(app-shdom-deltaM-tms-motivation)=
## TMS correction: motivation

Flux and hemispheric-integrated quantities are insensitive to the
angular detail of the phase function and are accurately recovered
by the $\delta$-M-scaled SHDOM solution alone: the scaled phase
$P'$ has the same *low-order* Legendre moments as $P$ for
$l = 0, \ldots, L$ (by construction), and these low-order moments
control integral quantities.

The situation is different for *radiance* at a specified viewing
direction, particularly in or near the solar aureole. The
single-scattered solar contribution depends on the full angular
shape of $P(\cos\Theta)$ evaluated at the actual scattering angle
$\Theta$ between $\hat{\Omega}_0$ and $\hat{\Omega}$, not just on
its low-order Legendre moments. The $\delta$-M-scaled solution
replaces $P$ by $P'$ everywhere, losing this angular detail.

The TMS method of {cite:t}`nakajima1988tms` resolves this by
splitting the output radiance into multiply-scattered and
singly-scattered pieces and applying $\delta$-M scaling only to
the former. The singly-scattered solar piece is evaluated using
the un-truncated phase function $P(\cos\Theta)$ at the exact
scattering angle.

(app-shdom-deltaM-tms)=
## TMS composite formula: Eq. (A27)

Let $J'_{lm}$ denote the spherical-harmonic source-function
coefficients produced by the $\delta$-M-scaled SHDOM solver.
According to the TMS prescription, the source function used for
radiance output at a specified viewing direction $(\mu, \phi)$ is

$$
\begin{aligned}
  J(\mu, \phi) \;=\;
  &\underbrace{\sum_{l, m}\,J'_{lm}\,Y_{lm}(\mu, \phi)}_{\text{(a) full $\delta$-M solved source}}
  \;-\; \underbrace{\sum_{l, m}\,
        \frac{F_0}{\mu_0}\,e^{-\tau_s}\,
        \frac{\omega'\,\chi'_l}{2l+1}\,Y_{lm}(\mu_0, \phi_0)\,Y_{lm}(\mu, \phi)}_{\text{(b) subtract $\delta$-M single-scatter solar}} \\[4pt]
  &\;+\; \underbrace{\frac{F_0}{\mu_0}\,\frac{e^{-\tau_s}}{1 - f\omega}\,\omega\,
        \sum_{l=0}^{\infty}\chi_l\,P_l(\cos\Theta)}_{\text{(c) add un-truncated single-scatter solar}},
\end{aligned}
$$ (eq-A27)

where $\cos\Theta = \hat{\Omega}\cdot\hat{\Omega}_0$ is the exact
scattering angle between the solar direction $(\mu_0, \phi_0)$ and
the viewing direction $(\mu, \phi)$, $\tau_s$ is the un-scaled
solar optical path to the point of evaluation, and the sum in (c)
ranges over *all* coefficients of the input phase function (not
just $l \le L$). This is Eq. (A27) of Evans (1998).

Term-by-term interpretation:

- (a) The full $\delta$-M-scaled source function from the SHDOM
  solver, evaluated at the viewing direction. This term includes
  all multiple-scattering and the *$\delta$-M-scaled*
  single-scatter solar contribution.
- (b) The $\delta$-M-scaled single-scatter solar contribution
  (obtained by specializing Eq. (A7) with primed quantities) is
  subtracted off, leaving only the multiple-scatter piece of (a).
- (c) The un-truncated single-scatter solar contribution,
  evaluated using the full original phase function
  $P(\cos\Theta) = \sum_l \chi_l P_l(\cos\Theta)$ at the exact
  scattering angle $\Theta$, is added back in. The prefactor
  $1/(1 - f\omega)$ rescales the direct-beam attenuation from the
  $\delta$-M-scaled optical path used internally by the solver to
  the un-scaled physical solar optical path. This is the central
  insight of {cite:t}`nakajima1988tms`: keep $\delta$-M for the
  multiply-scattered field but restore the true angular dependence
  of single-scatter for radiance output.

The $\delta$-M + TMS pairing delivers radiance accuracy in the
forward-peak and aureole regions comparable to the un-scaled
solution with $L$ several times larger, which is why Evans reports
that SHDOM's radiance accuracy "is assured by using the delta-M
scaling method along with the un-truncated phase function
single-scattering solution (Nakajima and Tanaka, 1988)".
