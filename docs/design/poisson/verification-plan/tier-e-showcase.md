(tier-e)=
# Tier E — End-to-end example: supercell thunderstorm

Tier E is a narrative demonstration that exercises the complete Phase 1
diagnostic pipeline: dynamics, microphysics, the electrification stub
of [](../governing-equations.md), the Poisson solve, and the diagnostic
output of $\boldsymbol{E}$ inside a dynamically evolving supercell run.
The configuration is the standard Klemp-Wilhelmson idealized supercell
test case bundled with MPAS-A, using the Kessler warm-rain microphysics
suite. Since graupel and cloud ice are absent in this suite, the stub
uses rain water as the positive proxy and cloud water as the negative
proxy:

$$
\rho_{\mathrm{stub}}
  = \alpha\, w_{\mathrm{mid}}\, q_{\mathrm{r}}
    - \beta\, q_{\mathrm{c}},
$$ (eq-stub-kessler)

where $w_{\mathrm{mid}}$ is vertical velocity averaged from interfaces
to layer midpoints. The run is one-way coupled: the MPAS state fills
$\rho_{\mathrm{stub}}$, the Poisson solve writes $\varphi$ and
$\boldsymbol{E}$, and no electric feedback is applied to dynamics,
microphysics, chemistry, or lightning.

A calibrated 2 h dynamic sequence used
$\alpha = \beta = 10^{-9}$ and wrote electrostatic diagnostics every
30 min to the compact output stream. The diagnostic figures show a
10 km half-width mean through the storm core. Across the four plotted
frames, the peak electric-field magnitude is
$4.17 \times 10^4$ V m$^{-1}$ at 60 min; the 2 h frame has peak
$|\boldsymbol{E}| = 2.59 \times 10^4$ V m$^{-1}$, maximum potential
$102$ MV relative to the grounded lower boundary, and charge-density
range $[-0.120, 0.475]$ nC m$^{-3}$. The final
preconditioned conjugate-gradient (PCG) residuals remain below
$10^{-10}$ throughout the sequence.

| Time | Peak $|\boldsymbol{E}|$ | $\varphi$ range | $\rho_{\mathrm{stub}}$ range |
|---:|---:|---:|---:|
| 30 min | 39.6 kV m$^{-1}$ | 0.03 to 178.9 MV | -0.095 to 0.548 nC m$^{-3}$ |
| 60 min | 41.7 kV m$^{-1}$ | -0.34 to 167.0 MV | -0.114 to 0.632 nC m$^{-3}$ |
| 90 min | 24.0 kV m$^{-1}$ | -1.02 to 88.0 MV | -0.078 to 0.459 nC m$^{-3}$ |
| 120 min | 25.9 kV m$^{-1}$ | 0.05 to 102.4 MV | -0.120 to 0.475 nC m$^{-3}$ |

```{figure} ../figures/tier_E_charge_coupled_030min.png
:alt: Charge-coupled supercell diagnostic at 30 minutes.
:width: 100%

Tier E charge-coupled supercell diagnostic at 30 min, shown as a
10 km half-width mean through the storm core. The left panel shows
liquid water content with negative and positive charge-density
contours; the right panel shows electrostatic potential with
electric-field streamlines.
```

```{figure} ../figures/tier_E_charge_coupled_060min.png
:alt: Charge-coupled supercell diagnostic at 60 minutes.
:width: 100%

As above, but at 60 min. This frame has the largest electric-field
magnitude in the plotted sequence, $4.17 \times 10^4$ V m$^{-1}$.
```

```{figure} ../figures/tier_E_charge_coupled_090min.png
:alt: Charge-coupled supercell diagnostic at 90 minutes.
:width: 100%

As above, but at 90 min, after the initial core has split into a
broader charge and condensate structure.
```

```{figure} ../figures/tier_E_charge_coupled_120min.png
:alt: Charge-coupled supercell diagnostic at 120 minutes.
:width: 100%

As above, but at 120 min. The diagnostic sequence shows that the
one-way source, Poisson solve, and electric-field reconstruction remain
coupled to the evolving convective state.
```

Tier E makes no predictive electrification claim. The source is a
hydrometeor-based diagnostic stub, not a validated noninductive
charging scheme. The solve omits lightning discharge, conductivity
relaxation, ion screening, and charge leakage. Potential is reported as
a modeled potential difference relative to the configured ground
boundary. The purpose of Tier E is to show that the solver, diagnostic
source path, time-sampled output, and host dynamical core compose
correctly on a realistic flow.
