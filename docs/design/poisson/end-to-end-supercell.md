(sec-showcase)=
# End-to-end example: supercell thunderstorm

The end-to-end example is a narrative demonstration that exercises the
complete diagnostic pipeline — dynamics, microphysics, the
electrification stub of [](governing-equations.md), the Poisson solve,
and the diagnostic output of $\boldsymbol{E}$ — inside a dynamically
evolving supercell run. The configuration is the standard
Klemp–Wilhelmson idealized supercell test case bundled with MPAS-A,
using the Kessler warm-rain microphysics suite. Since graupel and cloud
ice are absent in this suite, the stub uses rain water as the positive
proxy and cloud water as the negative proxy:

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

A calibrated 2 h dynamic sequence used $\alpha = \beta = 10^{-9}$ and
wrote electrostatic diagnostics every 30 min to a dedicated output
stream. The diagnostic figures show a 10 km half-width mean through the
storm core. Across the four plotted frames, the peak electric-field
magnitude is $4.17 \times 10^4$ V m$^{-1}$ at 60 min; the 2 h frame has
peak $|\boldsymbol{E}| = 2.59 \times 10^4$ V m$^{-1}$, maximum potential
$102$ MV relative to the grounded lower boundary, and charge-density
range $[-0.120,\,0.475]$ nC m$^{-3}$. The final PCG residuals remain
below $10^{-10}$ throughout the sequence.

```{figure} figures/tier_E_charge_coupled_030min.png
:name: fig-sc-030
:alt: Charge-coupled supercell diagnostic at 30 min.
:width: 100%

Charge-coupled supercell diagnostic at 30 min, shown as a 10 km
half-width mean through the storm core. The left panel shows liquid
water content with negative (red) and positive (black/gray)
charge-density contours; the right panel shows electrostatic potential
with electric-field streamlines.
```

```{figure} figures/tier_E_charge_coupled_060min.png
:name: fig-sc-060
:alt: Charge-coupled supercell diagnostic at 60 min.
:width: 100%

As in {numref}`fig-sc-030`, but at 60 min. This frame has the largest
electric-field magnitude in the plotted sequence,
$4.17 \times 10^4$ V m$^{-1}$.
```

```{figure} figures/tier_E_charge_coupled_090min.png
:name: fig-sc-090
:alt: Charge-coupled supercell diagnostic at 90 min.
:width: 100%

As in {numref}`fig-sc-030`, but at 90 min, after the initial core has
split into a broader charge and condensate structure.
```

```{figure} figures/tier_E_charge_coupled_120min.png
:name: fig-sc-120
:alt: Charge-coupled supercell diagnostic at 120 min.
:width: 100%

As in {numref}`fig-sc-030`, but at 120 min. The diagnostic sequence
shows that the one-way source, Poisson solve, and electric-field
reconstruction remain coupled to the evolving convective state.
```

The peak field magnitude is below common thunderstorm breakdown-scale
values, so the end-to-end example makes no predictive electrification
claim. The stub is hydrometeor-based diagnostic forcing, not a validated
noninductive charging scheme; the solve omits lightning discharge,
conductivity relaxation, ion screening, and charge leakage; and
potential is reported as a modeled potential difference relative to the
configured ground boundary. The purpose of the end-to-end example is to
show that the solver, diagnostic source path, time-sampled output, and
host dynamical core compose correctly on a realistic flow.
