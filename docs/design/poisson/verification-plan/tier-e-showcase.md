(tier-e)=
# Tier E — End-to-end example: supercell thunderstorm

Tier E is a narrative demonstration that exercises the complete
Phase 1 pipeline — dynamics, microphysics, the electrification
stub of [](../governing-equations.md), the Poisson solve, and the
diagnostic output of $\boldsymbol{E}$ — inside a realistic
end-to-end run. The configuration is the standard Klemp–Wilhelmson
idealized supercell test case bundled with MPAS-A, run for $12$
hours of simulated time with the electrification stub enabled and
the Poisson solver configured per the Phase 1 defaults of
[](../solver.md).

Four snapshot diagnostics are planned at two representative times
(mature storm, $t \approx 90$ min; dissipating storm,
$t \approx 4$ hr):

- Horizontal section of the charge density $\rho$ at $z = 6$ km
  altitude, showing the spatial structure of the separated charge
  produced by the stub.
- Vertical cross-section of the electrostatic potential $\varphi$
  through the storm updraft.
- Isosurface of $|\boldsymbol{E}| = 100$ kV m$^{-1}$ in three
  dimensions.
- Time series, over the full $12$-hour run, of the stored
  electrostatic energy
  $U_E(t) = \tfrac{\varepsilon_0}{2} \int |\boldsymbol{E}|^2 \, dV$,
  integrated over the domain.

Tier E makes no quantitative claim: the electrification stub is
not a validated electrification scheme, and the simulated fields
must not be read as a prediction of storm electrostatics. The
purpose of Tier E is narrative — to illustrate that the solver,
the stub, and the host dynamical core compose correctly on a
realistic flow and produce field structures physically consistent
with an electrified convective storm.
