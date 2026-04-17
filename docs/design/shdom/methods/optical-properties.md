(sec-optical-properties)=
# Optical-property pipeline via RRTMG

SHDOM requires, for each Cartesian cell and each active band, the
triple $(\beta, \omega, g)$ of extinction, single-scattering
albedo, and asymmetry. RRTMG computes the physical quantities
underlying this triple on its way through each column's solve, but
does not expose them: its public entry points return fluxes and
heating rates, not the intermediate optical properties, and its
band-internal $g$-point loop is entangled with the Monte Carlo
independent-column draw over sub-column cloud realizations
{cite:p}`pincus2003mcica`.

The framework solves this by *addition* rather than modification. A
new entry point `rrtmg_sw_compute_optical_properties` is appended
to the existing RRTMG shortwave module. The new routine replays
the band loop and the cell loop of the standard driver, calls the
unchanged cloud-optical-property kernel (`cldprmc_sw`) and
gas-absorption kernel (`taumol_sw`) and aerosol contribution,
skips the radiative transfer solver, and returns the assembled
optical properties to the caller. No line of the standard RRTMG
solver path is edited. The RRTMG subtree remains byte-identical to
its upstream release except for the appended new entry point, and
future upstream resyncs apply as ordinary merges. The exact
line-range footprint is documented in
[](../appendices/rrtmg-patches.md).

Two design choices inside the new entry point bear explicit
mention. First, the MCICA stochastic sub-column generator is
bypassed: SHDOM operates on the explicit 3D cloud field and has no
need for the stochastic overlap sampling that MCICA provides for a
1D column. The new routine calls `cldprmc_sw` on the deterministic
cloud fraction and condensate fields directly. This has the side
effect of making SHDOM's optical-property inputs exactly the
fields the dynamical core sees, which keeps SHDOM-vs-RRTMG
comparisons physically interpretable (both solvers see the same
instantaneous cloud). It also means that the RRTMG-internal
sub-column averaging and the SHDOM-internal 3D transport operate
on the same cloud representation, so any difference between their
outputs reflects only the 3D transport physics and not sampling
noise. Second, the output is *component-decomposed*: the routine
returns, per cell per level per band, nine scalars — the three
quantities $\beta$, $\beta_s = \beta\omega$, and $\beta_s g$ for
each of the three components (gas absorption, cloud, aerosol), and
not the combined $(\beta, \omega, g)$. The rationale is twofold.
Mixing the components destroys the information needed to apply the
bound-preserving combination of Eq. {eq}`eq-v2c-c2v` on the target
mesh; conversely, transporting the three intensive components
through the area-weighted map and combining afterwards preserves
the physical bounds by construction (see
[](mesh-mapping.md)). The combination step itself is deferred to
the SHDOM interface on the Cartesian side, where the combined
optical properties are consumed immediately by the solver.

A second add-only routine,
`rrtmg_sw_compute_perband_fluxes`, exposes RRTMG's per-band fluxes
for the purpose of the energy-budget comparison in the eventual
validation programme. RRTMG's standard driver sums band-by-band
fluxes into broadband totals and does not retain the per-band
decomposition; the validation plan requires the per-band RRTMG
flux in order to compare like with like against SHDOM's
single-visible-band output. The routine shares the band-loop
scaffolding with the existing driver but preserves the per-band
flux contributions; it is again appended at the end of the module
with no edits to existing lines.
