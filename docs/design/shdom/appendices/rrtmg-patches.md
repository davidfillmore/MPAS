(app-rrtmg-patches)=
# RRTMG patches

Two add-only routines are appended to the bundled RRTMG shortwave
module `module_ra_rrtmg_sw.F`. The first,
`rrtmg_sw_compute_optical_properties`, duplicates the
cell-and-band scaffolding of the existing `rrtmg_swrad` driver,
calls the unchanged cloud-optical-property kernel `cldprmc_sw`
and gas-absorption kernel `taumol_sw` and aerosol contribution,
skips the solver path (`reftra_sw` and the remainder of
`spcvmc_sw`), and returns the assembled optical properties as
component-decomposed $(\beta_{\mathrm{gas}},
\beta_{\mathrm{cloud}}, \omega_{\mathrm{cloud}},
g_{\mathrm{cloud}}, \beta_{\mathrm{aero}}, \omega_{\mathrm{aero}},
g_{\mathrm{aero}})$ arrays for each cell, level, and band
requested. The second, `rrtmg_sw_compute_perband_fluxes`,
similarly duplicates the band-loop scaffolding and exposes the
up/down/direct fluxes per band that the standard driver otherwise
sums into broadband totals.

Three invariants hold on the patch.

1. The existing RRTMG routines — `rrtmg_swrad`, `rrtmg_sw`,
   `spcvmc_sw`, `reftra_sw`, `cldprmc_sw`, `taumol_sw`,
   `module_ra_rrtmg_sw_aerosols` — are byte-identical to their
   upstream form. The patch adds new routines in an appended block
   at end-of-file and makes no edits to any existing line.
   Verification is by `git diff` against the vendored upstream
   baseline.
2. The new optical-property routine bypasses the MCICA sub-column
   generator, which is appropriate because SHDOM has the explicit
   3D cloud field and does not need stochastic overlap sampling.
   The deterministic cloud field (`cldfrac`, `clwp`, `rei`, `rel`)
   is passed directly into `cldprmc_sw`.
3. The per-band flux routine preserves the per-band contributions
   that the standard broadband-only driver discards; it does not
   otherwise alter the radiative transfer.

The line-range footprint of both patches, and their `git log`
trails against the vendored baseline, are recorded in `BUILD.md`
of the MPAS fork; the exact byte boundaries are omitted here as
they will drift with upstream RRTMG releases.
