(sec-framework)=
# Hybrid 1D/3D RT framework

The framework is assembled from four single-purpose modules that
sit between the MPAS-A physics driver and the existing RRTMG entry
points. The *dispatcher* orchestrates a single radiation step: it
invokes RRTMG as it would ordinarily, and, when SHDOM is enabled,
also invokes the 3D solver on a designated subset of cells for a
designated subset of bands. The *optical-property extractor*
exposes RRTMG's internal $(\beta, \beta_s, \beta_s g)$ by component
(gas, cloud, aerosol) without disturbing the SW/LW solver path
(see [](optical-properties.md)). The *mesh mapper* carries those
optical properties between MPAS-A's Voronoi cells and SHDOM's
Cartesian grid and returns SHDOM's fluxes back to the unstructured
side (see [](mesh-mapping.md)). The *SHDOM interface* mixes the
component-decomposed optical properties into SHDOM's combined
$(\beta, \omega, g)$ representation on the Cartesian side, calls
SHDOM, and exposes the returned fluxes and heating rates with a
clean Fortran-array API that contains no MPAS-A types.

## Architectural invariants

Four invariants govern how these modules compose across all phases
of the roadmap.

### Authoritative solver per (cell, band)

The dynamics-facing heating rate is, for every cell and every band,
the contribution of whichever solver has been designated
authoritative for that (cell, band) pair. In Phase 1 that
designation is RRTMG everywhere, and SHDOM runs in a purely
diagnostic mode: SHDOM-derived heating rates and fluxes are written
to output but are not fed back to the dynamical core. Later phases
flip individual (cell, band) entries from RRTMG-authoritative to
SHDOM-authoritative by adjusting a mask; no structural change to
the dispatcher is required, and the per-band diagnostic outputs of
both solvers remain available for cross-comparison throughout.

### Domain boundaries use 1D RT

SHDOM is never asked to resolve a lateral domain boundary. On the
doubly-periodic supercell configuration of Phase 1 no such boundary
exists. On later configurations in which SHDOM runs on a
sub-domain of an MPAS-A mesh, the sub-domain is necessarily
surrounded by cells running 1D RT, and the horizontal boundary of
the SHDOM Cartesian box is an *interior* 3D/1D interface — never a
domain edge. The consequence is that the only SHDOM boundary
condition that needs to be designed within the framework is the
interior 3D/1D interface (deferred to Phase 3); the periodic
boundary of Phase 1 is handled natively by SHDOM.

### Diagnostic SHDOM outputs exist wherever SHDOM runs

Even after a (cell, band) pair becomes SHDOM-authoritative in a
later phase, the corresponding RRTMG output continues to be
written, and vice versa. Cross-solver comparison therefore remains
available as an ongoing validation channel rather than being closed
off by the authority flip.

### Clean module boundaries

The mesh mapper has no knowledge of SHDOM or of RRTMG; the SHDOM
interface has no MPAS-A types in its public signature; the
optical-property extractor keeps the details of the RRTMG patch out
of the dispatcher. Each module can be unit-tested in isolation
against its interface, and in particular the mesh mapper — the
framework's central numerical contribution — is exercised against
synthetic fields on a standalone test driver (see
[](mesh-mapping.md)) without ever invoking RRTMG or SHDOM.

## Phase 1 dispatch flow

The Phase 1 shortwave dispatch flow, expressed in pseudocode, is:

```text
rad_dispatcher_sw:
  1. rrtmg_swrad(all_cells, all_bands) -> heating_dynamics
  2. if shdom_enabled and time_for_shdom:
       props = optical_props_rrtmg(cells_3d_mask, shdom_bands)
       exchange_halo(props)
       cart_props  = mesh_map.v2c(props, cartesian_box)
       cart_fluxes = shdom_interface(cart_props, bc, solar_geom)
       v_fluxes    = mesh_map.c2v(cart_fluxes, cells_3d_mask)
       heating_diag_shdom[cells_3d_mask, shdom_bands]
           = fluxdiv_to_heating(v_fluxes)
```

The two optional inputs — `cells_3d_mask` and `shdom_bands` — carry
trivial values in Phase 1 (all cells; a single visible band) and
become the principal configuration knobs in the later phases,
without any change to the dispatch structure itself.

Two runtime controls select the pipeline: setting
`config_radt_sw_scheme = 'shdom_hybrid'` activates the dispatcher
in place of the bare RRTMG call, and `config_shdom_cadence_secs`
throttles SHDOM to run only at intervals longer than the physics
timestep. When the cadence is longer than a single radiation step,
held steps expose the last-written SHDOM diagnostic fluxes and
record a flag; the dynamics-facing heating rate, served from
RRTMG, is unaffected. The full namelist surface — cell-mask
selection, Cartesian-box dimensions, boundary-condition mode, band
selection — is specified in the companion design document
available alongside the source code.
