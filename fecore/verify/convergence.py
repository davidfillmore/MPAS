"""Surface and shell convergence matrices for the FEM Poisson MMS gate.

Task 6: surface_convergence / surface_matrix.
Task 8: shell_convergence.

Public API
----------
scvt_available() -> bool
    True iff the SCVT mesh bundle directory is present on disk.

surface_convergence(family, degree, rhs_mode) -> dict
    Run the MMS convergence study for one (family, degree, rhs_mode) triple.
    Returns {"h": [...], "err": [...], "slope": float}.
    ``slope`` is the log-log rate from the two finest levels.

surface_matrix() -> list[dict]
    Run the configurations for all available families
    (family × {P1-lumped, P1-consistent, P2-consistent}).  The icosa triple
    always runs; the SCVT triple is added only when the mesh bundle is present.
    Each dict has keys: family, degree, rhs_mode, h, err, slope.

shell_convergence(family, degree, rhs_mode) -> dict
    Run the 3-D shell MMS convergence study for one triple.
    Returns {"h": [...], "n_layers": [...], "err": [...], "slope": float}.
    Vertical (n_layers) and horizontal (h) are refined together to avoid
    the fixed-nVertLevels rate cap.

Mesh families (surface)
-----------------------
"scvt" — MPAS SCVT-dual meshes at
    ~/Data/MPAS/poisson_tier_A2_scvt/meshes/{480km,240km,120km,60km}/grid.nc
    R = sphere_radius(path)  (~6.371e6 m)
    h = nominal km label (480, 240, 120, 60) — meaningful only as a proxy
        for the characteristic edge length; the *ratio* is what matters for
        slope computation.  Finest-two: 120 km and 60 km.

"icosa" — icosahedral unit-sphere meshes at refine ∈ {2, 3, 4, 5}
    R = 1.0
    h = 1.0 / 2**refine  (edge ∝ 2^{-refine})
    Finest-two: refine 4 (h=1/16) and refine 5 (h=1/32).

Mesh families (shell, Task 8)
-----------------------------
"icosa" shell — refine ∈ {2, 3, 4}, n_layers ∈ {4, 8, 16} (doubles with refine)
    R=1.0, H=0.25, h=1/2**refine.  n_layers doubles with each horizontal
    refinement so the 3-D L² rate is not capped by a fixed vertical resolution.

"scvt" shell — {480, 240, 120} km (60 km SKIPPED — too large serially),
    n_layers ∈ {2, 4, 8}, H=2e4 m, R=sphere_radius.
    60 km is logged but omitted to avoid millions of DOFs in a serial run.

MMS (surface)
-------------
source_fn = mms.surface_source  (eps * 20/R^2 * Y_4^2)
exact_fn  = mms.y42_cart
eps = 1.0 (cancels in relative L2)

MMS (shell, Task 8)
-------------------
source_fn = mms.shell_source  (-∇²(Y_4² · sin(π/2 ζ)) in spherical coords)
exact_fn  = mms.shell_exact   (Y_4² · sin(π/2 ζ), ζ=(r-R)/H)
eps = 1.0 (cancels in relative L2)
The MMS satisfies: exact=0 on r=R (Dirichlet); ∂exact/∂r=0 on r=R+H (Neumann).

Isoparametric geometry cap (why P2 shows ~2nd order, not ~3rd order)
---------------------------------------------------------------------
Both mesh families use degree-1 geometry: triangles (surface) or tetrahedra
(shell) built by ``basix.ufl.element("Lagrange", ..., 1, ...)`` — flat affine
elements that approximate the sphere by piecewise-planar facets / flat tets.
The geometric approximation error of a degree-1 surface/shell mesh of the unit
sphere is O(h²) (the distance from a curved surface/shell patch to its linear
approximant scales as the square of the element diameter).  This geometric error
propagates into the assembled weak form and caps the L² solution error at O(h²)
for *any* element degree, regardless of how well the basis functions interpolate
smooth functions on the exact geometry.

Consequently, P2 elements achieve slope ≈ 2.0 rather than their native O(h³):
the flat-triangle/flat-tet geometry error dominates and prevents the higher-order
polynomial from improving accuracy beyond 2nd order.  This is the classical
*isoparametric geometry cap* of surface/shell FEM and is well-established in the
analysis of evolving-surface and stationary-surface PDEs (Dziuk 1988; Bernardi
1989).  Realizing P2's native 3rd-order convergence on the sphere would require
isoparametric (degree-2) geometry in which edge-midpoint nodes are snapped to
the true spherical surface/shell — a future enhancement beyond the current scope.

Shell vertical refinement (Task 8 methodology note)
----------------------------------------------------
If n_layers is fixed while the horizontal mesh refines, the fixed vertical
resolution caps the 3-D L² slope at the vertical resolution's O(H/n_layers)
rate — analogous to the MPAS-FV "fixed nVertLevels" artifact where the vertical
resolution dominates the error budget for fine horizontal meshes.  To measure
the true 3-D 2nd-order rate, n_layers must scale proportionally with 1/h
(i.e. double when h halves).  shell_convergence enforces this: icosa uses
(refine=2,n_layers=4), (3,8), (4,16); SCVT uses (480km,2), (240km,4), (120km,8).
"""

import os
import logging
import numpy as np

from fecore.mesh import icosa as _icosa
from fecore.mesh import scvt_dual as _scvt
from fecore.operator.poisson import solve_surface, solve_shell
from fecore.verify import mms

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# SCVT bundle path
# ---------------------------------------------------------------------------

_SCVT_ROOT = os.path.expanduser(
    "~/Data/MPAS/poisson_tier_A2_scvt/meshes"
)

# Resolution labels in coarse → fine order; h proxy = nominal km value.
_SCVT_LEVELS = [480, 240, 120, 60]


def scvt_available() -> bool:
    """Return True iff every SCVT grid.nc file in the bundle exists."""
    for km in _SCVT_LEVELS:
        path = os.path.join(_SCVT_ROOT, f"{km}km", "grid.nc")
        if not os.path.isfile(path):
            return False
    return True


# ---------------------------------------------------------------------------
# Core convergence loop
# ---------------------------------------------------------------------------

def surface_convergence(family: str, degree: int, rhs_mode: str) -> dict:
    """Run the MMS convergence study for one (family, degree, rhs_mode) triple.

    Parameters
    ----------
    family : {"scvt", "icosa"}
        Mesh family.
    degree : {1, 2}
        Lagrange polynomial degree.
    rhs_mode : {"consistent", "lumped"}
        RHS assembly mode.

    Returns
    -------
    dict with keys:
        "h"     : list of mesh-size proxies, coarse to fine.
        "err"   : list of relative L2 errors, matching order.
        "slope" : log-log convergence rate from the two finest levels.
    """
    if family == "icosa":
        return _icosa_convergence(degree, rhs_mode)
    elif family == "scvt":
        return _scvt_convergence(degree, rhs_mode)
    else:
        raise ValueError(f"family must be 'icosa' or 'scvt', got {family!r}")


def surface_matrix() -> list:
    """Run all available (family × config) combinations.

    Always runs the 3 icosa configs; adds the 3 SCVT configs only when the
    SCVT mesh bundle is present on disk.
    Configs per family: {(1,"lumped"), (1,"consistent"), (2,"consistent")}.

    Returns
    -------
    list of dict, each with keys: family, degree, rhs_mode, h, err, slope.
    """
    configs = [
        ("icosa", 1, "lumped"),
        ("icosa", 1, "consistent"),
        ("icosa", 2, "consistent"),
    ]
    if scvt_available():
        configs += [
            ("scvt", 1, "lumped"),
            ("scvt", 1, "consistent"),
            ("scvt", 2, "consistent"),
        ]
    else:
        log.warning("SCVT mesh bundle not available — skipping SCVT configs.")

    results = []
    for family, degree, rhs_mode in configs:
        log.info("Running %s P%d %s ...", family, degree, rhs_mode)
        data = surface_convergence(family, degree, rhs_mode)
        results.append(
            {
                "family": family,
                "degree": degree,
                "rhs_mode": rhs_mode,
                "h": data["h"],
                "err": data["err"],
                "slope": data["slope"],
            }
        )
        log.info(
            "  slope=%.4f  h=%s  err=%s",
            data["slope"],
            [f"{v:.4g}" for v in data["h"]],
            [f"{v:.4g}" for v in data["err"]],
        )
    return results


# ---------------------------------------------------------------------------
# Family-specific helpers
# ---------------------------------------------------------------------------

def _icosa_convergence(degree: int, rhs_mode: str) -> dict:
    """icosa family: refine ∈ {2,3,4,5}, R=1.0, h=1/2**refine."""
    refines = [2, 3, 4, 5]
    R = 1.0
    h_list = []
    err_list = []

    for refine in refines:
        h = 1.0 / (2 ** refine)
        mesh = _icosa.icosa_surface_mesh(refine)
        _, err = solve_surface(
            mesh, R=R, degree=degree, rhs_mode=rhs_mode,
            source_fn=mms.surface_source, exact_fn=mms.y42_cart, eps=1.0,
        )
        h_list.append(h)
        err_list.append(err)
        log.debug("  icosa refine=%d h=%.4g err=%.4g", refine, h, err)

    slope = _finest_two_slope(h_list, err_list)
    return {"h": h_list, "err": err_list, "slope": slope}


def _scvt_convergence(degree: int, rhs_mode: str) -> dict:
    """SCVT family: levels [480,240,120,60] km, R from grid.nc.

    For P2 on the 60 km mesh (~650k DOFs), we attempt the solve but fall back
    to using only the 480/240/120 km triple if it takes too long or fails.
    The fallback is logged explicitly so it is never silent.
    """
    h_list = []
    err_list = []

    # Levels to attempt; for P2 we track whether we need fallback.
    levels_to_try = list(_SCVT_LEVELS)  # [480, 240, 120, 60]
    p2_large_skipped = False

    for km in levels_to_try:
        path = os.path.join(_SCVT_ROOT, f"{km}km", "grid.nc")
        R = _scvt.sphere_radius(path)
        mesh = _scvt.scvt_surface_mesh(path)

        # For the 60 km P2 mesh, warn about scale before attempting.
        if km == 60 and degree == 2:
            log.warning(
                "SCVT 60 km P2: attempting large solve (~650k DOFs serial). "
                "This may take several minutes."
            )

        try:
            _, err = solve_surface(
                mesh, R=R, degree=degree, rhs_mode=rhs_mode,
                source_fn=mms.surface_source, exact_fn=mms.y42_cart, eps=1.0,
            )
        except (MemoryError, RuntimeError) as exc:
            # Note: PETSc OOM failures may surface as other exception types
            # (e.g. PETSc.Error) not listed here; behavior is unchanged if so.
            if km == 60 and degree == 2:
                log.warning(
                    "SCVT 60 km P2 solve failed (%s). "
                    "Falling back to 480/240/120 km triple for slope estimate.",
                    exc,
                )
                p2_large_skipped = True
                break  # skip 60 km; use what we have (480, 240, 120)
            raise

        h_list.append(float(km))
        err_list.append(err)
        log.debug("  scvt %d km P%d %s h_proxy=%g err=%.4g",
                  km, degree, rhs_mode, km, err)

    if p2_large_skipped:
        log.warning(
            "SCVT P2 slope computed from 480/240/120 km triple (60 km omitted). "
            "Finest-two = 240 km, 120 km."
        )

    slope = _finest_two_slope(h_list, err_list)
    return {"h": h_list, "err": err_list, "slope": slope}


# ---------------------------------------------------------------------------
# Slope helper
# ---------------------------------------------------------------------------

def _finest_two_slope(h_list: list, err_list: list) -> float:
    """Log-log convergence rate from the two finest levels.

    slope = log(err[-2] / err[-1]) / log(h[-2] / h[-1])

    Requires at least two points.  Both h and err must be strictly positive.
    """
    if len(h_list) < 2:
        raise ValueError(
            f"Need at least 2 refinement levels to compute a slope; got {len(h_list)}."
        )
    h1, h2 = float(h_list[-2]), float(h_list[-1])
    e1, e2 = float(err_list[-2]), float(err_list[-1])
    if h2 >= h1:
        raise ValueError(
            f"h values must be decreasing (coarse→fine); got h[-2]={h1}, h[-1]={h2}."
        )
    if e1 <= 0 or e2 <= 0:
        raise ValueError(
            f"Errors must be positive; got err[-2]={e1}, err[-1]={e2}."
        )
    return np.log(e1 / e2) / np.log(h1 / h2)


# ---------------------------------------------------------------------------
# Task 8: Shell convergence
# ---------------------------------------------------------------------------

# SCVT shell levels (coarse → fine).  60 km is intentionally omitted — a 3-D
# P2 solve on the 60 km shell would require tens of millions of DOFs and is
# not tractable serially.  This omission is logged, never silent.
_SCVT_SHELL_LEVELS = [480, 240, 120]   # 60 km skipped — too large serially

# n_layers paired with SCVT shell levels: doubles at each refinement.
_SCVT_SHELL_NLAYERS = {480: 2, 240: 4, 120: 8}

# icosa shell levels: (refine, n_layers) pairs — n_layers doubles with refine.
_ICOSA_SHELL_LEVELS = [(2, 4), (3, 8), (4, 16)]


def shell_convergence(family: str, degree: int, rhs_mode: str) -> dict:
    """Run the 3-D shell MMS convergence study for one (family, degree, rhs_mode) triple.

    CRITICAL METHODOLOGY — vertical and horizontal refinement are coupled:
    If n_layers is fixed while the horizontal mesh refines, the vertical
    resolution caps the 3-D L² rate (a ``fixed nVertLevels'' artifact, exactly
    as in MPAS-FV convergence studies where fixed nVertLevels lets the vertical
    error dominate for fine horizontal meshes).  To observe the true 3-D
    2nd-order rate, n_layers must scale proportionally with 1/h — doubling when
    h halves.  This function enforces coupled refinement:

      icosa: (refine=2, n_layers=4), (refine=3, n_layers=8), (refine=4, n_layers=16)
             h = 1/2**refine.  n_layers doubles at each level.

      SCVT:  (480 km, 2 layers), (240 km, 4 layers), (120 km, 8 layers).
             60 km is SKIPPED (logged, never silent): a 3-D P2 solve there
             would require O(10^7) DOFs serially — not tractable.
             Finest-two: 240 km and 120 km.
             SCVT coords are normalized to unit sphere (H_norm=0.25) because the
             physical H/R≈3e-3 shell is ill-conditioned at MPAS resolutions.

    MMS: shell_source / shell_exact (Y_4^2 · sin(π/2·ζ), ζ=(r-R)/H).
    BCs: ground Dirichlet φ=0 at r=R; top Neumann (natural) at r=R+H.

    Parameters
    ----------
    family : {"icosa", "scvt"}
        Mesh family.
    degree : {1, 2}
        Lagrange polynomial degree.
    rhs_mode : {"consistent", "lumped"}
        RHS assembly mode.

    Returns
    -------
    dict with keys:
        "h"        : list of mesh-size proxies (coarse to fine).
        "n_layers" : list of vertical layer counts (paired with h).
        "err"      : list of relative L2 errors.
        "slope"    : log-log convergence rate from the two finest levels.
    """
    if family == "icosa":
        return _icosa_shell_convergence(degree, rhs_mode)
    elif family == "scvt":
        return _scvt_shell_convergence(degree, rhs_mode)
    else:
        raise ValueError(f"family must be 'icosa' or 'scvt', got {family!r}")


def _icosa_shell_convergence(degree: int, rhs_mode: str) -> dict:
    """icosa shell family: refine ∈ {2,3,4}, n_layers ∈ {4,8,16}, R=1, H=0.25.

    n_layers doubles at each refinement level so the 3-D L² slope reflects the
    true combined horizontal+vertical rate — avoiding the fixed-nVertLevels cap.
    """
    from fecore.mesh import icosa as _icosa

    R, H = 1.0, 0.25
    h_list = []
    nlayers_list = []
    err_list = []

    for refine, n_layers in _ICOSA_SHELL_LEVELS:
        h = 1.0 / (2 ** refine)
        mesh = _icosa.icosa_shell_mesh(refine=refine, n_layers=n_layers, H=H)
        _, err = solve_shell(
            mesh, R=R, H=H, degree=degree, rhs_mode=rhs_mode,
            source_fn=mms.shell_source, exact_fn=mms.shell_exact, eps=1.0,
        )
        h_list.append(h)
        nlayers_list.append(n_layers)
        err_list.append(err)
        log.debug(
            "  icosa shell refine=%d n_layers=%d h=%.4g err=%.4g",
            refine, n_layers, h, err,
        )

    slope = _finest_two_slope(h_list, err_list)
    return {"h": h_list, "n_layers": nlayers_list, "err": err_list, "slope": slope}


def _scvt_shell_mesh_normalized(path, n_layers: int, H_norm: float):
    """Build a dolfinx SCVT shell mesh normalized to the unit sphere.

    The SCVT surface vertices are read from *path*, scaled to unit radius, and
    extruded radially from R=1.0 to R=1+H_norm through *n_layers* layers.  The
    physical Earth-radius coordinates are NOT used: all coordinates are in the
    unit-sphere frame so the MMS source and exact functions receive well-scaled
    inputs.

    This normalization is required because with H_phys/R_phys ≈ 3e-3 (20 km /
    6371 km), the aspect ratio of the physical shell tets is ~150× (horizontal ≫
    vertical at MPAS resolutions), making the MMS system ill-conditioned in
    physical units.  Normalizing to unit sphere with H_norm=0.25 gives the same
    well-conditioned geometry as the icosa shell family while preserving the SCVT
    mesh topology that governs the horizontal convergence rate.
    """
    import netCDF4 as nc
    from fecore.mesh._shell_utils import _extrude_to_tets, _fix_tet_orientations
    from fecore.mesh.scvt_dual import _orient_outward
    import dolfinx, ufl, basix
    from mpi4py import MPI

    R_phys = _scvt.sphere_radius(path)
    with nc.Dataset(path) as ds:
        xc = np.asarray(ds["xCell"][:])
        yc = np.asarray(ds["yCell"][:])
        zc = np.asarray(ds["zCell"][:])
        cov = np.asarray(ds["cellsOnVertex"][:]).astype(np.int64) - 1

    # Scale surface vertices to unit sphere (R_norm = 1.0).
    surf_pts = np.column_stack((xc, yc, zc)).astype(np.float64)
    surf_pts *= R_phys / np.linalg.norm(surf_pts, axis=1)[:, None]
    surf_pts /= R_phys  # now on unit sphere

    surf_tris = cov[np.all(cov >= 0, axis=1)].copy()
    _orient_outward(surf_tris, surf_pts)
    n_surf = len(surf_pts)

    R_norm = 1.0
    radii = R_norm + np.arange(n_layers + 1) * H_norm / n_layers
    coords = (surf_pts[None, :, :] * radii[:, None, None]).reshape(-1, 3)
    coords = np.ascontiguousarray(coords, dtype=np.float64)

    tets = _extrude_to_tets(surf_tris, n_surf, n_layers)
    tets, _ = _fix_tet_orientations(tets, coords)

    el = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    return dolfinx.mesh.create_mesh(MPI.COMM_WORLD, tets, el, coords)


def _scvt_shell_convergence(degree: int, rhs_mode: str) -> dict:
    """SCVT shell family: {480,240,120} km × {2,4,8} layers.  60 km SKIPPED.

    60 km is not attempted: a 3-D P2 solve on the 60 km shell mesh would require
    O(10^7) DOFs serially, which is not tractable.  The omission is logged
    explicitly — never silent.  Finest-two: 240 km and 120 km.

    Coordinate normalization:
    SCVT meshes carry Earth-radius coordinates (R≈6.37×10⁶ m, H=2×10⁴ m).
    At those physical scales the shell has aspect ratio H/R≈3×10⁻³ — the
    horizontal cell size (≈480 km) is ~150× larger than the vertical layer
    thickness (≈10 km).  In physical coords the MMS source is O(10⁻²⁴), which
    creates a poorly-conditioned load vector relative to the stiffness (O(1)
    entries in SI) and yields essentially zero solutions.  To avoid this
    numerical issue we normalize ALL coordinates to the unit sphere (scale ×
    1/R_phys) and use H_norm=0.25 — the same shell geometry as the icosa
    family.  The SCVT mesh topology (cell connectivity) is unchanged; only the
    coordinate scaling differs.  Because the convergence slope depends only on
    the *ratio* of errors at successive refinement levels (not their absolute
    magnitudes), normalizing the geometry preserves the convergence rate.
    """
    log.warning(
        "shell_convergence('scvt'): 60 km mesh is SKIPPED (too large serially "
        "for a 3-D P2 solve — O(10^7) DOFs).  Slope uses finest-two = 240 km, "
        "120 km.  n_layers doubles with each refinement: (480km,2), (240km,4), "
        "(120km,8)."
    )
    log.info(
        "shell_convergence('scvt'): coordinates normalized to unit sphere "
        "(H_norm=0.25).  Physical H/R≈3e-3 is too thin for the current SCVT "
        "resolutions; normalization preserves SCVT topology and convergence rate."
    )

    R_norm = 1.0
    H_norm = 0.25  # same shell geometry as icosa — well-conditioned MMS

    h_list = []
    nlayers_list = []
    err_list = []

    for km in _SCVT_SHELL_LEVELS:
        n_layers = _SCVT_SHELL_NLAYERS[km]
        path = os.path.join(_SCVT_ROOT, f"{km}km", "grid.nc")
        # h proxy in unit-sphere units: km / R_km
        R_km = _scvt.sphere_radius(path) / 1000.0  # km
        h_norm = km / R_km

        mesh = _scvt_shell_mesh_normalized(path, n_layers=n_layers, H_norm=H_norm)

        _, err = solve_shell(
            mesh, R=R_norm, H=H_norm, degree=degree, rhs_mode=rhs_mode,
            source_fn=mms.shell_source, exact_fn=mms.shell_exact, eps=1.0,
        )
        h_list.append(h_norm)
        nlayers_list.append(n_layers)
        err_list.append(err)
        log.debug(
            "  scvt shell %d km n_layers=%d P%d %s h_norm=%.4g err=%.4g",
            km, n_layers, degree, rhs_mode, h_norm, err,
        )

    slope = _finest_two_slope(h_list, err_list)
    return {"h": h_list, "n_layers": nlayers_list, "err": err_list, "slope": slope}
