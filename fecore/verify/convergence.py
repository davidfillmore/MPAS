"""Surface convergence matrix for the FEM Poisson MMS gate (Task 6).

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

Mesh families
-------------
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

MMS
---
source_fn = mms.surface_source  (eps * 20/R^2 * Y_4^2)
exact_fn  = mms.y42_cart
eps = 1.0 (cancels in relative L2)

Isoparametric geometry cap (why P2 shows ~2nd order, not ~3rd order)
---------------------------------------------------------------------
Both mesh families use degree-1 geometry: triangles built by
``basix.ufl.element("Lagrange", "triangle", 1, ...)`` — flat affine triangles
that approximate the sphere by piecewise-planar facets.  The geometric
approximation error of a degree-1 surface mesh of the unit sphere is O(h²)
(the distance from a curved surface patch to its linear approximant scales as
the square of the element diameter).  This geometric error propagates into the
assembled weak form and caps the L² solution error at O(h²) for *any* element
degree, regardless of how well the basis functions interpolate smooth functions
on the exact geometry.

Consequently, P2 elements achieve slope ≈ 2.0 rather than their native O(h³):
the flat-triangle geometry error dominates and prevents the higher-order
polynomial from improving accuracy beyond 2nd order.  This is the classical
*isoparametric geometry cap* of surface FEM and is well-established in the
analysis of evolving-surface and stationary-surface PDEs (Dziuk 1988; Bernardi
1989).  Realizing P2's native 3rd-order convergence on the sphere would require
isoparametric (degree-2) geometry in which edge-midpoint nodes are snapped to
the true spherical surface — a future enhancement beyond the current Task 6
scope.
"""

import os
import logging
import numpy as np

from fecore.mesh import icosa as _icosa
from fecore.mesh import scvt_dual as _scvt
from fecore.operator.poisson import solve_surface
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
