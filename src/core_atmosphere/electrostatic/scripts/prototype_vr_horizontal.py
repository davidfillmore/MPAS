#!/usr/bin/env python3
"""Variable-resolution horizontal characterization: baseline vs cotangent (Task B4).

Axis 2 = Tier C: compare the baseline diagonal-Hodge operator and the cotangent
(DEC/FEM) operator on perturbed (non-well-centred) planar Voronoi meshes.

Measure: SOLUTION ERROR for both operators (fair common yardstick).
  - Baseline FV operator: solution-error order ≈ operator-residual order.
  - Cotangent operator: consistent in SOLUTION error (~2nd order) but NOT in
    pointwise/operator-residual terms; solution error is the honest measure.

Run:
    ~/miniconda3/envs/mpas/bin/python -m unittest \
        src.core_atmosphere.electrostatic.tests.test_vr_horizontal -v
"""
from __future__ import annotations

import importlib.util
import math
import pathlib
from types import SimpleNamespace

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import Delaunay, Voronoi

_MFD_PATH = pathlib.Path(__file__).resolve().parent / "mfd_operator.py"


def _mfd():
    """Load mfd_operator module by file path (avoids package-path ambiguity)."""
    spec = importlib.util.spec_from_file_location("mfd_operator", _MFD_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Re-export spd_min_eig so tests can call m.spd_min_eig(...)
def spd_min_eig(A_csr):
    """Smallest algebraic eigenvalue — delegates to mfd_operator.spd_min_eig."""
    return _mfd().spd_min_eig(A_csr)


# ---------------------------------------------------------------------------
# Mesh builder
# ---------------------------------------------------------------------------

def perturbed_mesh(n_side, perturb_frac, seed=0, L=1.0):
    """Build a perturbed planar Voronoi mesh on [0,L]^2.

    Interior points from an (n_side x n_side) square lattice are displaced by
    Gaussian noise with std = perturb_frac * spacing.  Boundary ring is fixed.
    Deterministic via numpy.random.default_rng(seed).

    Returns SimpleNamespace with:
        xy          (n, 2) point coordinates
        tris        (T, 3) Delaunay triangle vertex indices
        edge_pairs  list of (a, b) unique undirected edges (canonical a < b)
        dcEdge      dict (a,b) -> primal edge length
        dvEdge      dict (a,b) -> dual (Voronoi) edge length (0 for boundary)
        area        (n,) Voronoi cell areas (NaN-filled for open boundary cells)
        interior    (n,) bool mask: True for strictly interior cells
        L, n, n_side
    """
    dx = L / n_side
    xs, ys = np.meshgrid(
        (np.arange(n_side) + 0.5) * dx,
        (np.arange(n_side) + 0.5) * dx,
    )
    pts = np.column_stack((xs.ravel(), ys.ravel()))

    # Boundary ring: outermost ring of lattice points
    near_edge = (
        (pts[:, 0] < dx) | (pts[:, 0] > L - dx) |
        (pts[:, 1] < dx) | (pts[:, 1] > L - dx)
    )

    rng = np.random.default_rng(seed)
    jitter = rng.normal(scale=perturb_frac * dx, size=pts.shape)
    pts = pts.copy()
    pts[~near_edge] += jitter[~near_edge]

    # Delaunay triangulation (used to assemble cotangent Hodge)
    tri = Delaunay(pts)
    tris = tri.simplices  # (T, 3) int

    # Voronoi diagram (used for dcEdge / dvEdge / cell areas)
    vor = Voronoi(pts)

    # Build edge dictionaries
    dc = {}
    dv = {}
    edge_pairs = []
    for (a, b), ridge_verts in zip(vor.ridge_points, vor.ridge_vertices):
        a, b = int(a), int(b)
        key = (min(a, b), max(a, b))
        dc[key] = float(np.linalg.norm(pts[a] - pts[b]))
        if -1 in ridge_verts:
            dv[key] = 0.0  # open (boundary) Voronoi edge
        else:
            v0 = vor.vertices[ridge_verts[0]]
            v1 = vor.vertices[ridge_verts[1]]
            dv[key] = float(np.linalg.norm(v0 - v1))
        edge_pairs.append(key)

    # Voronoi cell areas (NaN for open/boundary cells)
    area = np.full(len(pts), np.nan)
    for p_idx, reg_idx in enumerate(vor.point_region):
        verts = vor.regions[reg_idx]
        if not verts or -1 in verts:
            continue
        poly = vor.vertices[verts]
        x_p, y_p = poly[:, 0], poly[:, 1]
        area[p_idx] = 0.5 * abs(
            np.dot(x_p, np.roll(y_p, 1)) - np.dot(y_p, np.roll(x_p, 1))
        )

    interior = ~near_edge

    return SimpleNamespace(
        xy=pts, tris=tris, edge_pairs=edge_pairs,
        dcEdge=dc, dvEdge=dv, area=area, interior=interior,
        L=L, n=len(pts), n_side=n_side,
    )


# ---------------------------------------------------------------------------
# Incidence matrix
# ---------------------------------------------------------------------------

def _incidence(mesh):
    """Signed incidence d0: (n_edges x n_nodes), orientation a->b."""
    n_edges = len(mesh.edge_pairs)
    rows, cols, data = [], [], []
    for e, (a, b) in enumerate(mesh.edge_pairs):
        rows += [e, e]
        cols += [a, b]
        data += [-1.0, 1.0]
    return sp.csr_matrix((data, (rows, cols)), shape=(n_edges, mesh.n))


# ---------------------------------------------------------------------------
# Operator assembly
# ---------------------------------------------------------------------------

def vr_operator(mesh, variant):
    """Assemble the scalar Poisson operator A = d0^T diag(w) d0.

    variant="baseline"  : w_e = dvEdge[e] / dcEdge[e]   (TPFA / diagonal Hodge)
    variant="cotangent" : w_e = sum of cotangent half-weights over triangles
                          sharing edge e.  Uses mfd_operator.cotangent_hodge_weights
                          (returns (weights, local_pairs) — NOT a dict).

    Returns ungrounded symmetric CSR matrix (annihilates constants in interior).
    """
    d0 = _incidence(mesh)
    n_edges = len(mesh.edge_pairs)
    w = np.zeros(n_edges)

    # Map canonical edge key -> index in edge_pairs
    edge_index = {pair: e for e, pair in enumerate(mesh.edge_pairs)}

    if variant == "baseline":
        for e, key in enumerate(mesh.edge_pairs):
            dc_e = mesh.dcEdge[key]
            dv_e = mesh.dvEdge[key]
            w[e] = (dv_e / dc_e) if dc_e > 1e-300 else 0.0

    elif variant == "cotangent":
        mfd = _mfd()
        xyz = np.column_stack((mesh.xy, np.zeros(mesh.n)))  # embed z=0
        for tri_verts in mesh.tris:
            tri_xyz = xyz[tri_verts]  # (3,3)
            weights, local_pairs = mfd.cotangent_hodge_weights(tri_xyz)
            # local_pairs = [(0,1),(1,2),(0,2)]; weights[m] = half-cot of opposite angle
            for m, (li, lj) in enumerate(local_pairs):
                gi = int(tri_verts[li])
                gj = int(tri_verts[lj])
                key = (min(gi, gj), max(gi, gj))
                if key in edge_index:
                    w[edge_index[key]] += weights[m]

    else:
        raise ValueError(f"Unknown variant: {variant!r} (expected 'baseline' or 'cotangent')")

    A = (d0.T @ sp.diags(w) @ d0).tocsr()
    return 0.5 * (A + A.T)  # symmetrize to kill floating-point skew


# ---------------------------------------------------------------------------
# Solution-error MMS
# ---------------------------------------------------------------------------

def _solution_error(mesh, variant):
    """Relative L2 solution error for phi_exact = sin(pi x/L) sin(pi y/L).

    Recipe:
      - Assemble A = vr_operator(mesh, variant)
      - Lumped load: b_i = 2*(pi/L)^2 * phi_exact_i * area_i  (area from Voronoi)
      - Pin boundary cells to phi_exact; solve interior system with spsolve
      - Return sqrt(sum_I area*(phi_num - phi_exact)^2) / sqrt(sum_I area*phi_exact^2)
    """
    A = vr_operator(mesh, variant)
    x, y = mesh.xy[:, 0], mesh.xy[:, 1]
    k = math.pi / mesh.L
    phi_exact = np.sin(k * x) * np.sin(k * y)

    # Cell areas: use Voronoi areas; fill NaN (boundary/open cells) with median of interior
    area = mesh.area.copy()
    int_mask = mesh.interior
    median_area = np.nanmedian(area[int_mask])
    area[np.isnan(area)] = median_area

    # Clamp negative/zero areas (degenerate slivers) to median
    area = np.where(area > 0, area, median_area)

    # Lumped load
    b = 2.0 * (k ** 2) * phi_exact * area

    # Split into interior (I) and boundary (bd) sets
    I = np.where(int_mask)[0]
    bd = np.where(~int_mask)[0]

    # RHS for interior dofs: b[I] - A[I, bd] @ phi_exact[bd]
    A_II = A[np.ix_(I, I)]
    A_Ibd = A[np.ix_(I, bd)]
    rhs = b[I] - A_Ibd @ phi_exact[bd]

    phi_num_I = spla.spsolve(A_II.tocsc(), rhs)

    phi_num = phi_exact.copy()
    phi_num[I] = phi_num_I

    # Relative L2 error over interior cells
    err2 = float(np.sum(area[I] * (phi_num[I] - phi_exact[I]) ** 2))
    ref2 = float(np.sum(area[I] * phi_exact[I] ** 2))
    return math.sqrt(err2 / ref2) if ref2 > 0 else float("inf")


def vr_solution_error(perturb_frac, variant, sides=(16, 24, 36), seed=0, L=1.0):
    """Compute finest-two log-log convergence slope in solution error.

    Runs the MMS solution error at each mesh size in `sides` and returns
    the slope log(err[-2]/err[-1]) / log(sides[-1]/sides[-2]).
    """
    errs = []
    for s in sides:
        mesh = perturbed_mesh(s, perturb_frac, seed=seed, L=L)
        errs.append(_solution_error(mesh, variant))
    # Finest-two slope
    slope = math.log(errs[-2] / errs[-1]) / math.log(sides[-1] / sides[-2])
    return slope


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------

def vr_sweep(perturb_fracs=(0.0, 0.1, 0.25, 0.5)):
    """Run the distortion characterization sweep.

    For each perturbation fraction, compute solution-error slopes for both
    operators and return list of dicts:
        {perturb, baseline_slope, cotangent_slope}

    Scientific question recorded: does the cotangent operator degrade LESS
    than baseline under distortion (a robustness win), or do both degrade
    similarly (distortion sensitivity is not a Hodge problem)?
    """
    rows = []
    for pf in perturb_fracs:
        baseline_slope = vr_solution_error(pf, "baseline")
        cotangent_slope = vr_solution_error(pf, "cotangent")
        rows.append({
            "perturb": float(pf),
            "baseline_slope": baseline_slope,
            "cotangent_slope": cotangent_slope,
        })
    return rows
