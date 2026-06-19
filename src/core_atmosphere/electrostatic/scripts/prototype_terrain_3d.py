#!/usr/bin/env python3
"""3-D terrain prototype: triangular-lattice horizontal x terrain-following layers."""
from __future__ import annotations
import math
from collections import deque
from types import SimpleNamespace
import importlib.util as _ilu
import pathlib as _pl
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import Delaunay, cKDTree

_MFD = _pl.Path(__file__).resolve().parent / "mfd_operator.py"
def _mfd():
    s = _ilu.spec_from_file_location("mfd_operator", _MFD)
    mod = _ilu.module_from_spec(s); s.loader.exec_module(mod); return mod

def _tri_lattice(n_side, L):
    """Regular triangular lattice of points on [0,L]^2 (rows offset by half-step)."""
    dx = L / n_side
    pts = []
    for j in range(n_side):
        y = (j + 0.5) * dx * math.sqrt(3.0) / 2.0
        xoff = 0.0 if j % 2 == 0 else 0.5 * dx
        for i in range(n_side):
            pts.append([(i + 0.5) * dx + xoff, y])
    pts = np.array(pts)
    pts[:, 1] *= L / pts[:, 1].max()          # rescale y to fill [0,L]
    return pts                                  # (nCellsH, 2)

def terrain_mesh_3d(n_side, nz, hill_fraction, L=1.0, H=1.0):
    xy = _tri_lattice(n_side, L)               # (nCellsH, 2)
    nCellsH = xy.shape[0]
    # smooth interior bump vanishing on the boundary: sin(pi x/L) sin(pi y/L)
    bump = np.sin(math.pi * xy[:, 0] / L) * np.sin(math.pi * xy[:, 1] / L)
    z_s = hill_fraction * H * bump             # (nCellsH,)
    zeta = (np.arange(nz) + 0.5) * (H / nz)    # (nz,)
    # terrain-following: z = z_s*(1 - zeta/H) + zeta
    zcol = z_s[:, None] * (1.0 - zeta[None, :] / H) + zeta[None, :]   # (nCellsH, nz)
    xyz = np.zeros((nCellsH * nz, 3))
    for ih in range(nCellsH):
        for k in range(nz):
            xyz[ih * nz + k] = (xy[ih, 0], xy[ih, 1], zcol[ih, k])
    # interior mask: drop boundary columns (touching the domain edge) and k=0,nz-1
    edge = (np.abs(xy[:, 0]) < 1e-9) | (np.abs(xy[:, 0] - L) < 1e-9) | \
           (np.abs(xy[:, 1]) < 1e-9) | (np.abs(xy[:, 1] - L) < 1e-9)
    # mark a column as boundary if it is within ~1 spacing of the edge
    dx = L / n_side
    near_edge = (xy[:, 0] < dx) | (xy[:, 0] > L - dx) | (xy[:, 1] < dx) | (xy[:, 1] > L - dx)
    interior = np.zeros(nCellsH * nz, dtype=bool)
    for ih in range(nCellsH):
        for k in range(nz):
            interior[ih * nz + k] = (not near_edge[ih]) and (0 < k < nz - 1)
    def cell_index(ih, k): return ih * nz + k
    return SimpleNamespace(xyz=xyz, nCellsH=nCellsH, nz=nz, xy=xy, zcol=zcol,
                           interior=interior, cell_index=cell_index, L=L, H=H,
                           z_s=z_s)


# ---------------------------------------------------------------------------
# B2: 3-D terrain cotangent operator (Delaunay tetrahedralization + P1 stiffness)
# ---------------------------------------------------------------------------

def terrain_operator_ungrounded(mesh):
    """Assemble A = Σ_tet K_tet (P1 stiffness) over Delaunay tets; symmetric PSD,
    annihilates constants (K@1=0 per tet, so A@1=0).

    Two defenses against structured-lattice Delaunay pathologies:
    1. Volume filter: skip tets with vol < 1e-12 (degenerate slivers that crash inv).
    2. Horizontal-reach filter: skip tets where any horizontal column pair is more
       than 2*dx apart. This removes the convex-hull boundary tets that Delaunay on
       a structured point cloud always creates — these bridge corner points many grid
       spacings apart and are physically spurious (they do not represent neighboring
       terrain columns). Threshold 2*dx retains genuine near-diagonal cross-column
       couplings (~sqrt(2)*dx diagonal) while rejecting long-range artifacts.
       This is a deterministic, geometry-based filter (no jitter, no randomness).
    """
    mfd = _mfd()
    N = mesh.xyz.shape[0]
    nz = mesh.nz
    dx = mesh.L / int(round(math.sqrt(mesh.nCellsH)))
    h_thresh = 2.0 * dx          # max allowed horizontal span within one tet
    tri = Delaunay(mesh.xyz)
    rows, cols, data = [], [], []
    for tet in tri.simplices:
        # Filter degenerate tets BEFORE calling tet_p1_stiffness (which inverts M).
        p = mesh.xyz[tet]
        M = np.column_stack((np.ones(4), p))
        vol = abs(np.linalg.det(M)) / 6.0
        if vol < 1e-12:
            continue
        # Filter long-range convex-hull artifacts: skip tets spanning distant columns.
        ihs = [int(i) // nz for i in tet]
        max_h = max(np.linalg.norm(mesh.xy[ihs[a]] - mesh.xy[ihs[b]])
                    for a in range(4) for b in range(a + 1, 4))
        if max_h > h_thresh:
            continue
        K = mfd.tet_p1_stiffness(p)
        for a in range(4):
            for b in range(4):
                rows.append(int(tet[a])); cols.append(int(tet[b])); data.append(K[a, b])
    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N))
    return 0.5 * (A + A.T)


def terrain_operator(mesh, ground_cell=0):
    """Ground the ungrounded operator: zero row/col for ground_cell, preserve diagonal."""
    A = terrain_operator_ungrounded(mesh).tolil()
    d = A[ground_cell, ground_cell]
    A[ground_cell, :] = 0.0
    A[:, ground_cell] = 0.0
    A[ground_cell, ground_cell] = d
    return A.tocsr()


def spd_min_eig_of(A):
    """Smallest algebraic eigenvalue via ARPACK (SA shift-invert)."""
    return float(spla.eigsh(A, k=1, which="SA", return_eigenvectors=False)[0])


def _horizontal_adjacency(mesh):
    """Nearest-neighbor column adjacency from the triangular lattice (~6 neighbors)."""
    dx = mesh.L / int(round(math.sqrt(mesh.nCellsH)))
    tree = cKDTree(mesh.xy)
    adj = [[] for _ in range(mesh.nCellsH)]
    for ih in range(mesh.nCellsH):
        idx = tree.query_ball_point(mesh.xy[ih], r=1.5 * dx)
        adj[ih] = [j for j in idx if j != ih]
    return adj


def halo_rings(mesh):
    """Max horizontal cell-ring coupled by the operator (Fortran halo-width feasibility)."""
    A = terrain_operator_ungrounded(mesh).tocoo()
    adj = _horizontal_adjacency(mesh)

    def hdist(src, dst, cap=6):
        if src == dst:
            return 0
        seen = {src}; q = deque([(src, 0)])
        while q:
            c, dd = q.popleft()
            if dd >= cap:
                continue
            for nb in adj[c]:
                if nb == dst:
                    return dd + 1
                if nb not in seen:
                    seen.add(nb); q.append((nb, dd + 1))
        return cap

    worst = 0
    nz = mesh.nz
    for r, c, v in zip(A.row, A.col, A.data):
        if r == c or abs(v) < 1e-14:
            continue
        worst = max(worst, hdist(r // nz, c // nz))
    return worst
