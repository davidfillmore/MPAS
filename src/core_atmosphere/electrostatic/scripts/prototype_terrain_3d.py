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
from scipy.spatial import Delaunay

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
# B2: 3-D terrain cotangent operator (STRUCTURED conforming prism->tet mesh)
# ---------------------------------------------------------------------------

def _structured_prism_tets(mesh):
    """Return ALL tets of a STRUCTURED, CONFORMING prism->tet mesh of the column
    extrusion — no Delaunay-in-3-D, no volume filter, no horizontal-reach filter.

    Construction
    ------------
    1. Horizontal triangulation: ``scipy.spatial.Delaunay(mesh.xy)`` of the column
       centres (a 2-D Delaunay of the lattice columns is well behaved) gives the
       column-triangle connectivity only.
    2. For every horizontal triangle (columns a,b,c) and every vertical gap between
       consecutive cell-centre levels L and L+1 (L = 0 .. nz-2), the six nodes
       (a,L),(b,L),(c,L),(a,L+1),(b,L+1),(c,L+1) form a triangular PRISM.  Each
       prism is split into 3 tetrahedra.

    Conformity (CRITICAL)
    ---------------------
    The quadrilateral side faces shared between horizontally-adjacent prisms must
    be split by the SAME diagonal from both sides, or the mesh is non-conforming
    (gaps/overlaps) and the operator is inconsistent.  We use the canonical
    sorted-vertex Freudenthal-Kuhn prism split:

        * Order the three columns by GLOBAL column index -> (g0 < g1 < g2) with
          local bottom labels (s0,s1,s2) and tops (t0,t1,t2) = (s0+? ...) directly
          above.
        * Apply the fixed 3-tet template to the SORTED prism:
              (s0,s1,s2,t2), (s0,s1,t2,t1), (s0,t1,t2,t0).

    Because the global column index is a TOTAL order, any two prisms that share a
    vertical quad face (i.e. share two columns) sort those two columns identically,
    so they pick the IDENTICAL face diagonal.  The split therefore tiles the column
    extrusion with no gaps or overlaps (verified: the 3 tet volumes sum exactly to
    the prism volume for every column ordering, and adjacent prisms agree on the
    shared-face diagonal).

    Node indexing matches ``mesh.cell_index(ih, k) = ih*nz + k``.

    Returns list of ``(tet, vol)`` where ``tet`` is a length-4 int array of node
    indices and ``vol`` is the physical tet volume.  ``terrain_operator_ungrounded``,
    ``_cell_volumes`` and ``halo_rings`` all consume this single helper so the
    operator (A @ phi) and the RHS (source * V_cell) integrate over identical,
    fully-covering supports.
    """
    nz = mesh.nz
    xyz = mesh.xyz

    def node(ih, k):
        return ih * nz + k

    def tet_vol(t):
        p = xyz[list(t)]
        M = np.column_stack((np.ones(4), p))
        return abs(np.linalg.det(M)) / 6.0

    tri = Delaunay(mesh.xy)                      # 2-D Delaunay: column connectivity
    tets = []
    for simplex in tri.simplices:
        cols = sorted(int(c) for c in simplex)   # ascending GLOBAL column index
        c0, c1, c2 = cols
        for kL in range(nz - 1):
            kU = kL + 1
            # local bottom verts (s0,s1,s2) sorted by global column, tops above
            s0, s1, s2 = node(c0, kL), node(c1, kL), node(c2, kL)
            t0, t1, t2 = node(c0, kU), node(c1, kU), node(c2, kU)
            for t in ((s0, s1, s2, t2), (s0, s1, t2, t1), (s0, t1, t2, t0)):
                v = tet_vol(t)
                if v <= 0.0:
                    continue                     # only skip exactly-degenerate tets
                tets.append((np.array(t, dtype=int), v))
    return tets


def terrain_operator_ungrounded(mesh):
    """Assemble A = Sum_tet K_tet (P1 stiffness) over the STRUCTURED conforming
    prism->tet mesh; symmetric PSD, annihilates constants (K@1 = 0 per tet, so
    A@1 = 0).

    The tet set comes from ``_structured_prism_tets`` (full coverage of the column
    extrusion, no filters).  P1 stiffness on a conforming full-coverage tet mesh is
    a consistent discrete Laplacian, and SPD + constant-annihilation are automatic.
    """
    mfd = _mfd()
    N = mesh.xyz.shape[0]
    rows, cols, data = [], [], []
    for tet, _vol in _structured_prism_tets(mesh):
        p = mesh.xyz[tet]
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
    """Nearest-neighbor column adjacency = the EDGES of the same 2-D Delaunay
    triangulation the operator's prisms are built from.

    Deriving adjacency from the triangulation (rather than a fixed KD-tree radius)
    makes the halo measurement consistent with the actual operator support: two
    columns are 1-ring neighbours iff they share a horizontal-triangle edge, which
    is exactly when their prism columns can couple in A.  A fixed 1.5*dx radius
    misses some genuine edge neighbours on the staggered lattice's domain edge
    (offset rows sit ~2.1*dx apart) and would spuriously report those true
    neighbours as a 2nd ring.
    """
    tri = Delaunay(mesh.xy)
    adj = [set() for _ in range(mesh.nCellsH)]
    for simplex in tri.simplices:
        a, b, c = (int(i) for i in simplex)
        for u, v in ((a, b), (b, c), (a, c)):
            adj[u].add(v); adj[v].add(u)
    return [sorted(s) for s in adj]


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


# ---------------------------------------------------------------------------
# B3: MMS consistency gate (Axis 1) -- honest pointwise truncation
# ---------------------------------------------------------------------------

def _cell_volumes(mesh):
    """Per-cell lumped-mass volume: 1/4 of the sum of tet volumes incident to the cell.

    Uses ``_structured_prism_tets`` (the SAME conforming, full-coverage tet set the
    operator assembles over) rather than a fresh Delaunay, so (A @ phi) and
    b = eps*lambda*phi*V integrate over identical supports — required for a
    self-consistent P1 MMS residual.
    """
    vol = np.zeros(mesh.xyz.shape[0])
    for tet, v in _structured_prism_tets(mesh):
        for a in tet:
            vol[a] += v / 4.0
    return vol


def coverage_fraction(mesh):
    """Fraction of the column-extrusion volume covered by the structured tets.

    Denominator is the ANALYTIC column-extrusion volume: the 2-D triangulated
    column-hull, vertically extruded between the bottom (k=0) and top (k=nz-1)
    cell-centre levels (terrain-following column thicknesses included).  A
    conforming, gap-free tiling returns ~1.0.  (The structured tets exactly tile
    this region by construction.)
    """
    nz = mesh.nz
    tri = Delaunay(mesh.xy)
    # analytic extrusion volume = sum over horizontal triangles of (base area *
    # mean column extent from level 0 to level nz-1).
    analytic = 0.0
    for simplex in tri.simplices:
        a, b, c = (int(i) for i in simplex)
        p = mesh.xy[[a, b, c]]
        base = 0.5 * abs((p[1, 0] - p[0, 0]) * (p[2, 1] - p[0, 1])
                         - (p[2, 0] - p[0, 0]) * (p[1, 1] - p[0, 1]))
        extent = np.mean([mesh.zcol[a, nz - 1] - mesh.zcol[a, 0],
                          mesh.zcol[b, nz - 1] - mesh.zcol[b, 0],
                          mesh.zcol[c, nz - 1] - mesh.zcol[c, 0]])
        analytic += base * extent
    tet_total = sum(v for _t, v in _structured_prism_tets(mesh))
    return float(tet_total / analytic)


def _terrain_residual_details(mesh, eps):
    """Pointwise + L2 MMS truncation diagnostics for the structured operator.

    Manufactured solution phi = sin(kx x) sin(ky y) sin(kz z) (zero on the box
    boundary).  Consistent lumped-mass FV RHS b_i = eps*(kx^2+ky^2+kz^2)*phi_i*V_i
    with V_i the lumped-mass cell volume from the structured tet set.  Residual
    r = A_ung @ phi - b.

    Returns a dict over INTERIOR cells (mesh.interior) restricted further to cells
    with |b_i| > 1e-3 * max|b| (drop near-nodal cells where b ~ 0 and the relative
    measure is ill-posed):
        pointwise_max : max_i |r_i / b_i|   (the honest, normalization-independent
                        per-cell relative truncation; ~decreasing-with-h for a
                        consistent operator)
        pointwise_rms : RMS of r_i / b_i
        rel_l2        : ||r||_2 / ||b||_2   (the convention the 2-D prototype's
                        _mms_residual_relative_l2 uses; cross-check)
    No ||phi|| normalization anywhere.
    """
    A = terrain_operator_ungrounded(mesh)
    x, y, z = mesh.xyz[:, 0], mesh.xyz[:, 1], mesh.xyz[:, 2]
    kx = ky = math.pi / mesh.L
    kz = math.pi / mesh.H
    phi = np.sin(kx * x) * np.sin(ky * y) * np.sin(kz * z)
    V = _cell_volumes(mesh)
    b = eps * (kx**2 + ky**2 + kz**2) * phi * V
    r = A @ phi - b

    mask = mesh.interior & (np.abs(b) > 1e-3 * np.max(np.abs(b)))
    rb = r[mask] / b[mask]
    return {
        "pointwise_max": float(np.max(np.abs(rb))),
        "pointwise_rms": float(np.sqrt(np.mean(rb**2))),
        "rel_l2": float(np.linalg.norm(r[mask]) / np.linalg.norm(b[mask])),
        "n_active": int(mask.sum()),
    }


def terrain_mms_order(hill_fraction, sides=(8, 12, 18), nz_of=None,
                      L=1.0, H=1.0, eps=1.0, return_details=False):
    """Estimate MMS convergence order for the 3-D terrain cotangent operator using
    the HONEST pointwise relative truncation ``max_i |r_i / b_i|``.

    For each n_side in ``sides`` builds the terrain mesh, samples the manufactured
    solution phi = sin(pi x/L) sin(pi y/L) sin(pi z/H), forms the consistent FV RHS
    b = eps*((pi/L)^2+(pi/L)^2+(pi/H)^2)*phi*V_cell with V_cell the lumped-mass
    volume from the structured tet set, and measures the pointwise relative
    truncation over interior cells with |b_i| > 1e-3*max|b|.  Returns the finest-two
    log-log slope of ``pointwise_max`` as the convergence-order estimate.  A
    consistent 2nd-order operator gives slope ~2 (error decreasing with refinement).

    With return_details=True also returns (sides, list-of-detail-dicts).
    """
    if nz_of is None:
        nz_of = lambda s: s
    details = []
    for s in sides:
        g = terrain_mesh_3d(s, nz_of(s), hill_fraction, L, H)
        details.append(_terrain_residual_details(g, eps))
    em = [d["pointwise_max"] for d in details]
    slope = math.log(em[-2] / em[-1]) / math.log(sides[-1] / sides[-2])
    if return_details:
        return slope, list(sides), details
    return slope
