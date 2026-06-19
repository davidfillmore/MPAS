#!/usr/bin/env python3
"""Mimetic (MFD) core for the globally-consistent Poisson operator.

A = d0ᵀ H d0 with H a per-cell CONSISTENT SPD block (the mimetic Hodge /
transmissibility) replacing the diagonal TPFA Hodge. Consistency is structural:
    T_E @ C_E == eps * N_E,   using the geometric identity  N_Eᵀ @ C_E = vol·I,
where N_E are area-weighted outward face normals and C_E = face_centroid - centroid.
The stability term is SPD on the complement of range(C_E) and vanishes on range(C_E),
so it never disturbs consistency. A = d0ᵀ H d0 is SPD whenever every block is SPD.
"""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def mimetic_hodge_block(face_normals, face_centroids, cell_centroid,
                        cell_volume, eps, gamma=1.0):
    """Consistent SPD mimetic Hodge block for one cell (couples its faces)."""
    N = np.asarray(face_normals, dtype=float)                 # (nf, d) area*normal
    C = np.asarray(face_centroids, dtype=float) - np.asarray(cell_centroid, float)
    nf, d = N.shape
    K = float(eps) * np.eye(d)
    # Consistent part: reproduces eps*grad of linear fields exactly (rank d).
    T0 = (N @ K @ N.T) / float(cell_volume)
    # Stability: SPD on range(C)^perp, zero on range(C) (keeps consistency).
    Pc = C @ np.linalg.solve(C.T @ C, C.T)                    # projector onto range(C)
    scale = gamma * np.trace(T0) / nf
    T = T0 + scale * (np.eye(nf) - Pc)
    return 0.5 * (T + T.T)


def assemble_dT_H_d(n_cells, faces, blocks, face_cells, ground_cell=0):
    """Assemble grounded A = d0ᵀ H d0.

    n_cells   : int
    faces     : list[list[int]]  faces[c] = ordered global face ids of cell c
                (same order as the rows/cols of blocks[c]).
    blocks    : list[np.ndarray] blocks[c] = T_E for cell c (len(faces[c]) square).
    face_cells: (n_faces, 2) int  cellsOnFace; interior face has both >= 0, a
                boundary face has the absent side = -1. Sign of d0 on face f is
                +1 for face_cells[f,1], -1 for face_cells[f,0].
    ground_cell: int cell to ground (Dirichlet gauge), or None to return the raw
                un-grounded d0ᵀ H d0 (the honest object for interior consistency
                measures, which a boundary condition would contaminate).
    Returns CSR (n_cells x n_cells), grounded unless ground_cell is None.
    """
    n_faces = face_cells.shape[0]
    # Signed incidence d0 (faces x cells): (d0 phi)_f = phi[hi] - phi[lo].
    dr, dc, dd = [], [], []
    for f in range(n_faces):
        lo, hi = int(face_cells[f, 0]), int(face_cells[f, 1])
        if lo >= 0:
            dr.append(f); dc.append(lo); dd.append(-1.0)
        if hi >= 0:
            dr.append(f); dc.append(hi); dd.append(1.0)
    d0 = sp.csr_matrix((dd, (dr, dc)), shape=(n_faces, n_cells))
    # Global face-Hodge H: sum each cell's block into its faces' (row,col) slots.
    hr, hc, hd = [], [], []
    for c in range(n_cells):
        fs = faces[c]; T = blocks[c]
        for a in range(len(fs)):
            for b in range(len(fs)):
                v = T[a, b]
                if v != 0.0:
                    hr.append(fs[a]); hc.append(fs[b]); hd.append(v)
    H = sp.csr_matrix((hd, (hr, hc)), shape=(n_faces, n_faces))
    A = (d0.T @ H @ d0).tolil()
    A = 0.5 * (A + A.T)
    if ground_cell is None:
        return A.tocsr()
    diag_g = A[ground_cell, ground_cell]
    A[ground_cell, :] = 0.0; A[:, ground_cell] = 0.0
    A[ground_cell, ground_cell] = diag_g
    return A.tocsr()


def cotangent_hodge_weights(triangle_xyz):
    """Per-triangle contribution to the DIAGONAL scalar edge-Hodge ``⋆₁``.

    This is the lowest-order discrete-exterior-calculus (DEC) / lowest-order
    Whitney scalar Hodge that pairs with a PLAIN ±1 incidence ``d0`` to give the
    consistent scalar Poisson operator ``A = d0ᵀ ⋆₁ d0`` (the cotangent / P1
    nodal stiffness).  For an interior primal edge ``e`` shared by two triangles
    the accumulated weight equals ``½(cot α + cot β) = l_e / d_e`` (dual-edge
    length over primal-edge length), exactly the Voronoi/Delaunay scalar Hodge.

    Returns ``(weights, local_pairs)`` with ``local_pairs = [(0,1),(1,2),(0,2)]``
    and ``weights[m]`` the cotangent half-weight of local edge ``local_pairs[m]``
    (the cotangent of the angle at the opposite vertex, times ½).  Sum these over
    all triangles into a per-edge diagonal ``⋆₁``.

    NOTE — why NOT the full Whitney 1-form edge-mass.  The lowest-order Whitney
    *1-form* mass ``M₁ = ∫ w_e·w_e'`` (``prototype_tier_A2``'s
    ``whitney_triangle_mass``) is a DIFFERENT object: its DIAGONAL is
    ``5/6·l_e/d_e`` (proportional to ⋆₁) but its OFF-DIAGONALS make
    ``d0ᵀ M₁ d0`` NOT the scalar Laplacian — assembling it everywhere on a
    structured terrain mesh fails the flat 2nd-order guard (slope ≈ −0.02). The
    scalar Poisson operator needs the DIAGONAL Hodge built here. (On the sphere
    side ``M₁`` is used only as a localized few-edge defect enrichment, with the
    diagonal ``l_e/d_e`` retained on every bulk edge — consistent with this.)
    """
    p = np.asarray(triangle_xyz, dtype=float)
    if p.shape[1] == 3:
        # Work in the triangle's own plane (z≈0 embedding is fine, but be general).
        e0 = p[1] - p[0]
        e1 = p[2] - p[0]
        u = e0 / np.linalg.norm(e0)
        w = e1 - (e1 @ u) * u
        v = w / np.linalg.norm(w)
        p = np.column_stack(((p - p[0]) @ u, (p - p[0]) @ v))
    local_pairs = [(0, 1), (1, 2), (0, 2)]
    # cot of the angle at vertex i is the weight of the OPPOSITE edge (j,k).
    weight_by_edge = {}
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        a = p[j] - p[i]
        b = p[k] - p[i]
        cross = abs(a[0] * b[1] - a[1] * b[0])
        cot = (a @ b) / cross if cross > 1e-300 else 0.0
        weight_by_edge[frozenset((j, k))] = 0.5 * cot
    weights = np.array([weight_by_edge[frozenset(pair)] for pair in local_pairs])
    return weights, local_pairs


def tet_p1_stiffness(tet_xyz):
    """P1 (linear-FEM) local stiffness for one tetrahedron — the 3-D cotangent
    Hodge contribution. tet_xyz: (4,3). K_ij = vol * (grad lambda_i . grad lambda_j),
    symmetric PSD, K @ 1 = 0. Reduces, when assembled, to the cotangent Laplacian."""
    p = np.asarray(tet_xyz, dtype=float)
    M = np.column_stack((np.ones(4), p))        # rows: [1, x, y, z] per vertex
    Minv = np.linalg.inv(M)
    grads = Minv[1:4, :].T                        # (4,3): grad of each barycentric basis
    vol = abs(np.linalg.det(M)) / 6.0
    K = vol * (grads @ grads.T)
    return 0.5 * (K + K.T)


def spd_min_eig(A_csr):
    """Smallest algebraic eigenvalue (SPD check) of a grounded operator."""
    return float(spla.eigsh(A_csr, k=1, which="SA", return_eigenvectors=False)[0])
