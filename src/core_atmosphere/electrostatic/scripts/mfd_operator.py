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


def spd_min_eig(A_csr):
    """Smallest algebraic eigenvalue (SPD check) of a grounded operator."""
    return float(spla.eigsh(A_csr, k=1, which="SA", return_eigenvectors=False)[0])
