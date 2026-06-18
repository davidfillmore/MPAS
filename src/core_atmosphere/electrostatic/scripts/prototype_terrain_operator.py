#!/usr/bin/env python3
"""
2-D (x–z) terrain-following coordinate prototype — Task A1.

Provides `terrain_grid`, the analytic grid and metric foundation that
Tasks A2 (slope-corrected operator) and A3 (MMS order gate) build on.

Coordinate system
-----------------
Physical domain: x in [0, L) (periodic), z in [0, H].
Surface hill:

    z_s(x)     = hill_height * cos(2*pi*x/L)

Terrain-following coordinate zeta in [0, H]:

    z(x, zeta) = zeta + z_s(x) * (1 - zeta/H)

so z=z_s(x) at zeta=0 (surface) and z=H at zeta=H (top, flat).

Analytic metrics
----------------
    dz/dzeta              = 1 - z_s(x)/H
    zz  = dzeta/dz        = 1 / (1 - z_s(x)/H)          [vertical metric]
    z_s'(x)               = -hill_height*(2*pi/L)*sin(2*pi*x/L)
    dzdx(x, zeta)         = z_s'(x) * (1 - zeta/H)       [terrain slope metric]

Note: dzdx → 0 at the top (zeta = H) and everywhere when hill_height = 0.

Grid layout
-----------
Cell centres:
    x[i]       = (i + 0.5) * L/nx,      i = 0 .. nx-1
    zeta[k]    = (k + 0.5) * H/nz,      k = 0 .. nz-1

Interfaces:
    x_int[i]   = i * L/nx,              i = 0 .. nx      (not stored; not needed)
    zeta_int[k]= k * H/nz,              k = 0 .. nz
"""

from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np
import scipy.sparse as sp


def terrain_grid(nx: int, nz: int, hill_height: float, L: float, H: float):
    """Build the analytic 2-D terrain-following coordinate grid and metrics.

    Parameters
    ----------
    nx, nz : int
        Number of cells in the x and zeta (vertical) directions.
    hill_height : float
        Amplitude of the cosine surface hill.  Must satisfy
        |hill_height| < H to keep the Jacobian positive everywhere.
    L, H : float
        Domain length in x and total height.

    Returns
    -------
    SimpleNamespace with attributes
    --------------------------------
    dx     : float           uniform x spacing  = L/nx
    dzeta  : float           uniform zeta spacing = H/nz
    x      : (nx,)           cell-centre x positions (periodic)
    zeta   : (nz,)           cell-centre zeta values
    zeta_int : (nz+1,)       interface zeta values (0, dzeta, …, H)
    z      : (nx, nz)        physical height at cell centres
    z_int  : (nx, nz+1)      physical height at zeta interfaces
    dz     : (nx, nz)        layer thickness in physical z (z_int diff)
    dzdx   : (nx, nz)        terrain slope metric at cell centres
    dzdx_int : (nx, nz+1)    terrain slope metric at zeta interfaces
    zz     : (nx, nz)        vertical metric dzeta/dz at cell centres
    """
    dx = L / nx
    dzeta = H / nz

    # --- 1-D coordinate arrays -----------------------------------------------
    x = (np.arange(nx) + 0.5) * dx                  # (nx,)
    zeta = (np.arange(nz) + 0.5) * dzeta            # (nz,)
    zeta_int = np.arange(nz + 1) * dzeta            # (nz+1,)

    # --- Surface hill and its x-derivative ------------------------------------
    # z_s(x) = hill_height * cos(2*pi*x/L)
    z_s = hill_height * np.cos(2.0 * math.pi * x / L)           # (nx,)
    # z_s'(x) = -hill_height * (2*pi/L) * sin(2*pi*x/L)
    dz_s = -hill_height * (2.0 * math.pi / L) * np.sin(2.0 * math.pi * x / L)  # (nx,)

    # --- Physical heights -----------------------------------------------------
    # z(x, zeta) = zeta + z_s(x) * (1 - zeta/H)
    # Broadcast: z_s is (nx,), zeta is (nz,) → result (nx, nz)
    # (1 - zeta/H) factor at cell centres
    fac_c = 1.0 - zeta / H          # (nz,)
    fac_i = 1.0 - zeta_int / H      # (nz+1,)  — exactly 0 at k=nz

    z = zeta[np.newaxis, :] + z_s[:, np.newaxis] * fac_c[np.newaxis, :]       # (nx, nz)
    z_int = zeta_int[np.newaxis, :] + z_s[:, np.newaxis] * fac_i[np.newaxis, :]  # (nx, nz+1)

    # Layer thickness (difference of interface z values)
    dz = np.diff(z_int, axis=1)    # (nx, nz)

    # --- Metric fields --------------------------------------------------------
    # dzdx = z_s'(x) * (1 - zeta/H)   [terrain slope metric]
    dzdx = dz_s[:, np.newaxis] * fac_c[np.newaxis, :]          # (nx, nz)
    dzdx_int = dz_s[:, np.newaxis] * fac_i[np.newaxis, :]      # (nx, nz+1)

    # zz = dzeta/dz = 1 / (1 - z_s(x)/H)   [vertical metric, column-constant]
    zz_col = 1.0 / (1.0 - z_s / H)          # (nx,)
    zz = np.broadcast_to(zz_col[:, np.newaxis], (nx, nz)).copy()  # (nx, nz)

    return SimpleNamespace(
        dx=dx,
        dzeta=dzeta,
        x=x,
        zeta=zeta,
        zeta_int=zeta_int,
        z=z,
        z_int=z_int,
        dz=dz,
        dzdx=dzdx,
        dzdx_int=dzdx_int,
        zz=zz,
    )


# =============================================================================
# Task A2 — slope-corrected gradient G, weight W, operator A = Gᵀ W G
# =============================================================================
#
# Stored potential phi[i, k] lives on the (x, zeta) grid (i = column, k = level)
# with the flat cell index c(i, k) = i*nz + k.  The PHYSICAL gradient on the
# terrain-following grid follows from the chain rule through z(x, zeta):
#
#     grad_x = d(phi)/dx|_z = d_x phi  -  dzdx * zz * d_zeta phi   (horizontal)
#     grad_z = d(phi)/dz    = zz * d_zeta phi                      (vertical)
#
# G samples grad_x on x-faces (between columns i and i+1, periodic) and grad_z
# on interior z-faces (between levels k and k+1).  The slope correction in
# grad_x is the cross-derivative term: it couples an x-face to cells at level k
# AND k±1 in both adjacent columns, with dzdx and zz interpolated to the face.
#
# The operator is built LITERALLY as A = Gᵀ W G with W diagonal and strictly
# positive, so A is symmetric positive-semidefinite by construction; grounding
# one cell (Dirichlet) removes the constant null mode and makes A SPD.  When the
# slope is zero (hill_height = 0) the cross-derivative term vanishes (dzdx = 0)
# and A reduces exactly to the standard finite-volume 5-point Laplacian.


def _dzeta_stencil(k, nz, dzeta):
    """Stencil entries [(level, coeff), ...] for d_zeta phi at cell-centre level k.

    Centred in the interior; one-sided at the top/bottom boundary where the
    centred neighbour is missing.  Coefficients already include 1/dzeta.
    """
    if nz < 2:
        return []
    if k == 0:
        return [(1, 1.0 / dzeta), (0, -1.0 / dzeta)]
    if k == nz - 1:
        return [(nz - 1, 1.0 / dzeta), (nz - 2, -1.0 / dzeta)]
    return [(k + 1, 0.5 / dzeta), (k - 1, -0.5 / dzeta)]


def terrain_gradient(grid):
    """Slope-corrected physical-gradient operator G on the terrain grid.

    Returns a CSR sparse matrix of shape (n_faces, n_cells) whose rows are the
    x-faces (first n_x = nx*nz rows) followed by the interior z-faces
    (next n_z = nx*(nz-1) rows).  Columns are cells c(i, k) = i*nz + k.

    x-face row fx(i, k) = i*nz + k       (face between columns i and i+1, periodic)
    z-face row fz(i, k) = nx*nz + i*(nz-1) + k   (interface between levels k, k+1)
    """
    nx = grid.x.size
    nz = grid.zeta.size
    dx = grid.dx
    dzeta = grid.dzeta
    dzdx = grid.dzdx
    zz = grid.zz

    n_cells = nx * nz
    n_xf = nx * nz
    n_zf = nx * (nz - 1)

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    def add(r, c, v):
        rows.append(r)
        cols.append(c)
        data.append(v)

    def cell(i, k):
        return i * nz + k

    # --- x-faces: grad_x = d_x phi - dzdx*zz*d_zeta phi ------------------------
    for i in range(nx):
        ip = (i + 1) % nx
        for k in range(nz):
            r = i * nz + k
            # Horizontal centred difference across the face.
            add(r, cell(ip, k), 1.0 / dx)
            add(r, cell(i, k), -1.0 / dx)
            # Slope (cross-derivative) correction, interpolated to the face.
            dzdx_f = 0.5 * (dzdx[i, k] + dzdx[ip, k])
            zz_f = 0.5 * (zz[i, k] + zz[ip, k])
            s = -dzdx_f * zz_f
            if s != 0.0:
                # Average the centred d_zeta phi of the two adjacent columns.
                for j in (i, ip):
                    for kk, coef in _dzeta_stencil(k, nz, dzeta):
                        add(r, cell(j, kk), 0.5 * s * coef)

    # --- z-faces: grad_z = zz * d_zeta phi ------------------------------------
    base = n_xf
    for i in range(nx):
        for k in range(nz - 1):
            r = base + i * (nz - 1) + k
            zz_f = 0.5 * (zz[i, k] + zz[i, k + 1])
            add(r, cell(i, k + 1), zz_f / dzeta)
            add(r, cell(i, k), -zz_f / dzeta)

    G = sp.coo_matrix(
        (data, (rows, cols)), shape=(n_xf + n_zf, n_cells)
    ).tocsr()
    return G


def terrain_weight(grid, eps):
    """Diagonal SPD face weight W (n_faces, n_faces) for A = Gᵀ W G.

    Each weight is eps times the physical control volume sampling that face's
    gradient: dx * dzeta / zz_face (computational cell volume dx*dzeta times the
    Jacobian dz/dzeta = 1/zz).  Strictly positive because zz > 0 and eps > 0, so
    A = Gᵀ W G is symmetric positive-semidefinite.  Face ordering matches
    terrain_gradient (x-faces then z-faces).
    """
    nx = grid.x.size
    nz = grid.zeta.size
    dx = grid.dx
    dzeta = grid.dzeta
    zz = grid.zz
    vol = eps * dx * dzeta

    n_xf = nx * nz
    n_zf = nx * (nz - 1)
    diag = np.empty(n_xf + n_zf)

    # x-faces
    for i in range(nx):
        ip = (i + 1) % nx
        for k in range(nz):
            zz_f = 0.5 * (zz[i, k] + zz[ip, k])
            diag[i * nz + k] = vol / zz_f

    # z-faces
    base = n_xf
    for i in range(nx):
        for k in range(nz - 1):
            zz_f = 0.5 * (zz[i, k] + zz[i, k + 1])
            diag[base + i * (nz - 1) + k] = vol / zz_f

    return sp.diags(diag, format="csr")


def assemble_terrain_operator(grid, eps, ground_cell=0):
    """Assemble the metric-corrected SPD operator A = Gᵀ W G (CSR).

    Built literally as G.T @ W @ G (hence symmetric, positive-semidefinite),
    then grounded at one cell to remove the constant null mode: that cell's row
    and column are zeroed and its (positive) diagonal is kept, decoupling it and
    leaving a Dirichlet Laplacian on the remaining cells => SPD.
    """
    G = terrain_gradient(grid)
    W = terrain_weight(grid, eps)
    A = (G.T @ W @ G).tolil()

    diag_g = A[ground_cell, ground_cell]
    A[ground_cell, :] = 0.0
    A[:, ground_cell] = 0.0
    A[ground_cell, ground_cell] = diag_g
    return A.tocsr()
