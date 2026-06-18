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


# =============================================================================
# Task A3 — operator-only MMS convergence gate (THE GATE)
# =============================================================================
#
# Measures whether A = Gᵀ W G is CONSISTENT (~2nd order) on the sloped mesh by
# applying the operator to a manufactured field and comparing to the consistent
# finite-volume right-hand side on INTERIOR cells.
#
# Manufactured solution (function of PHYSICAL x, z):
#     phi_exact(x, z) = sin(kx*x) * sin(kz*z),     kx = 2*pi/L, kz = pi/H
# Continuous operator (constant eps):
#     -div(eps grad phi) = eps*(kx^2 + kz^2) * phi_exact
# The discrete operator A integrates -div(eps grad phi) over each physical cell
# volume V_cell = dx * dz_physical (dz_physical = dx*dzeta/zz = grid.dz), so the
# CONSISTENT FV right-hand side is
#     b_c = eps*(kx^2 + kz^2) * phi_sampled_c * V_cell_c .
# This V_cell weighting matches the W = eps * dx*dzeta/zz face weighting A2 uses.
#
# Residual r = A @ phi_sampled - b, measured as the relative L2 norm over
# INTERIOR cells only — the top/bottom boundary rows (k = 0, k = nz-1) carry the
# one-sided d_zeta stencil and are excluded so the bulk truncation order is
# isolated, not the boundary stencil.
#
# Note on grounding: the residual uses the RAW (un-grounded) A = Gᵀ W G, built
# from the same A2 terrain_gradient/terrain_weight blocks.  Grounding is a
# single-cell Dirichlet boundary condition (zeroes one row AND column); zeroing
# the column injects an O(1) defect into that cell's interior neighbour, which
# would contaminate — and at fine resolution dominate — the truncation measure.
# The interior-stencil consistency that this gate measures is a property of
# Gᵀ W G independent of that boundary condition, so the ungrounded operator is
# the honest object to measure.


def _mms_residual_relative_l2(grid, eps, kx, kz):
    """Relative L2 norm of the operator-only MMS residual over interior cells.

    r = (Gᵀ W G) @ phi_sampled - b, with b the consistent FV right-hand side
    eps*(kx^2+kz^2)*phi*V_cell.  Interior excludes the k=0 and k=nz-1 boundary
    rows (one-sided d_zeta stencil).  Cell ordering c(i,k) = i*nz + k.
    """
    nx = grid.x.size
    nz = grid.zeta.size

    G = terrain_gradient(grid)
    W = terrain_weight(grid, eps)
    A = (G.T @ W @ G).tocsr()

    # Manufactured field sampled at physical cell centres (nx, nz).
    phi = np.sin(kx * grid.x)[:, np.newaxis] * np.sin(kz * grid.z)

    # Physical cell volume dx * dz_physical (dz = dx*dzeta/zz column thickness).
    V_cell = grid.dx * grid.dz
    b = eps * (kx ** 2 + kz ** 2) * phi * V_cell

    r = A @ phi.reshape(-1) - b.reshape(-1)

    # Interior mask: drop top/bottom boundary rows (one-sided d_zeta stencil).
    mask = np.zeros((nx, nz), dtype=bool)
    mask[:, 1:nz - 1] = True
    mask = mask.reshape(-1)

    r_norm = np.linalg.norm(r[mask])
    b_norm = np.linalg.norm(b.reshape(-1)[mask])
    return r_norm / b_norm


def terrain_mms_order(
    hill_fraction,
    grids=((32, 24), (64, 48), (128, 96)),
    L=1.0,
    H=1.0,
    eps=1.0,
    return_details=False,
):
    """Finest-two-mesh log-log convergence slope of the operator-only residual.

    Refines through ``grids`` (each doubling nx, nz), measures the interior
    relative-L2 residual of A = Gᵀ W G against the consistent FV right-hand
    side, and fits the slope of the finest two meshes (h ~ 1/nx).

    hill_fraction : surface hill amplitude as a fraction of H (0 => flat).
    return_details : if True, also return (grids, errors) for reporting.
    """
    hill_height = hill_fraction * H
    kx = 2.0 * math.pi / L
    kz = math.pi / H

    errs = []
    for nx, nz in grids:
        g = terrain_grid(nx, nz, hill_height, L, H)
        errs.append(_mms_residual_relative_l2(g, eps, kx, kz))

    # Finest two meshes: h halves between them, so log ratio over log(2).
    h_coarse = 1.0 / grids[-2][0]
    h_fine = 1.0 / grids[-1][0]
    slope = math.log(errs[-2] / errs[-1]) / math.log(h_coarse / h_fine)

    if return_details:
        return slope, list(grids), errs
    return slope


# =============================================================================
# MFD globally-consistent operator on the x-z terrain mesh (Task A1 bake-off)
# =============================================================================
#
# The prior naive slope-corrected GᵀWG (above) is SPD but inconsistent on slopes
# (MMS slope ~ -0.52): its explicit column-averaged cross-derivative term is the
# defect.  The MFD route removes that term entirely.  The discrete gradient d0 is
# the PLAIN coordinate-difference between cells; ALL the metric/slope coupling
# lives in the per-cell mimetic Hodge block built from the cell's REAL physical
# (x, z) corner geometry.  The operator is A = d0ᵀ H d0.
#
# Each cell (i, k) is a physical quad whose corners come from the interface
# heights z_int (x at the interfaces i*dx and (i+1)*dx).  The slope therefore
# enters geometrically — through the sloped-quad face normals/centroids fed to
# mimetic_hodge_block — not through an explicit gradient cross-term.
#
# Global face numbering (matches the face lists handed to assemble_dT_H_d):
#   x-faces: fx(i,k) = i*nz + k,              i=0..nx-1 (periodic), k=0..nz-1
#   z-faces: fz(i,k) = nx*nz + i*(nz+1) + k,  i=0..nx-1, k=0..nz (incl. boundary)
# The MMS phi = sin(kx*x)*sin(kz*z) with kz = pi/H is zero at BOTH z=0 and z=H, so
# the ground (k=0) and top (k=nz) faces are consistently treated as Dirichlet-0
# (single-sided d0 rows, face_cells = -1 on the absent side); periodic x-faces are
# interior.  The gauge for the grounded operator comes from grounding one cell.

import importlib.util as _ilu
import pathlib as _pl

_MFD = _pl.Path(__file__).resolve().parent / "mfd_operator.py"


def _mfd_mod():
    spec = _ilu.spec_from_file_location("mfd_operator", _MFD)
    m = _ilu.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def terrain_cell_geometry(grid):
    """Physical-quad geometry per cell -> (faces, face_cells, geom_per_cell).

    For each cell returns its 4 global face ids, the global cellsOnFace table,
    and the per-cell ``(N, fc, centroid, vol)`` built from the PHYSICAL sloped-quad
    corners — this is where the terrain slope enters, geometrically.

    Global face numbering:
        x-faces: fx(i,k) = i*nz + k,                 i=0..nx-1 (periodic), k=0..nz-1
        z-faces: fz(i,k) = nx*nz + i*(nz+1) + k,     i=0..nx-1, k=0..nz (incl bndry)
    """
    nx, nz, dx = grid.x.size, grid.zeta.size, grid.dx
    z_int = grid.z_int                       # (nx, nz+1) physical interface heights
    n_xf = nx * nz
    n_zf = nx * (nz + 1)
    n_faces = n_xf + n_zf

    def fx(i, k):
        return (i % nx) * nz + k

    def fz(i, k):
        return n_xf + (i % nx) * (nz + 1) + k

    def cell(i, k):
        return (i % nx) * nz + k

    # cellsOnFace (lo, hi); sign convention: +hi, -lo (matches d0 = phi_hi - phi_lo)
    face_cells = np.full((n_faces, 2), -1, dtype=int)
    for i in range(nx):
        for k in range(nz):
            face_cells[fx(i, k)] = (cell(i - 1, k), cell(i, k))   # left of (i,k)
    for i in range(nx):
        for k in range(nz + 1):
            lo = cell(i, k - 1) if k > 0 else -1
            hi = cell(i, k) if k < nz else -1
            face_cells[fz(i, k)] = (lo, hi)

    faces, geom = [], []
    for i in range(nx):
        for k in range(nz):
            # physical corners of cell (i,k): x at interfaces i*dx, (i+1)*dx
            xl, xr = i * dx, (i + 1) * dx
            zbl, ztl = z_int[i, k], z_int[i, k + 1]                 # left column
            zbr, ztr = z_int[(i + 1) % nx, k], z_int[(i + 1) % nx, k + 1]
            corners = np.array([[xl, zbl], [xr, zbr], [xr, ztr], [xl, ztl]])
            centroid = corners.mean(0)
            vol = 0.5 * abs(
                (xr - xl) * (ztl + ztr - zbl - zbr)
            )  # trapezoid area = dx * mean(thickness)
            # face order: left, right, bottom, top  (must match `faces` list below)
            edges = [(3, 0), (1, 2), (0, 1), (2, 3)]   # left,right,bottom,top corners
            N, fc = [], []
            for a, b in edges:
                e = corners[b] - corners[a]
                outw = np.array([e[1], -e[0]])
                mid = 0.5 * (corners[a] + corners[b])
                if np.dot(outw, mid - centroid) < 0:
                    outw = -outw
                N.append(outw)
                fc.append(mid)
            faces.append([fx(i, k), fx(i + 1, k), fz(i, k), fz(i, k + 1)])
            geom.append((np.array(N), np.array(fc), centroid, vol))
    return faces, face_cells, geom


def _terrain_mfd_blocks(grid, eps):
    """Shared builder: (n_cells, faces, face_cells, mimetic blocks) for the MFD path.

    Single source of the per-cell geometry + Hodge blocks so the grounded operator
    and the un-grounded residual measure both feed the same numbers into
    ``mfd_operator.assemble_dT_H_d`` (no duplicated assembly).
    """
    mfd = _mfd_mod()
    n_cells = grid.x.size * grid.zeta.size
    faces, face_cells, geom = terrain_cell_geometry(grid)
    blocks = [mfd.mimetic_hodge_block(N, fc, c, vol, eps) for (N, fc, c, vol) in geom]
    return mfd, n_cells, faces, face_cells, blocks


# --- Whitney / P1 scalar-Hodge path (the bake-off PIVOT, the GO variant) -------
# Variant "mfd" above is mis-formulated: A0's mimetic_hodge_block satisfies
# T·C = eps·N with C the HALF-cell face-centroid offset (the mixed/face-pressure
# convention), so it is the wrong primitive for A = d0ᵀ H d0 (full cell-to-cell
# differences); its dense opposite-face coupling injects a spurious distance-2
# stencil and it diverges on the hill.  The consistent route is the lowest-order
# Whitney/DEC SCALAR edge-Hodge: A = d0ᵀ ⋆₁ d0 with a PLAIN ±1 incidence and ⋆₁
# the DIAGONAL cotangent Hodge (l_e/d_e), assembled per triangle from
# mfd_operator.cotangent_hodge_weights.  This is the P1 nodal stiffness on a
# triangulation of the PHYSICAL cell centres, so the terrain slope enters purely
# through the triangle geometry.  See that helper's docstring for why the full
# Whitney 1-form edge-mass (off-diagonals) is NOT the scalar operator.


def terrain_center_triangulation(grid):
    """Triangulate the physical (x, z) cell-centre grid -> (points, triangles).

    Cell centres are at physical ``(grid.x[i], grid.z[i,k])`` with flat index
    ``c(i,k)=i*nz+k``.  Each primal quad of centres
    ``{(i,k),(i+1,k),(i+1,k+1),(i,k+1)}`` (periodic in i, k=0..nz-2) is split by
    the FIXED diagonal (i,k)-(i+1,k+1) into two triangles.  The fixed diagonal is
    deliberate: alternating it (criss-cross) breaks consistency on the sloped
    mesh (verified: hill slope collapses to ~0).  The periodic seam in x is
    handled by unwrapping the right column's x by +L inside each seam triangle so
    triangle areas/angles are physical.

    Returns ``(points, triangles)``: ``points`` is (n_cells, 2) physical centre
    coordinates (un-unwrapped); ``triangles`` is a list of
    ``(cell_triple, triangle_xyz)`` with ``triangle_xyz`` the (3, 3) seam-unwrapped
    corner coordinates (z=0 third component) ready for the Hodge primitive.
    """
    nx, nz = grid.x.size, grid.zeta.size
    L = nx * grid.dx
    z = grid.z                                   # (nx, nz) physical centre heights

    def cell(i, k):
        return (i % nx) * nz + k

    points = np.empty((nx * nz, 2))
    for i in range(nx):
        for k in range(nz):
            points[cell(i, k)] = (grid.x[i], z[i, k])

    def xyz(c, xshift):
        return np.array([points[c, 0] + xshift, points[c, 1], 0.0])

    triangles = []
    for i in range(nx):
        ip = (i + 1) % nx
        xshift = L if ip < i else 0.0           # unwrap periodic seam
        for k in range(nz - 1):
            c00, c10 = cell(i, k), cell(ip, k)
            c11, c01 = cell(ip, k + 1), cell(i, k + 1)
            p00, p01 = xyz(c00, 0.0), xyz(c01, 0.0)
            p10, p11 = xyz(c10, xshift), xyz(c11, xshift)
            # Fixed diagonal c00-c11.
            triangles.append(((c00, c10, c11), np.array([p00, p10, p11])))
            triangles.append(((c00, c11, c01), np.array([p00, p11, p01])))
    return points, triangles


def _whitney_terrain_operator_ungrounded(grid, eps):
    """Ungrounded A = d0ᵀ ⋆₁ d0 (scalar diagonal Hodge) on the centre triangulation.

    ⋆₁ is the per-edge sum of cotangent half-weights from
    ``mfd_operator.cotangent_hodge_weights`` over the local triangles; ``d0`` is
    the plain ±1 cell-difference incidence ``(d0 φ)_e = φ[b]−φ[a]``.  Single
    source for both the grounded operator and the residual measure.
    """
    mfd = _mfd_mod()
    n_cells = grid.x.size * grid.zeta.size
    _, triangles = terrain_center_triangulation(grid)

    edge_index = {}
    star = {}

    def edge_id(a, b):
        key = (min(a, b), max(a, b))
        e = edge_index.get(key)
        if e is None:
            e = len(edge_index)
            edge_index[key] = e
        return e

    for triple, tri_xyz in triangles:
        weights, local_pairs = mfd.cotangent_hodge_weights(tri_xyz)
        for m, (la, lb) in enumerate(local_pairs):
            e = edge_id(triple[la], triple[lb])
            star[e] = star.get(e, 0.0) + eps * weights[m]

    n_edges = len(edge_index)
    star_diag = np.array([star[e] for e in range(n_edges)])
    dr, dc, dd = [], [], []
    for (a, b), e in edge_index.items():
        dr.extend((e, e)); dc.extend((a, b)); dd.extend((-1.0, 1.0))
    d0 = sp.csr_matrix((dd, (dr, dc)), shape=(n_edges, n_cells))
    H = sp.diags(star_diag)
    A = (d0.T @ H @ d0).tocsr()
    return 0.5 * (A + A.T)


def _terrain_operator_ungrounded(grid, eps, variant):
    """Ungrounded terrain operator for a given variant (shared by operator + measure)."""
    if variant == "mfd":
        mfd, n_cells, faces, face_cells, blocks = _terrain_mfd_blocks(grid, eps)
        return mfd.assemble_dT_H_d(n_cells, faces, blocks, face_cells, ground_cell=None)
    if variant == "whitney":
        return _whitney_terrain_operator_ungrounded(grid, eps)
    raise ValueError(
        f"terrain mesh supports variant 'mfd' or 'whitney', got {variant!r}"
    )


def _ground_operator(A, ground_cell):
    """Return A with one cell grounded (Dirichlet gauge): zero its row/col, keep diag."""
    A = A.tolil()
    diag_g = A[ground_cell, ground_cell]
    A[ground_cell, :] = 0.0
    A[:, ground_cell] = 0.0
    A[ground_cell, ground_cell] = diag_g
    return A.tocsr()


def mfd_terrain_operator(grid, eps, variant="mfd", ground_cell=0):
    """Grounded globally-consistent terrain operator A = d0ᵀ H d0.

    ``variant="mfd"``     -> A0 mimetic-block route (mis-formulated; see notes).
    ``variant="whitney"`` -> lowest-order Whitney/P1 scalar Hodge (the GO route).
    """
    A = _terrain_operator_ungrounded(grid, eps, variant)
    return _ground_operator(A, ground_cell)


def _mms_residual_relative_l2_mfd(grid, eps, kx, kz, variant):
    """Operator-only MMS residual (interior, relative L2) for a terrain variant.

    Same measure as ``_mms_residual_relative_l2``: the RAW (un-grounded) operator
    against the consistent right-hand side eps*(kx^2+kz^2)*phi*V_cell, over
    interior cells only (grounding would inject an O(1) defect at the gauge cell's
    neighbour and contaminate the order).  The FV point·V_cell load is also the
    leading (lumped-mass) Galerkin load for the P1 Whitney operator, so the same
    measure validates both variants -- and the flat guard confirms slope ~2.0 for
    Whitney with it (no measure change was needed for the pivot).
    """
    nx, nz = grid.x.size, grid.zeta.size
    A = _terrain_operator_ungrounded(grid, eps, variant)
    phi = np.sin(kx * grid.x)[:, np.newaxis] * np.sin(kz * grid.z)
    V_cell = grid.dx * grid.dz
    b = eps * (kx ** 2 + kz ** 2) * phi * V_cell
    r = A @ phi.reshape(-1) - b.reshape(-1)
    mask = np.zeros((nx, nz), dtype=bool)
    mask[:, 1:nz - 1] = True
    mask = mask.reshape(-1)
    return np.linalg.norm(r[mask]) / np.linalg.norm(b.reshape(-1)[mask])


def terrain_mms_order_mfd(hill_fraction, variant="mfd",
                          grids=((32, 24), (64, 48), (128, 96)),
                          L=1.0, H=1.0, eps=1.0, return_details=False):
    """Finest-two MMS slope of the terrain operator (mirrors terrain_mms_order)."""
    hill_height = hill_fraction * H
    kx, kz = 2.0 * math.pi / L, math.pi / H
    errs = []
    for nx, nz in grids:
        g = terrain_grid(nx, nz, hill_height, L, H)
        errs.append(_mms_residual_relative_l2_mfd(g, eps, kx, kz, variant))
    slope = math.log(errs[-2] / errs[-1]) / math.log(grids[-1][0] / grids[-2][0])
    if return_details:
        return slope, list(grids), errs
    return slope
