#!/usr/bin/env python3
"""3-D terrain prototype: triangular-lattice horizontal x terrain-following layers."""
from __future__ import annotations
import math
from types import SimpleNamespace
import numpy as np

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
