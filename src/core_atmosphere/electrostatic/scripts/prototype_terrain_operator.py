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
