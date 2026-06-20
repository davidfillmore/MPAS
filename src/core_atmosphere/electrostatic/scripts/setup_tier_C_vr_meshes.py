#!/usr/bin/env python3
"""Generate Tier C variable-resolution spherical MMS mesh bundles.

A circular refinement patch (fine core, smooth transition annulus, coarse
far-field) fed to MPAS-Tools build_spherical_mesh as a variable cellWidth.
Each bundle matches the Tier A.2 layout so run_tier_A2_sphere_mms.py runs
against it unchanged. Reuses the A.2 namelist/streams/init/partition plumbing.
"""
from __future__ import annotations

import argparse, json, math, os, pathlib, shutil, sys
import numpy as np
import netCDF4 as nc

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
import setup_tier_A2_sphere_meshes as a2  # reuse namelist/streams/init/partition

EARTH_RADIUS_M = 6371229.0

# Nested 4x VR triple: same region, (fine,coarse) halved each step.
VR_MESHES = {
    "vr_480_120": dict(fine_km=120.0, coarse_km=480.0),
    "vr_240_60":  dict(fine_km=60.0,  coarse_km=240.0),
    "vr_120_30":  dict(fine_km=30.0,  coarse_km=120.0),
}
DEFAULT_PATCH = dict(lat0=0.0, lon0=0.0, r_in_deg=14.0, r_out_deg=34.0)


def great_circle_deg(lat, lon, lat0, lon0):
    """Haversine angular distance in degrees; vectorized over lat/lon arrays."""
    la, lo = np.radians(lat), np.radians(lon)
    la0, lo0 = math.radians(lat0), math.radians(lon0)
    dlat, dlon = la - la0, lo - lo0
    a = np.sin(dlat / 2.0) ** 2 + np.cos(la0) * np.cos(la) * np.sin(dlon / 2.0) ** 2
    return np.degrees(2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))))


def smoothstep(d, r_in, r_out):
    """Cosine ramp 0->1 across [r_in, r_out]; flat outside."""
    t = np.clip((np.asarray(d, dtype=float) - r_in) / (r_out - r_in), 0.0, 1.0)
    return 0.5 * (1.0 - np.cos(np.pi * t))


def vr_cell_width(lat2d, lon2d, *, fine_km, coarse_km, lat0, lon0, r_in_deg, r_out_deg):
    """cellWidth(km) for the circular refinement patch."""
    d = great_circle_deg(lat2d, lon2d, lat0, lon0)
    return fine_km + (coarse_km - fine_km) * smoothstep(d, r_in_deg, r_out_deg)


def build_cell_width_grid(*, fine_km, coarse_km, lat0, lon0, r_in_deg, r_out_deg, dlat_deg=1.0):
    """Return (cell_width[nlat,nlon], lon_1d, lat_1d) for build_spherical_mesh."""
    lat = np.linspace(-90.0, 90.0, int(round(180.0 / dlat_deg)) + 1)
    lon = np.linspace(-180.0, 180.0, int(round(360.0 / dlat_deg)) + 1)
    lon2d, lat2d = np.meshgrid(lon, lat)
    cw = vr_cell_width(lat2d, lon2d, fine_km=fine_km, coarse_km=coarse_km,
                       lat0=lat0, lon0=lon0, r_in_deg=r_in_deg, r_out_deg=r_out_deg)
    return cw, lon, lat


def find_pentagons(grid_nc):
    """Return (idx, lat_deg, lon_deg) for cells with nEdgesOnCell == 5."""
    with nc.Dataset(grid_nc) as ds:
        ne = np.asarray(ds.variables["nEdgesOnCell"][:])
        lat = np.degrees(np.asarray(ds.variables["latCell"][:]))
        lon = np.degrees(np.asarray(ds.variables["lonCell"][:]))
    idx = np.where(ne == 5)[0]
    return idx, lat[idx], lon[idx]
