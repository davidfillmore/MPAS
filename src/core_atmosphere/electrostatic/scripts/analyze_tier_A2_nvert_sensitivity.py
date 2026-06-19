#!/usr/bin/env python3
"""
Tier A.2 nVertLevels sensitivity analysis for mms_sphere_horizontal.

Reads per-mesh output.nc files (produced by run_tier_A2_sphere_mms.py
--source mms_sphere_horizontal) and reports:

  1. Per-mesh L2 and per-level L2 spread (verifies horizontal isolation is exact).
  2. nVertLevels sensitivity at a fixed horizontal resolution (480km): L2 variation
     across nv4 / nv16 / nv32 runs that share the same horizontal mesh.

Usage example::

    python analyze_tier_A2_nvert_sensitivity.py \
        --run-root ~/Data/MPAS/poisson_tier_A2_scvt \
        --mesh-list 480km 240km 120km 60km \
        --nvert-variants 480km_nv4 480km_nv16 480km_nv32

The script reads output.nc from <run-root>/runs/<label>/output.nc.
"""

from __future__ import annotations

import argparse
import math
import pathlib

import netCDF4 as nc
import numpy as np


def y42(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    return np.cos(lat) ** 2 * (7.0 * np.sin(lat) ** 2 - 1.0) * np.cos(2.0 * lon)


def phi_exact(lat: np.ndarray, lon: np.ndarray, z: np.ndarray, z_top: float) -> np.ndarray:
    kz = 0.5 * np.pi / z_top
    return y42(lat, lon) * np.sin(kz * z)


def load_output(ncpath: pathlib.Path) -> dict:
    """Read phi, lat, lon, area, zgrid from output.nc; return analysis dict."""
    with nc.Dataset(ncpath) as ds:
        def strip_time(var):
            data = var[:]
            dims = list(var.dimensions)
            if dims and dims[0] == "Time":
                data = data[0]
                dims = dims[1:]
            return np.asarray(data), dims

        phi_arr, phi_dims = strip_time(ds.variables["phi"])
        lat_arr = np.asarray(ds.variables["latCell"][:]).ravel()
        lon_arr = np.asarray(ds.variables["lonCell"][:]).ravel()
        area_arr = np.asarray(ds.variables["areaCell"][:]).ravel()
        zg_arr, zg_dims = strip_time(ds.variables["zgrid"])

        # Normalise phi to (nVertLevels, nCells)
        if phi_dims[0] == "nCells":
            phi_arr = phi_arr.T

        nv = phi_arr.shape[0]
        nc_cells = phi_arr.shape[1]

        # Normalise zgrid to (nCells, nVertLevelsP1)
        if zg_arr.shape[0] == nc_cells:
            zmid = 0.5 * (zg_arr[:, :-1] + zg_arr[:, 1:]).T  # (nVertLevels, nCells)
        else:
            z1d = zg_arr.ravel()
            zmid = 0.5 * (z1d[:-1] + z1d[1:])[:, None]  # (nVertLevels, 1) broadcast

        z_top = float(np.max(zg_arr))

        # Compute per-level L2
        per_level_l2 = []
        for k in range(nv):
            z_k = zmid[k, :] if zmid.shape[1] > 1 else float(zmid[k, 0])
            exact_k = phi_exact(lat_arr, lon_arr, z_k, z_top)
            err = phi_arr[k, :] - exact_k
            num = float(np.sum(area_arr * err ** 2))
            den = float(np.sum(area_arr * exact_k ** 2))
            l2_k = math.sqrt(num / den) if den > 0 else math.sqrt(num)
            per_level_l2.append(l2_k)

        # Volume-weighted global L2
        if zmid.shape[1] > 1:
            dz = np.abs(np.diff(zg_arr, axis=1)).T  # (nVertLevels, nCells)
        else:
            z1d = zg_arr.ravel()
            dz = np.abs(np.diff(z1d))[:, None]  # (nVertLevels, 1)

        exact_all = phi_exact(lat_arr[None, :], lon_arr[None, :], zmid, z_top)
        vol = area_arr[None, :] * dz
        num_all = float(np.sum(vol * (phi_arr - exact_all) ** 2))
        den_all = float(np.sum(vol * exact_all ** 2))
        l2_global = math.sqrt(num_all / den_all) if den_all > 0 else math.sqrt(num_all)

        return {
            "nVertLevels": nv,
            "nCells": nc_cells,
            "z_top": z_top,
            "l2_global": l2_global,
            "per_level_l2": per_level_l2,
        }


def convergence_slope(entries: list[tuple[float, float]]) -> float:
    """OLS slope on the two finest entries of (h, l2) pairs sorted ascending by h."""
    usable = [(h, l2) for h, l2 in entries if math.isfinite(h) and l2 > 0]
    if len(usable) < 2:
        return math.nan
    usable = sorted(usable, key=lambda x: x[0])
    finest = usable[:2]
    return float(
        np.polyfit(
            np.log([r[0] for r in finest]),
            np.log([r[1] for r in finest]),
            1,
        )[0]
    )


def mesh_spacing_m(label: str) -> float:
    import re
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*km", label.lower())
    if m:
        return float(m.group(1)) * 1000.0
    return math.nan


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=pathlib.Path, required=True)
    parser.add_argument("--mesh-list", nargs="+", default=["480km", "240km", "120km", "60km"])
    parser.add_argument(
        "--nvert-variants",
        nargs="+",
        default=["480km_nv4", "480km_nv16", "480km_nv32"],
        help="Labels for nVertLevels sensitivity runs (must share horizontal mesh).",
    )
    args = parser.parse_args(argv)
    run_root = args.run_root.expanduser().resolve()

    # --- Section 1: convergence across horizontal resolutions ---
    print("=" * 64)
    print("Section 1: horizontal convergence (mms_sphere_horizontal)")
    print("=" * 64)
    h_pairs = []
    prev_h, prev_l2 = None, None
    for label in args.mesh_list:
        ncpath = run_root / "runs" / label / "output.nc"
        if not ncpath.exists():
            print(f"  {label}: output.nc not found, skip")
            continue
        info = load_output(ncpath)
        h = mesh_spacing_m(label)
        h_pairs.append((h, info["l2_global"]))
        lv = info["per_level_l2"]
        spread = max(lv) / min(lv) if min(lv) > 0 else math.nan
        step = (
            np.log(info["l2_global"] / prev_l2) / np.log(h / prev_h)
            if prev_h is not None
            else math.nan
        )
        print(
            f"  {label:12s}  nV={info['nVertLevels']:3d}  L2={info['l2_global']:.6e}  "
            f"per-level spread={spread:.4f}  step slope={step:.4f}"
        )
        prev_h, prev_l2 = h, info["l2_global"]

    slope = convergence_slope(h_pairs)
    print(f"  Finest-two L2 slope: {slope:.4f}")

    # --- Section 2: nVertLevels sensitivity ---
    print()
    print("=" * 64)
    print("Section 2: nVertLevels sensitivity (fixed 480km horizontal)")
    print("=" * 64)
    l2_vals = []
    for label in args.nvert_variants:
        ncpath = run_root / "runs" / label / "output.nc"
        if not ncpath.exists():
            print(f"  {label}: output.nc not found, skip")
            continue
        info = load_output(ncpath)
        l2_vals.append(info["l2_global"])
        print(
            f"  {label:15s}  nV={info['nVertLevels']:3d}  L2={info['l2_global']:.6e}"
        )
    if l2_vals:
        variation = (max(l2_vals) - min(l2_vals)) / (sum(l2_vals) / len(l2_vals))
        print(f"  L2 variation across nVertLevels variants: {variation * 100:.2f}%")
        if variation < 0.05:
            print("  => Horizontal isolation CONFIRMED: slope is vertical-independent.")
        else:
            print("  => WARNING: slope shows nVertLevels dependence — isolation not exact.")


if __name__ == "__main__":
    main()
