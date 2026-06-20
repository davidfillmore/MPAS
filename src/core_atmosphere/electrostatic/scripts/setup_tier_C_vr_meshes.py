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


def generate_vr_mesh(bundle_dir, work_dir, *, fine_km, coarse_km, patch, dlat_deg,
                     earth_radius, plot_cell_width, force):
    """Generate grid.nc + graph.info for one VR mesh via build_spherical_mesh."""
    from mpas_tools.mesh.creation.build_mesh import build_spherical_mesh
    grid, graph = bundle_dir / "grid.nc", bundle_dir / "graph.info"
    if grid.exists() and graph.exists() and not force:
        return
    work_dir.mkdir(parents=True, exist_ok=True)
    a2.clean_generated_mesh_files(work_dir)
    cw, lon, lat = build_cell_width_grid(fine_km=fine_km, coarse_km=coarse_km,
                                         dlat_deg=dlat_deg, **patch)
    prev = pathlib.Path.cwd()
    try:
        os.chdir(work_dir)
        build_spherical_mesh(cw, lon, lat, earth_radius=earth_radius,
                             out_filename="base_mesh.nc", plot_cellWidth=plot_cell_width)
    finally:
        os.chdir(prev)
    shutil.copy2(work_dir / "base_mesh.nc", grid)
    shutil.copy2(work_dir / "graph.info", graph)


def setup_vr_mesh(args, label):
    sizes = VR_MESHES[label]
    patch = dict(lat0=args.lat0, lon0=args.lon0, r_in_deg=args.r_in_deg, r_out_deg=args.r_out_deg)
    run_root = args.run_root.expanduser().resolve()
    bundle = run_root / "meshes" / label
    work = run_root / "mesh_scratch" / label
    bundle.mkdir(parents=True, exist_ok=True)
    print(f"{label}: generating VR mesh fine={sizes['fine_km']} coarse={sizes['coarse_km']} km")
    generate_vr_mesh(bundle, work, fine_km=sizes["fine_km"], coarse_km=sizes["coarse_km"],
                     patch=patch, dlat_deg=args.dlat_deg, earth_radius=args.earth_radius,
                     plot_cell_width=args.plot_cell_width, force=args.force)
    a2.write_bundle_inputs(bundle, args.nvertlevels, args.ztop, args.force)
    a2.partition_mesh(bundle, args.ranks, args.force)
    if not args.mesh_only:
        init_model = args.init_model.expanduser().resolve()
        if not init_model.exists():
            raise FileNotFoundError(f"init model does not exist: {init_model}")
        a2.run_init_atmosphere(bundle, init_model, args.ranks, args.mpiexec, args.force)
    # Pentagon-vs-patch audit
    pidx, plat, plon = find_pentagons(bundle / "grid.nc")
    pdist = great_circle_deg(plat, plon, patch["lat0"], patch["lon0"])
    with nc.Dataset(bundle / "grid.nc") as ds:
        n_cells = len(ds.dimensions["nCells"])
    return {
        "mesh": label, **sizes, **patch, "nCells": n_cells,
        "n_pentagons": int(len(pidx)),
        "pentagon_dist_deg": [round(float(x), 3) for x in np.sort(pdist)],
        "pentagons_in_transition": int(np.sum((pdist > patch["r_in_deg"]) & (pdist < patch["r_out_deg"]))),
        "bundle": str(bundle),
    }


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-root", type=pathlib.Path, required=True)
    p.add_argument("--mesh-list", nargs="+", default=["vr_480_120"], choices=list(VR_MESHES))
    p.add_argument("--ranks", type=int, default=8)
    p.add_argument("--init-model", type=pathlib.Path, default=pathlib.Path("./init_atmosphere_model"))
    p.add_argument("--mpiexec", default="mpiexec")
    p.add_argument("--nvertlevels", type=int, default=16)
    p.add_argument("--ztop", type=float, default=20000.0)
    p.add_argument("--earth-radius", type=float, default=EARTH_RADIUS_M)
    p.add_argument("--dlat-deg", type=float, default=1.0)
    p.add_argument("--lat0", type=float, default=DEFAULT_PATCH["lat0"])
    p.add_argument("--lon0", type=float, default=DEFAULT_PATCH["lon0"])
    p.add_argument("--r-in-deg", type=float, default=DEFAULT_PATCH["r_in_deg"])
    p.add_argument("--r-out-deg", type=float, default=DEFAULT_PATCH["r_out_deg"])
    p.add_argument("--mesh-only", action="store_true")
    p.add_argument("--plot-cell-width", action="store_true")
    p.add_argument("--force", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    run_root = args.run_root.expanduser().resolve()
    for sub in ("meshes", "runs", "results"):
        (run_root / sub).mkdir(parents=True, exist_ok=True)
    summaries = [setup_vr_mesh(args, label) for label in args.mesh_list]
    manifest = {"mesh_type": "variable_circular", "earth_radius_m": args.earth_radius,
                "ztop_m": args.ztop, "patch": {k: getattr(args, k) for k in
                ("lat0", "lon0", "r_in_deg", "r_out_deg", "dlat_deg")}, "meshes": summaries}
    path = run_root / "mesh_manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {path}")
    for s in summaries:
        print(f"{s['mesh']}: nCells={s['nCells']} pentagons={s['n_pentagons']} "
              f"in_transition={s['pentagons_in_transition']} dists={s['pentagon_dist_deg'][:3]}…")
    return 0


if __name__ == "__main__":
    sys.exit(main())
