#!/usr/bin/env python3
"""Terrain z-grid MMS gate for the electrostatic Poisson solver.

For each (mesh, K) pair, this runner builds one hill-independent flat bundle
(meshes/flat_<mesh>_K<K>) from the flat Tier A.1 bundle, regenerating init.nc
at K vertical levels once. Each hill amplitude then copies that flat bundle
and its init.nc into a hill-keyed bundle (meshes/terr_<mesh>_K<K>_h<h>),
rewrites zgrid to analytic cosine-hill terrain-following columns, runs the
zero-duration MMS solve with terrain_mode='zgrid', and reports interior
relative L2 errors.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import pathlib
import shutil
import sys

import netCDF4 as nc
import numpy as np


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
GATE_MIN_SLOPE = 1.8
BUFFER_LAYERS = 2
DEFAULT_SEQUENCE = ["15km:25", "7.5km:50", "3.75km:100"]
DEFAULT_HILLS = [1000.0, 6000.0]
BUNDLE_FILES = (
    "grid.nc",
    "namelist.init_atmosphere",
    "streams.init_atmosphere",
    "namelist.atmosphere",
    "streams.atmosphere",
    "stream_list.atmosphere.output",
)


def load_helper_module(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


a1 = load_helper_module("run_tier_A1_cartesian_mms")
a2 = load_helper_module("setup_tier_A2_sphere_meshes")


def phi_exact_terrain(x, y, z, x_period, y_period, z_top, hill_height):
    """Mirror electrostatic_mms_terrain_phi_exact."""
    terrain = hill_height * np.cos(2.0 * np.pi * np.asarray(x) / x_period) * np.cos(
        2.0 * np.pi * np.asarray(y) / y_period
    )
    ucoord = (np.asarray(z) - terrain) / (z_top - terrain)
    return np.sin(0.5 * np.pi * ucoord)


def parse_sequence(items):
    """Parse MESH:K entries."""
    sequence = []
    for item in items:
        mesh, separator, levels = item.partition(":")
        if not separator:
            raise SystemExit(f"--sequence entries must be MESH:K, got {item!r}")
        sequence.append((mesh, int(levels)))
    return sequence


def flat_bundle_dir(run_root, mesh, k_levels):
    """Hill-independent flat bundle location for one (mesh, K)."""
    return run_root / "meshes" / f"flat_{mesh}_K{k_levels}"


def seed_hill_bundle(flat_bundle, hill_bundle):
    """Copy run files and the flat init.nc into a hill-keyed bundle.

    rewrite_zgrid_terrain then imposes terrain on the copy; the shared flat
    init.nc is never modified, so every hill amplitude reuses one
    init_atmosphere job per (mesh, K).
    """
    hill_bundle.mkdir(parents=True, exist_ok=True)
    for name in BUNDLE_FILES:
        shutil.copy2(flat_bundle / name, hill_bundle / name)
    for pattern in ("graph.info*", "block.graph.info*"):
        for path in flat_bundle.glob(pattern):
            shutil.copy2(path, hill_bundle / path.name)
    shutil.copy2(flat_bundle / "init.nc", hill_bundle / "init.nc")


def build_variant_bundle(
    a1_mesh_dir,
    dst_dir,
    nvertlevels,
    ztop,
    init_model,
    ranks,
    mpiexec,
    force=False,
):
    """Build the hill-independent flat bundle: copy the flat A.1 bundle and
    regenerate init.nc at K levels and ztop.

    The bundle depends only on (mesh, K, ztop); the hill enters later, when
    each hill case copies this bundle (seed_hill_bundle) and rewrites the copy
    (rewrite_zgrid_terrain). init.nc here stays flat and is never modified.
    """
    if (dst_dir / "init.nc").exists() and not force:
        return

    dst_dir.mkdir(parents=True, exist_ok=True)
    for name in BUNDLE_FILES:
        shutil.copy2(a1_mesh_dir / name, dst_dir / name)

    for pattern in ("graph.info*", "block.graph.info*"):
        for path in a1_mesh_dir.glob(pattern):
            shutil.copy2(path, dst_dir / path.name)

    namelist = dst_dir / "namelist.init_atmosphere"
    text = namelist.read_text()
    text = a1.set_namelist_value(text, "dimensions", "config_nvertlevels", str(nvertlevels))
    text = a1.set_namelist_value(text, "vertical_grid", "config_ztop", f"{ztop:.1f}")
    text = a1.set_namelist_value(
        text,
        "vertical_grid",
        "config_hybrid_top_z",
        f"{ztop:.1f}",
    )
    namelist.write_text(text)

    init_model = pathlib.Path(init_model).resolve()
    if not init_model.exists():
        raise FileNotFoundError(f"init executable does not exist: {init_model}")

    a2.run_init_atmosphere(dst_dir, init_model, ranks, mpiexec, force)


def rewrite_zgrid_terrain(init_nc, *, hill_height, ztop):
    """Impose analytic terrain-following columns on init.nc, idempotently."""
    with nc.Dataset(init_nc, "r+") as dataset:
        zgrid = dataset.variables["zgrid"]
        dims = zgrid.dimensions
        level_axis = dims.index("nVertLevelsP1")
        nlevels_p1 = zgrid.shape[level_axis]
        xcell = np.asarray(dataset.variables["xCell"][:])
        ycell = np.asarray(dataset.variables["yCell"][:])
        x_period = a1.positive_attr(dataset, "x_period")
        y_period = a1.positive_attr(dataset, "y_period")
        if not math.isfinite(x_period):
            x_period = float(np.max(xcell) - np.min(xcell))
        if not math.isfinite(y_period):
            y_period = float(np.max(ycell) - np.min(ycell))

        terrain = hill_height * np.cos(2.0 * np.pi * xcell / x_period) * np.cos(
            2.0 * np.pi * ycell / y_period
        )
        zeta = np.linspace(0.0, ztop, nlevels_p1)

        if level_axis == 1:
            znew = zeta[None, :] + terrain[:, None] * (1.0 - zeta[None, :] / ztop)
        else:
            znew = zeta[:, None] + terrain[None, :] * (1.0 - zeta[:, None] / ztop)
        zgrid[:] = znew


def add_terrain_namelist_keys(run_dir, hill_height):
    namelist = run_dir / "namelist.atmosphere"
    text = namelist.read_text()
    text = a1.set_namelist_value(
        text,
        "electrostatic",
        "config_electrostatic_terrain_mode",
        "'zgrid'",
    )
    text = a1.set_namelist_value(
        text,
        "electrostatic",
        "config_electrostatic_hill_height",
        f"{hill_height:.3f}",
    )
    namelist.write_text(text)


def compute_interior_errors(output_nc, hill_height, buffer_layers=BUFFER_LAYERS):
    """Compute full and interior volume-weighted relative L2 errors."""
    with nc.Dataset(output_nc) as dataset:
        phi = a1.require_level_cell(dataset, "phi", "nVertLevels")
        zgrid = a1.require_level_cell(dataset, "zgrid", "nVertLevelsP1")
        area = a1.require_cell(dataset, "areaCell")
        xcell = a1.require_cell(dataset, "xCell")
        ycell = a1.require_cell(dataset, "yCell")
        x_period, y_period, z_top = a1.infer_periods(dataset, xcell, ycell, zgrid)
        cg_iter_count = int(a1.scalar_diagnostic(dataset, "cg_iter_count", 0))
        cg_residual_initial = float(a1.scalar_diagnostic(dataset, "cg_residual_initial"))
        cg_residual_final = float(a1.scalar_diagnostic(dataset, "cg_residual_final"))

    zmid = 0.5 * (zgrid[:-1, :] + zgrid[1:, :])
    dz = zgrid[1:, :] - zgrid[:-1, :]
    z_bottom = float(np.min(zgrid[0, :]))
    dz_p = (z_top - z_bottom) / phi.shape[0]
    terrain = hill_height * np.cos(2.0 * np.pi * xcell / x_period) * np.cos(
        2.0 * np.pi * ycell / y_period
    )
    exact = phi_exact_terrain(
        xcell[None, :],
        ycell[None, :],
        zmid,
        x_period,
        y_period,
        z_top,
        hill_height,
    )
    volume = area[None, :] * dz
    err2 = volume * (phi - exact) ** 2
    ref2 = volume * exact**2
    interior = zmid >= (terrain[None, :] + buffer_layers * dz_p)

    def relative_l2(mask):
        numerator = float(err2[mask].sum())
        denominator = float(ref2[mask].sum())
        return math.sqrt(numerator / denominator) if denominator > 0.0 else float("inf")

    return {
        "l2_interior": relative_l2(interior),
        "l2_full": relative_l2(np.ones_like(interior, dtype=bool)),
        "n_interior": int(interior.sum()),
        "cg_iter_count": cg_iter_count,
        "cg_residual_initial": cg_residual_initial,
        "cg_residual_final": cg_residual_final,
    }


def check_residual(label, residual, tol):
    """Abort unless the final CG residual is a finite value within tol.

    Written as `not (residual <= tol)` so NaN — from a diverged solve or a
    missing cg_residual_final diagnostic (scalar_diagnostic defaults to
    NaN) — fails the gate instead of slipping past `residual > tol`.
    """
    if not (residual <= tol):
        raise RuntimeError(
            f"{label}: final CG residual {residual:.3e} is not <= {tol:.1e}"
        )


def finest_two_slope(rows, key):
    usable = [
        row
        for row in rows
        if math.isfinite(row["h_m"]) and row[key] > 0.0 and math.isfinite(row[key])
    ]
    if len(usable) < 2:
        return math.nan

    usable = sorted(usable, key=lambda row: row["h_m"])
    finest = usable[:2]
    if finest[0]["h_m"] == finest[1]["h_m"]:
        raise SystemExit(
            "finest_two_slope: the two finest usable rows share h_m = "
            f"{finest[0]['h_m']:g} m "
            f"({finest[0].get('mesh', '?')} and {finest[1].get('mesh', '?')}); "
            "a combined-refinement slope needs two distinct horizontal "
            "spacings (the finest two usable rows share h_m, so no "
            "combined-refinement slope can be formed)"
        )
    return math.log(finest[1][key] / finest[0][key]) / math.log(
        finest[1]["h_m"] / finest[0]["h_m"]
    )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=pathlib.Path, required=True)
    parser.add_argument(
        "--a1-root",
        type=pathlib.Path,
        required=True,
        help="Tier A.1 root containing flat mesh bundles under meshes/.",
    )
    parser.add_argument("--model", type=pathlib.Path, required=True)
    parser.add_argument("--init-model", type=pathlib.Path, required=True)
    parser.add_argument("--ranks", type=int, required=True)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument(
        "--sequence",
        nargs="+",
        default=DEFAULT_SEQUENCE,
        help="Combined-refinement MESH:K pairs, coarse to fine.",
    )
    parser.add_argument(
        "--hill-heights",
        nargs="+",
        type=float,
        default=DEFAULT_HILLS,
        help="Hill amplitudes in meters; the first is the hard gate.",
    )
    parser.add_argument("--ztop", type=float, default=20000.0)
    parser.add_argument("--residual-tol", type=float, default=1.0e-8)
    parser.add_argument("--analysis-only", action="store_true")
    parser.add_argument("--force", action="store_true", help="Rebuild variant bundles.")
    return parser.parse_args(argv)


def run_case(args, run_root, mesh, k_levels, hill_height):
    label = f"terr_{mesh}_K{k_levels}_h{int(round(hill_height))}"
    bundle = run_root / "meshes" / label
    if not args.analysis_only:
        seed_hill_bundle(flat_bundle_dir(run_root, mesh, k_levels), bundle)
        rewrite_zgrid_terrain(bundle / "init.nc", hill_height=hill_height, ztop=args.ztop)

    run_dir, partition = a1.prepare_run_dir(
        run_root,
        label,
        args.ranks,
        args.model.expanduser().resolve(),
        "mms_terrain",
    )
    add_terrain_namelist_keys(run_dir, hill_height)
    print(f"{label}: using {partition.name} with {args.ranks} ranks")

    if not args.analysis_only:
        a1.run_mpas(run_dir, args.ranks, args.mpiexec)

    output_nc = run_dir / "output.nc"
    if not output_nc.exists():
        raise FileNotFoundError(f"{output_nc} was not created")

    result = compute_interior_errors(output_nc, hill_height)
    check_residual(label, result["cg_residual_final"], args.residual_tol)

    return {
        "mesh": mesh,
        "K": k_levels,
        "h_m": a1.mesh_spacing_m(mesh),
        "hill_m": hill_height,
        "run_dir": str(run_dir),
        **result,
    }


def write_rows(rows, run_root, hill_height):
    csv_path = run_root / "results" / f"terrain_zgrid_h{int(round(hill_height))}.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return csv_path


def main(argv=None):
    args = parse_args(argv)
    run_root = args.run_root.expanduser().resolve()
    a1_root = args.a1_root.expanduser().resolve()
    (run_root / "results").mkdir(parents=True, exist_ok=True)
    sequence = parse_sequence(args.sequence)
    gate_failed = False

    if not args.analysis_only:
        for mesh, k_levels in sequence:
            build_variant_bundle(
                a1_root / "meshes" / mesh,
                flat_bundle_dir(run_root, mesh, k_levels),
                k_levels,
                args.ztop,
                args.init_model.expanduser(),
                args.ranks,
                args.mpiexec,
                force=args.force,
            )

    for hill_index, hill_height in enumerate(args.hill_heights):
        rows = []
        for mesh, k_levels in sequence:
            row = run_case(args, run_root, mesh, k_levels, hill_height)
            rows.append(row)
            print(
                f"{row['mesh']}/K{row['K']}/h{hill_height:.0f}: "
                f"L2int={row['l2_interior']:.6e}, L2full={row['l2_full']:.6e}, "
                f"CG iter={row['cg_iter_count']}, "
                f"final residual={row['cg_residual_final']:.6e}"
            )

        slope = finest_two_slope(rows, "l2_interior")
        csv_path = write_rows(rows, run_root, hill_height)
        print(
            f"hill {hill_height:.0f} m: interior L2 slope (finest two) = "
            f"{slope:.3f}"
        )
        print(f"Wrote {csv_path}")

        if hill_index == 0:
            if slope >= GATE_MIN_SLOPE:
                print(f"GATE PASS: gentle-hill slope {slope:.3f} >= {GATE_MIN_SLOPE}")
            else:
                print(f"GATE FAIL: gentle-hill slope {slope:.3f} < {GATE_MIN_SLOPE}")
                gate_failed = True
        elif slope >= GATE_MIN_SLOPE:
            print(f"SOFT PASS: steep-hill slope {slope:.3f} >= {GATE_MIN_SLOPE}")
        else:
            print(
                f"SOFT MISS: steep-hill slope {slope:.3f} < {GATE_MIN_SLOPE}; "
                "recorded for the terrain-boundary upgrade discussion"
            )

    return 1 if gate_failed else 0


if __name__ == "__main__":
    sys.exit(main())
