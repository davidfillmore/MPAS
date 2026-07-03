#!/usr/bin/env python3
"""Prepare and run a terrain charge-coupled supercell case."""

from __future__ import annotations

import argparse
import pathlib
import shutil
import sys

import netCDF4 as nc
import numpy as np


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_charge_coupled_supercell as charge_supercell  # noqa: E402
import run_tier_A1_cartesian_mms as tier_a1  # noqa: E402
from run_tripole_supercell import seed_run_dir_from_template  # noqa: E402


DEFAULT_RUN_DIR = pathlib.Path(
    "~/Data/MPAS/poisson_charge_coupled_supercell_terrain_h1000/run"
)
DEFAULT_TEMPLATE_RUN_DIR = pathlib.Path("~/Data/MPAS/supercell")
DEFAULT_HILL_HEIGHT = 1000.0
DEFAULT_HILL_HALF_WIDTH_KM = 20.0


def cosine_bell_terrain(xcell, ycell, *, hill_height, half_width_m):
    x0 = 0.5 * (float(np.nanmin(xcell)) + float(np.nanmax(xcell)))
    y0 = 0.5 * (float(np.nanmin(ycell)) + float(np.nanmax(ycell)))
    radius = np.sqrt((xcell - x0) ** 2 + (ycell - y0) ** 2)
    terrain = np.zeros_like(xcell, dtype=float)
    inside = radius < half_width_m
    terrain[inside] = 0.5 * hill_height * (
        1.0 + np.cos(np.pi * radius[inside] / half_width_m)
    )
    return terrain


def rewrite_supercell_init_terrain(init_nc, *, hill_height, hill_half_width_m):
    """Rewrite supercell_init.nc with monotone terrain-following zgrid columns."""
    if not np.isfinite(hill_half_width_m) or hill_half_width_m <= 0.0:
        raise ValueError(
            f"terrain half-width must be finite and positive, got {hill_half_width_m}"
        )
    if not np.isfinite(hill_height):
        raise ValueError(f"terrain hill height must be finite, got {hill_height}")

    with nc.Dataset(init_nc, "r+") as dataset:
        xcell = np.asarray(dataset.variables["xCell"][:], dtype=float)
        ycell = np.asarray(dataset.variables["yCell"][:], dtype=float)
        zgrid = dataset.variables["zgrid"]
        dims = zgrid.dimensions
        if "nCells" not in dims or "nVertLevelsP1" not in dims:
            raise ValueError(f"zgrid has unsupported dimensions {dims}")

        cell_axis = dims.index("nCells")
        level_axis = dims.index("nVertLevelsP1")
        if cell_axis == level_axis:
            raise ValueError(f"zgrid has unsupported dimensions {dims}")

        n_cells = zgrid.shape[cell_axis]
        nlevels_p1 = zgrid.shape[level_axis]
        if n_cells != xcell.size or n_cells != ycell.size:
            raise ValueError(
                f"zgrid nCells={n_cells} does not match xCell/yCell size "
                f"{xcell.size}/{ycell.size}"
            )

        ztop = float(np.nanmax(np.asarray(zgrid[:], dtype=float)))
        if ztop <= 0.0:
            raise ValueError(f"zgrid top must be positive, got {ztop}")

        terrain = cosine_bell_terrain(
            xcell,
            ycell,
            hill_height=hill_height,
            half_width_m=hill_half_width_m,
        )
        terrain_min = float(np.nanmin(terrain))
        terrain_max = float(np.nanmax(terrain))
        if (
            not np.isfinite(terrain_min)
            or not np.isfinite(terrain_max)
            or terrain_min < 0.0
            or terrain_max >= ztop
        ):
            raise ValueError(
                "terrain range must be finite and satisfy "
                f"0 <= terrain <= terrain_max < ztop; got "
                f"[{terrain_min}, {terrain_max}] with ztop {ztop}"
            )

        zeta = np.linspace(0.0, ztop, nlevels_p1)
        znew = terrain[:, None] + zeta[None, :] * (ztop - terrain[:, None]) / ztop

        if cell_axis == 0 and level_axis == 1:
            zgrid[:] = znew
        elif level_axis == 0 and cell_axis == 1:
            zgrid[:] = znew.T
        else:
            raise ValueError(f"zgrid has unsupported dimensions {dims}")

        if "ter" in dataset.variables:
            dataset.variables["ter"][:] = terrain

    return {
        "nCells": int(n_cells),
        "ztop": float(ztop),
        "terrain_min": terrain_min,
        "terrain_max": terrain_max,
    }


def seed_terrain_run_dir(template_run_dir, run_dir):
    """Seed a terrain run directory, then force supercell_init.nc to be a copy."""
    template_run_dir = template_run_dir.expanduser().resolve()
    run_dir = run_dir.expanduser().resolve()
    if template_run_dir == run_dir:
        raise ValueError(
            f"template_run_dir and run_dir resolve to the same directory: {run_dir}"
        )

    seed_run_dir_from_template(template_run_dir, run_dir)
    source = template_run_dir / "supercell_init.nc"
    target = run_dir / "supercell_init.nc"
    if target.exists() or target.is_symlink():
        target.unlink()
    shutil.copy2(source, target)


def add_terrain_namelist_keys(run_dir, hill_height):
    """Enable terrain zgrid electrostatics in namelist.atmosphere."""
    namelist = run_dir / "namelist.atmosphere"
    text = namelist.read_text()
    text = tier_a1.set_namelist_value(
        text,
        "electrostatic",
        "config_electrostatic_terrain_mode",
        "'zgrid'",
    )
    text = tier_a1.set_namelist_value(
        text,
        "electrostatic",
        "config_electrostatic_hill_height",
        f"{hill_height:.3f}",
    )
    namelist.write_text(text)


def prepare_run_dir(
    *,
    template_run_dir,
    run_dir,
    model,
    ranks,
    run_duration,
    interval,
    stub_alpha,
    stub_beta,
    poisson_tol,
    poisson_max_iter,
    solve_at_init,
    output_interval,
    hill_height,
    hill_half_width_km,
):
    """Create and configure an isolated terrain dynamic-coupling run directory."""
    seed_terrain_run_dir(template_run_dir, run_dir)
    run_dir = run_dir.expanduser().resolve()
    stats = rewrite_supercell_init_terrain(
        run_dir / "supercell_init.nc",
        hill_height=hill_height,
        hill_half_width_m=hill_half_width_km * 1000.0,
    )
    tier_a1.symlink_model(model.expanduser(), run_dir)

    namelist = run_dir / "namelist.atmosphere"
    streams = run_dir / "streams.atmosphere"
    if not namelist.exists() or not streams.exists():
        raise FileNotFoundError(
            f"{run_dir} needs namelist.atmosphere and streams.atmosphere; "
            f"seed it from a complete supercell run directory"
        )

    charge_supercell.configure_namelist(
        namelist,
        run_duration=run_duration,
        interval=interval,
        stub_alpha=stub_alpha,
        stub_beta=stub_beta,
        poisson_tol=poisson_tol,
        poisson_max_iter=poisson_max_iter,
        solve_at_init=solve_at_init,
    )
    add_terrain_namelist_keys(run_dir, hill_height)
    charge_supercell.configure_streams(streams, output_interval)
    charge_supercell.ensure_coupled_output_stream_list(
        run_dir / "stream_list.atmosphere.output"
    )
    partition = tier_a1.find_partition_file(run_dir, ranks, namelist.read_text())
    if (run_dir / "supercell_init.nc").is_symlink():
        raise RuntimeError(f"{run_dir / 'supercell_init.nc'} must be a real copied file")
    return run_dir, partition, stats


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--template-run-dir",
        type=pathlib.Path,
        default=DEFAULT_TEMPLATE_RUN_DIR,
        help="Existing supercell run directory to seed inputs from.",
    )
    parser.add_argument(
        "--run-dir",
        type=pathlib.Path,
        default=DEFAULT_RUN_DIR,
        help="Isolated terrain run directory to create or reuse.",
    )
    parser.add_argument(
        "--model",
        type=pathlib.Path,
        default=pathlib.Path("./atmosphere_model"),
        help="Path to the atmosphere_model executable.",
    )
    parser.add_argument("--ranks", type=int, default=8, help="MPI ranks to use.")
    parser.add_argument("--mpiexec", default="mpiexec", help="MPI launcher executable.")
    parser.add_argument(
        "--run-duration",
        default="00_00:10:00",
        help="MPAS config_run_duration value.",
    )
    parser.add_argument(
        "--electrostatic-interval",
        type=float,
        default=60.0,
        help="Seconds between diagnostic electrostatic solves.",
    )
    parser.add_argument("--stub-alpha", type=float, default=1.0e-8)
    parser.add_argument("--stub-beta", type=float, default=1.0e-8)
    parser.add_argument("--poisson-tol", type=float, default=1.0e-10)
    parser.add_argument("--poisson-max-iter", type=int, default=5000)
    parser.add_argument(
        "--output-interval",
        default="00:10:00",
        help="MPAS output stream interval for output.nc.",
    )
    parser.add_argument(
        "--solve-at-init",
        action="store_true",
        help="Also run a diagnostic electrostatic solve during initialization.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare the run directory but do not run MPAS.",
    )
    parser.add_argument("--hill-height", type=float, default=DEFAULT_HILL_HEIGHT)
    parser.add_argument(
        "--hill-half-width-km",
        type=float,
        default=DEFAULT_HILL_HALF_WIDTH_KM,
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    run_dir, partition, stats = prepare_run_dir(
        template_run_dir=args.template_run_dir,
        run_dir=args.run_dir,
        model=args.model,
        ranks=args.ranks,
        run_duration=args.run_duration,
        interval=args.electrostatic_interval,
        stub_alpha=args.stub_alpha,
        stub_beta=args.stub_beta,
        poisson_tol=args.poisson_tol,
        poisson_max_iter=args.poisson_max_iter,
        solve_at_init=args.solve_at_init,
        output_interval=args.output_interval,
        hill_height=args.hill_height,
        hill_half_width_km=args.hill_half_width_km,
    )
    print(f"Prepared {run_dir}")
    print(f"Using partition {partition}")
    print(
        "Terrain range "
        f"{stats['terrain_min']:.3f} to {stats['terrain_max']:.3f} m; "
        f"ztop {stats['ztop']:.3f} m"
    )
    if args.prepare_only:
        return 0
    charge_supercell.run_mpas(run_dir, args.ranks, args.mpiexec)
    print(f"Wrote {run_dir / 'output.nc'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
