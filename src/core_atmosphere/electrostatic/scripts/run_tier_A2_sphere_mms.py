#!/usr/bin/env python3
"""
Tier A.2 spherical-harmonic MMS convergence runner for MPAS electrostatics.

The script prepares one run directory per global mesh under a run root, runs
atmosphere_model, computes potential-error norms from output.nc, and writes a
CSV and log-log convergence plot under run-root/results.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import sys

import netCDF4 as nc
import numpy as np


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_tier_A1_cartesian_mms as tier_a1  # noqa: E402


SOURCES = ("mms_sphere", "mms_sphere_horizontal")
OUTPUT_FIELDS = (
    "latCell",
    "lonCell",
    "zgrid",
    "areaCell",
    "phi",
    "rho_charge",
    "E_normal",
    "E_vector",
    "cg_iter_count",
    "cg_residual_initial",
    "cg_residual_final",
)


def y42(lat, lon):
    """Return the unnormalized real Y_4^2-like horizontal factor."""
    return np.cos(lat) ** 2 * (7.0 * np.sin(lat) ** 2 - 1.0) * np.cos(2.0 * lon)


def phi_exact(lat, lon, z, z_top):
    """Return the boundary-compatible Tier A.2 exact potential."""
    kz = 0.5 * np.pi / z_top
    return y42(lat, lon) * np.sin(kz * z)


def ensure_output_stream_list(stream_list_path):
    """Ensure the output stream list contains fields needed by Tier A.2."""
    if stream_list_path.exists():
        existing = stream_list_path.read_text().splitlines()
    else:
        existing = []

    present = {line.strip() for line in existing if line.strip() and not line.startswith("#")}
    lines = list(existing)
    for field in OUTPUT_FIELDS:
        if field not in present:
            lines.append(field)

    stream_list_path.write_text("\n".join(lines).rstrip() + "\n")


def require_level_cell(dataset, name, level_dim):
    """Read a level-by-cell variable independent of stored dimension order."""
    return tier_a1.require_level_cell(dataset, name, level_dim)


def require_cell(dataset, name):
    """Read a cell-centered variable."""
    return tier_a1.require_cell(dataset, name)


def compute_errors(output_nc):
    """Compute relative L2 and absolute Linf errors from one MPAS output file."""
    with nc.Dataset(output_nc) as dataset:
        phi = require_level_cell(dataset, "phi", "nVertLevels")
        lat = require_cell(dataset, "latCell")
        lon = require_cell(dataset, "lonCell")
        area = require_cell(dataset, "areaCell")
        zgrid = require_level_cell(dataset, "zgrid", "nVertLevelsP1")

        zmid = 0.5 * (zgrid[:-1, :] + zgrid[1:, :])
        dz = zgrid[1:, :] - zgrid[:-1, :]
        z_top = float(np.max(zgrid))
        if z_top <= 0.0:
            z_top = 1.0

        exact = phi_exact(lat[None, :], lon[None, :], zmid, z_top)
        volume = area[None, :] * dz

        numerator = float(np.sum(volume * (phi - exact) ** 2))
        denominator = float(np.sum(volume * exact**2))
        l2 = math.sqrt(numerator / denominator) if denominator > 0.0 else math.sqrt(numerator)
        linf = float(np.max(np.abs(phi - exact)))

        return {
            "l2": l2,
            "linf": linf,
            "cg_iter_count": int(tier_a1.scalar_diagnostic(dataset, "cg_iter_count", 0)),
            "cg_residual_initial": float(tier_a1.scalar_diagnostic(dataset, "cg_residual_initial")),
            "cg_residual_final": float(tier_a1.scalar_diagnostic(dataset, "cg_residual_final")),
            "nCells": int(phi.shape[1]),
            "nVertLevels": int(phi.shape[0]),
            "z_top": z_top,
        }


def result_stem(source):
    """Return result-file stem for a source mode."""
    if source == "mms_sphere":
        return "tier_A2_sphere_convergence"
    return "tier_A2_sphere_horizontal_convergence"


def write_csv(rows, results_dir, source):
    csv_path = results_dir / f"{result_stem(source)}.csv"
    fieldnames = (
        "mesh",
        "source",
        "h_m",
        "ranks",
        "nCells",
        "nVertLevels",
        "l2",
        "linf",
        "cg_iter_count",
        "cg_residual_initial",
        "cg_residual_final",
        "run_dir",
    )
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})
    return csv_path


def write_plot(rows, results_dir, source, l2_slope, linf_slope):
    usable = [
        row
        for row in rows
        if math.isfinite(row["h_m"]) and row["l2"] > 0.0 and row["linf"] > 0.0
    ]
    if len(usable) < 2:
        return None

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    usable = sorted(usable, key=lambda row: row["h_m"], reverse=True)
    h = np.array([row["h_m"] for row in usable])
    l2 = np.array([row["l2"] for row in usable])
    linf = np.array([row["linf"] for row in usable])

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.loglog(h, l2, "o-", label=f"L2 slope={l2_slope:.2f}")
    ax.loglog(h, linf, "s-", label=f"Linf slope={linf_slope:.2f}")
    h_ref = np.array([np.min(h), np.max(h)])
    ax.loglog(h_ref, l2[-1] * (h_ref / h[-1]) ** 2, "k--", label="2nd order")
    ax.set_xlabel("h (m)")
    ax.set_ylabel("relative error")
    ax.set_title(f"Tier A.2 spherical MMS ({source})")
    ax.legend()
    ax.grid(True, which="both", linestyle=":", linewidth=0.5)

    plot_path = results_dir / f"{result_stem(source)}.png"
    fig.savefig(plot_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    return plot_path


def build_rows(args):
    run_root = args.run_root.expanduser().resolve()
    model = args.model.expanduser().resolve()
    results_dir = run_root / "results"
    rows = []

    for mesh_name in args.mesh_list:
        run_dir, partition = tier_a1.prepare_run_dir(
            run_root,
            mesh_name,
            args.ranks,
            model,
            args.source,
        )
        namelist = run_dir / "namelist.atmosphere"
        text = tier_a1.set_namelist_value(
            namelist.read_text(),
            "electrostatic",
            "config_poisson_tol",
            f"{args.poisson_tol:.1e}",
        )
        text = tier_a1.set_namelist_value(
            text,
            "electrostatic",
            "config_poisson_max_iter",
            str(args.poisson_max_iter),
        )
        namelist.write_text(text)
        ensure_output_stream_list(run_dir / "stream_list.atmosphere.output")
        print(f"{mesh_name}: using {partition.name} with {args.ranks} ranks")

        if args.prepare_only:
            continue

        if not args.analysis_only:
            tier_a1.run_mpas(run_dir, args.ranks, args.mpiexec)

        output = run_dir / "output.nc"
        if not output.exists():
            raise FileNotFoundError(f"{output} was not created")

        result = compute_errors(output)
        row = {
            "mesh": mesh_name,
            "source": args.source,
            "h_m": tier_a1.mesh_spacing_m(mesh_name),
            "ranks": args.ranks,
            "run_dir": str(run_dir),
            **result,
        }
        rows.append(row)

        print(
            f"{mesh_name}: L2={row['l2']:.6e}, Linf={row['linf']:.6e}, "
            f"CG iter={row['cg_iter_count']}, final residual={row['cg_residual_final']:.6e}"
        )
        if row["cg_residual_final"] > args.residual_tol:
            raise RuntimeError(
                f"{mesh_name}: final CG residual {row['cg_residual_final']} "
                f"exceeds {args.residual_tol}"
            )

    return rows, results_dir


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-root",
        type=pathlib.Path,
        required=True,
        help="Root containing meshes/, runs/, and results/ directories.",
    )
    parser.add_argument("--mesh-list", nargs="+", required=True, help="Mesh labels to run.")
    parser.add_argument("--ranks", type=int, required=True, help="MPI ranks to use.")
    parser.add_argument(
        "--model",
        type=pathlib.Path,
        required=True,
        help="Path to the atmosphere_model executable.",
    )
    parser.add_argument("--mpiexec", default="mpiexec", help="MPI launcher executable.")
    parser.add_argument(
        "--source",
        choices=SOURCES,
        default="mms_sphere_horizontal",
        help=(
            "MMS source mode. mms_sphere is the coupled 3D smoke test; "
            "mms_sphere_horizontal uses the discrete vertical operator in the RHS "
            "to isolate horizontal truncation error."
        ),
    )
    parser.add_argument(
        "--analysis-only",
        action="store_true",
        help="Skip MPAS execution and analyze existing output.nc files.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare run directories but do not run or analyze output.",
    )
    parser.add_argument(
        "--residual-tol",
        type=float,
        default=1.0e-10,
        help="Maximum accepted final PCG residual.",
    )
    parser.add_argument(
        "--poisson-tol",
        type=float,
        default=1.0e-12,
        help="PCG tolerance written to namelist.atmosphere.",
    )
    parser.add_argument(
        "--poisson-max-iter",
        type=int,
        default=5000,
        help="PCG iteration cap written to namelist.atmosphere.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    rows, results_dir = build_rows(args)
    if args.prepare_only:
        print(f"Prepared {len(args.mesh_list)} run directories under {args.run_root.expanduser()}")
        return 0

    csv_path = write_csv(rows, results_dir, args.source)
    l2_slope = tier_a1.convergence_slope(rows, "l2")
    linf_slope = tier_a1.convergence_slope(rows, "linf")
    plot_path = write_plot(rows, results_dir, args.source, l2_slope, linf_slope)

    print(f"Wrote {csv_path}")
    if plot_path is not None:
        print(f"Wrote {plot_path}")
    if math.isfinite(l2_slope):
        print(f"L2 convergence slope (finest two): {l2_slope:.3f}")
    if math.isfinite(linf_slope):
        print(f"Linf convergence slope (finest two): {linf_slope:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
