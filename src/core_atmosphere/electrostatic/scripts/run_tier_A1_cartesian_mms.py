#!/usr/bin/env python3
"""
Tier A.1 Cartesian MMS convergence runner for MPAS electrostatics.

The script prepares one run directory per mesh under a run root, runs
atmosphere_model, computes potential-error norms from output.nc, and writes a
CSV and log-log convergence plot under run-root/results.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

import netCDF4 as nc
import numpy as np


ELECTROSTATIC_CONFIG = (
    ("config_electrostatic_enable", ".true."),
    ("config_electrostatic_solve_at_init", ".true."),
    ("config_electrostatic_source", "'mms_cart'"),
    ("config_poisson_preconditioner", "'jacobi'"),
    ("config_poisson_tol", "1.0e-8"),
    ("config_poisson_max_iter", "2000"),
    ("config_electrostatic_bc_ground", "0.0"),
)

REQUIRED_OUTPUT_FIELDS = (
    "xCell",
    "yCell",
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

STALE_PATTERNS = (
    "output.nc",
    "run.out",
    "log.atmosphere.*.out",
    "log.atmosphere.*.err",
)


def phi_exact(x, y, z, x_period, y_period, z_top):
    """Return the corrected Cartesian MMS potential."""
    kx = 2.0 * np.pi / x_period
    ky = 2.0 * np.pi / y_period
    kz = 0.5 * np.pi / z_top
    return np.sin(kx * x) * np.sin(ky * y) * np.sin(kz * z)


def mesh_spacing_m(mesh_name):
    """Infer nominal mesh spacing in meters from labels such as 15km."""
    name = mesh_name.lower().replace("_", ".")
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*km", name)
    if match:
        return float(match.group(1)) * 1000.0
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*m", name)
    if match:
        return float(match.group(1))
    return math.nan


def configure_namelist_text(text):
    """Return namelist text configured for the zero-duration MMS solve."""
    text = set_namelist_value(text, "nhyd_model", "config_run_duration", "'00_00:00:00'")
    return replace_namelist_block(text, "electrostatic", ELECTROSTATIC_CONFIG)


def set_namelist_value(text, block_name, key, value):
    """Set one namelist key, preserving the rest of the block when possible."""
    pattern = re.compile(rf"(?ms)^&{re.escape(block_name)}\b.*?^\s*/\s*$")
    match = pattern.search(text)
    if not match:
        return text.rstrip() + f"\n\n&{block_name}\n    {key} = {value}\n/\n"

    lines = match.group(0).splitlines()
    key_pattern = re.compile(rf"^\s*{re.escape(key)}\s*=")
    filtered = [line for line in lines if not key_pattern.match(line)]
    slash_index = next(
        (index for index, line in enumerate(filtered) if line.strip() == "/"),
        len(filtered),
    )
    filtered.insert(slash_index, f"    {key} = {value}")
    return text[: match.start()] + "\n".join(filtered) + text[match.end() :]


def replace_namelist_block(text, block_name, entries):
    """Replace or append a complete namelist block."""
    block = "\n".join(
        [f"&{block_name}"] + [f"    {key} = {value}" for key, value in entries] + ["/"]
    )
    pattern = re.compile(rf"(?ms)^&{re.escape(block_name)}\b.*?^\s*/\s*$")
    if pattern.search(text):
        return pattern.sub(block, text, count=1)
    return text.rstrip() + "\n\n" + block + "\n"


def ensure_streams_netcdf(streams_path):
    """Force all streams to netCDF I/O and make output.nc an initial-only stream."""
    tree = ET.parse(streams_path)
    root = tree.getroot()

    for elem in root:
        if elem.tag in {"stream", "immutable_stream"}:
            elem.set("io_type", "netcdf")

    output = root.find("./stream[@name='output']")
    if output is None:
        output = ET.SubElement(root, "stream", {"name": "output", "type": "output"})

    output.set("type", "output")
    output.set("io_type", "netcdf")
    output.set("filename_template", "output.nc")
    output.set("filename_interval", "none")
    output.set("output_interval", "initial_only")

    if output.find("./file[@name='stream_list.atmosphere.output']") is None:
        ET.SubElement(output, "file", {"name": "stream_list.atmosphere.output"})

    ET.indent(tree, space="  ")
    tree.write(streams_path, encoding="unicode")


def ensure_output_stream_list(stream_list_path):
    """Ensure the output stream list contains the variables needed for analysis."""
    if stream_list_path.exists():
        existing = stream_list_path.read_text().splitlines()
    else:
        existing = []

    present = {line.strip() for line in existing if line.strip() and not line.startswith("#")}
    lines = list(existing)
    for field in REQUIRED_OUTPUT_FIELDS:
        if field not in present:
            lines.append(field)

    stream_list_path.write_text("\n".join(lines).rstrip() + "\n")


def configure_namelist(path):
    path.write_text(configure_namelist_text(path.read_text()))


def seed_run_dir_from_mesh_bundle(mesh_dir, run_dir):
    """Symlink mesh-bundle contents into the run directory without overwriting inputs."""
    if not mesh_dir.exists():
        return

    for source in mesh_dir.iterdir():
        if source.name == "atmosphere_model" or source.name in {"output.nc", "run.out"}:
            continue
        if source.name.startswith("log.atmosphere."):
            continue

        target = run_dir / source.name
        if target.exists() or target.is_symlink():
            continue
        os.symlink(source.resolve(), target)


def symlink_model(model, run_dir):
    target = run_dir / "atmosphere_model"
    model = model.resolve()
    if not model.exists():
        raise FileNotFoundError(f"model executable does not exist: {model}")

    if target.is_symlink() or not target.exists():
        if target.exists() or target.is_symlink():
            target.unlink()
        os.symlink(model, target)
    elif target.resolve() != model:
        raise FileExistsError(f"{target} exists and is not a symlink to {model}")


def parse_namelist_value(text, key):
    match = re.search(rf"^\s*{re.escape(key)}\s*=\s*([^!\n]+)", text, flags=re.MULTILINE)
    if not match:
        return None
    return match.group(1).strip().strip("'\"")


def find_partition_file(run_dir, ranks, namelist_text):
    prefix = parse_namelist_value(namelist_text, "config_block_decomp_file_prefix")
    if prefix:
        candidate = run_dir / f"{prefix}{ranks}"
        if candidate.exists():
            return candidate

    candidates = sorted(run_dir.glob(f"*.graph.info.part.{ranks}"))
    if candidates:
        return candidates[0]

    raise FileNotFoundError(
        f"No partition file matching rank count {ranks} in {run_dir}; "
        f"expected something like *.graph.info.part.{ranks}"
    )


def prepare_run_dir(run_root, mesh_name, ranks, model):
    """Create and configure the run directory for one mesh."""
    run_dir = run_root / "runs" / mesh_name
    mesh_dir = run_root / "meshes" / mesh_name
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_root / "results").mkdir(parents=True, exist_ok=True)

    seed_run_dir_from_mesh_bundle(mesh_dir, run_dir)
    symlink_model(model, run_dir)

    namelist = run_dir / "namelist.atmosphere"
    streams = run_dir / "streams.atmosphere"
    if not namelist.exists() or not streams.exists():
        raise FileNotFoundError(
            f"{run_dir} needs namelist.atmosphere and streams.atmosphere; "
            f"place them there or in {mesh_dir}"
        )

    configure_namelist(namelist)
    ensure_streams_netcdf(streams)
    ensure_output_stream_list(run_dir / "stream_list.atmosphere.output")

    partition = find_partition_file(run_dir, ranks, namelist.read_text())
    return run_dir, partition


def remove_stale_outputs(run_dir):
    for pattern in STALE_PATTERNS:
        for path in run_dir.glob(pattern):
            if path.is_file() or path.is_symlink():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)


def run_mpas(run_dir, ranks, mpiexec):
    """Run MPAS in the run directory and capture stdout/stderr in run.out."""
    remove_stale_outputs(run_dir)
    command = [mpiexec, "-n", str(ranks), "./atmosphere_model"]
    with (run_dir / "run.out").open("w") as log:
        result = subprocess.run(
            command,
            cwd=run_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"MPAS failed in {run_dir}; see {run_dir / 'run.out'}")


def read_without_time(dataset, name):
    if name not in dataset.variables:
        raise KeyError(f"output.nc does not contain required variable {name}")

    variable = dataset.variables[name]
    data = variable[:]
    dims = list(variable.dimensions)
    if dims and dims[0] == "Time":
        data = data[0]
        dims = dims[1:]
    return np.asarray(data), dims


def require_level_cell(dataset, name, level_dim):
    data, dims = read_without_time(dataset, name)
    if dims == [level_dim, "nCells"]:
        return data
    if dims == ["nCells", level_dim]:
        return data.T
    raise ValueError(f"{name} has unsupported dimensions {dims}")


def require_cell(dataset, name):
    data, dims = read_without_time(dataset, name)
    if dims != ["nCells"]:
        raise ValueError(f"{name} has unsupported dimensions {dims}")
    return data


def scalar_diagnostic(dataset, name, default=math.nan):
    if name not in dataset.variables:
        return default
    data, _ = read_without_time(dataset, name)
    array = np.asarray(data)
    if array.shape == ():
        return array.item()
    return array.flat[0].item()


def infer_periods(x, y, zgrid):
    x_period = float(np.max(x) - np.min(x))
    y_period = float(np.max(y) - np.min(y))
    z_top = float(np.max(zgrid))
    return (
        x_period if x_period > 0.0 else 1.0,
        y_period if y_period > 0.0 else 1.0,
        z_top if z_top > 0.0 else 1.0,
    )


def compute_errors(output_nc):
    """Compute relative L2 and absolute Linf errors from one MPAS output file."""
    with nc.Dataset(output_nc) as dataset:
        phi = require_level_cell(dataset, "phi", "nVertLevels")
        x = require_cell(dataset, "xCell")
        y = require_cell(dataset, "yCell")
        area = require_cell(dataset, "areaCell")
        zgrid = require_level_cell(dataset, "zgrid", "nVertLevelsP1")

        for field in ("rho_charge", "E_normal", "E_vector"):
            if field not in dataset.variables:
                raise KeyError(f"output.nc does not contain required variable {field}")

        zmid = 0.5 * (zgrid[:-1, :] + zgrid[1:, :])
        dz = zgrid[1:, :] - zgrid[:-1, :]
        x_period, y_period, z_top = infer_periods(x, y, zgrid)

        exact = phi_exact(
            x[None, :],
            y[None, :],
            zmid,
            x_period,
            y_period,
            z_top,
        )
        volume = area[None, :] * dz

        numerator = float(np.sum(volume * (phi - exact) ** 2))
        denominator = float(np.sum(volume * exact**2))
        l2 = math.sqrt(numerator / denominator) if denominator > 0.0 else math.sqrt(numerator)
        linf = float(np.max(np.abs(phi - exact)))

        return {
            "l2": l2,
            "linf": linf,
            "cg_iter_count": int(scalar_diagnostic(dataset, "cg_iter_count", 0)),
            "cg_residual_initial": float(scalar_diagnostic(dataset, "cg_residual_initial")),
            "cg_residual_final": float(scalar_diagnostic(dataset, "cg_residual_final")),
            "nCells": int(phi.shape[1]),
            "nVertLevels": int(phi.shape[0]),
            "x_period": x_period,
            "y_period": y_period,
            "z_top": z_top,
        }


def write_csv(rows, results_dir):
    csv_path = results_dir / "tier_A1_convergence.csv"
    fieldnames = (
        "mesh",
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


def convergence_slope(rows, error_key):
    usable = [
        row
        for row in rows
        if math.isfinite(row["h_m"]) and row[error_key] > 0.0 and math.isfinite(row[error_key])
    ]
    if len(usable) < 2:
        return math.nan
    usable = sorted(usable, key=lambda row: row["h_m"])
    finest = usable[:2]
    return float(
        np.polyfit(
            np.log([row["h_m"] for row in finest]),
            np.log([row[error_key] for row in finest]),
            1,
        )[0]
    )


def write_plot(rows, results_dir, l2_slope, linf_slope):
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
    ax.set_title("Tier A.1 Cartesian MMS convergence")
    ax.legend()
    ax.grid(True, which="both", linestyle=":", linewidth=0.5)

    plot_path = results_dir / "tier_A1_convergence.png"
    fig.savefig(plot_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    return plot_path


def build_rows(args):
    run_root = args.run_root.expanduser().resolve()
    model = args.model.expanduser().resolve()
    results_dir = run_root / "results"
    rows = []

    for mesh_name in args.mesh_list:
        run_dir, partition = prepare_run_dir(run_root, mesh_name, args.ranks, model)
        print(f"{mesh_name}: using {partition.name} with {args.ranks} ranks")

        if args.prepare_only:
            continue

        if not args.analysis_only:
            run_mpas(run_dir, args.ranks, args.mpiexec)

        output = run_dir / "output.nc"
        if not output.exists():
            raise FileNotFoundError(f"{output} was not created")

        result = compute_errors(output)
        row = {
            "mesh": mesh_name,
            "h_m": mesh_spacing_m(mesh_name),
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
        default=1.0e-8,
        help="Maximum accepted final PCG residual.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    rows, results_dir = build_rows(args)
    if args.prepare_only:
        print(f"Prepared {len(args.mesh_list)} run directories under {args.run_root.expanduser()}")
        return 0

    csv_path = write_csv(rows, results_dir)
    l2_slope = convergence_slope(rows, "l2")
    linf_slope = convergence_slope(rows, "linf")
    plot_path = write_plot(rows, results_dir, l2_slope, linf_slope)

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
