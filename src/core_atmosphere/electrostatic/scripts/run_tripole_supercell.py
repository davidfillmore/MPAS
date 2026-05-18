#!/usr/bin/env python3
"""
Synthetic thundercloud tripole runner for MPAS electrostatics.

The script prepares an isolated zero-duration supercell-mesh run, enables the
diagnostic electrostatic Poisson solve with config_electrostatic_source =
'tripole', and plots rho_charge, phi, and |E| from output.nc.
"""

from __future__ import annotations

import argparse
import math
import os
import pathlib
import re
import shutil
import sys

import netCDF4 as nc
import numpy as np


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_tier_A1_cartesian_mms as tier_a1  # noqa: E402


ELECTROSTATIC_CONFIG = (
    ("config_electrostatic_enable", ".true."),
    ("config_electrostatic_solve_at_init", ".true."),
    ("config_electrostatic_source", "'tripole'"),
    ("config_poisson_preconditioner", "'jacobi'"),
    ("config_poisson_tol", "1.0e-10"),
    ("config_poisson_max_iter", "5000"),
    ("config_electrostatic_bc_ground", "0.0"),
)

REQUIRED_OUTPUT_FIELDS = (
    "xCell",
    "yCell",
    "zgrid",
    "areaCell",
    "rho_charge",
    "phi",
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


def configure_namelist_text(text, poisson_tol=1.0e-10, poisson_max_iter=5000):
    """Return namelist text configured for a zero-duration tripole solve."""
    text = tier_a1.set_namelist_value(
        text,
        "nhyd_model",
        "config_run_duration",
        "'00_00:00:00'",
    )
    entries = tuple(
        (key, f"{poisson_tol:.1e}")
        if key == "config_poisson_tol"
        else (key, str(poisson_max_iter))
        if key == "config_poisson_max_iter"
        else (key, value)
        for key, value in ELECTROSTATIC_CONFIG
    )
    text = tier_a1.replace_namelist_block(text, "electrostatic", entries)
    text = re.sub(r"\n{2,}(&electrostatic)", r"\n\1", text)
    return text.rstrip() + "\n"


def configure_namelist(path, poisson_tol, poisson_max_iter):
    """Configure namelist.atmosphere in place."""
    path.write_text(
        configure_namelist_text(
            path.read_text(),
            poisson_tol=poisson_tol,
            poisson_max_iter=poisson_max_iter,
        )
    )


def ensure_output_stream_list(stream_list_path):
    """Ensure the output stream list contains fields needed for tripole plots."""
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


def should_skip_template_item(path):
    """Return True for generated files that should not seed a fresh run."""
    name = path.name
    return (
        name == "atmosphere_model"
        or name == "init_atmosphere_model"
        or name == "output.nc"
        or name == "run.out"
        or name.startswith("run_")
        or name.startswith("log.")
        or name == "restart_timestamp"
    )


def should_link_template_item(path):
    """Return True when a template item should be symlinked instead of copied."""
    name = path.name
    return path.is_dir() or path.suffix == ".nc" or "graph.info" in name


def seed_run_dir_from_template(template_run_dir, run_dir):
    """Copy small template inputs and link large mesh/output inputs."""
    template_run_dir = template_run_dir.expanduser().resolve()
    if not template_run_dir.exists():
        raise FileNotFoundError(f"template run directory does not exist: {template_run_dir}")

    run_dir.mkdir(parents=True, exist_ok=True)
    for source in template_run_dir.iterdir():
        if should_skip_template_item(source):
            continue

        target = run_dir / source.name
        if target.exists() or target.is_symlink():
            continue

        if should_link_template_item(source):
            os.symlink(source.resolve(), target)
        elif source.is_file():
            shutil.copy2(source, target)


def prepare_run_dir(
    template_run_dir,
    run_dir,
    model,
    ranks,
    poisson_tol=1.0e-10,
    poisson_max_iter=5000,
):
    """Create and configure the isolated tripole run directory."""
    run_dir = run_dir.expanduser().resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir.parent / "results").mkdir(parents=True, exist_ok=True)

    seed_run_dir_from_template(template_run_dir, run_dir)
    tier_a1.symlink_model(model.expanduser(), run_dir)

    namelist = run_dir / "namelist.atmosphere"
    streams = run_dir / "streams.atmosphere"
    if not namelist.exists() or not streams.exists():
        raise FileNotFoundError(
            f"{run_dir} needs namelist.atmosphere and streams.atmosphere; "
            f"seed it from a complete supercell run directory"
        )

    configure_namelist(namelist, poisson_tol, poisson_max_iter)
    tier_a1.ensure_streams_netcdf(streams)
    ensure_output_stream_list(run_dir / "stream_list.atmosphere.output")

    partition = tier_a1.find_partition_file(run_dir, ranks, namelist.read_text())
    return run_dir, partition


def read_e_vector(dataset):
    """Read E_vector as (R3,nVertLevels,nCells) independent of stored order."""
    data, dims = tier_a1.read_without_time(dataset, "E_vector")
    required = ("R3", "nVertLevels", "nCells")
    if sorted(dims) != sorted(required):
        raise ValueError(f"E_vector has unsupported dimensions {dims}")
    order = [dims.index(dim) for dim in required]
    return np.transpose(data, order)


def centerline_indices(x, y):
    """Return cell indices along a vertical x-slice through the domain center."""
    x0 = 0.5 * (float(np.min(x)) + float(np.max(x)))
    x_slice = x[int(np.argmin(np.abs(x - x0)))]
    indices = np.where(np.isclose(x, x_slice))[0]
    if len(indices) == 0:
        indices = np.array([int(np.argmin(np.abs(x - x0)))])
    return indices[np.argsort(y[indices])]


def read_plot_fields(output_nc):
    """Read tripole plotting fields from an MPAS output.nc file."""
    with nc.Dataset(output_nc) as dataset:
        x = tier_a1.require_cell(dataset, "xCell")
        y = tier_a1.require_cell(dataset, "yCell")
        zgrid = tier_a1.require_level_cell(dataset, "zgrid", "nVertLevelsP1")
        rho_charge = tier_a1.require_level_cell(dataset, "rho_charge", "nVertLevels")
        phi = tier_a1.require_level_cell(dataset, "phi", "nVertLevels")
        e_vector = read_e_vector(dataset)

    zmid = 0.5 * (zgrid[:-1, :] + zgrid[1:, :])
    e_mag = np.sqrt(np.sum(e_vector * e_vector, axis=0))
    mid_level = int(rho_charge.shape[0] // 2)
    indices = centerline_indices(x, y)

    return {
        "xCell": x,
        "yCell": y,
        "zgrid": zgrid,
        "zMid": zmid,
        "rho_charge": rho_charge,
        "phi": phi,
        "E_mag": e_mag,
        "mid_level": mid_level,
        "cross_section_cell_indices": indices,
    }


def centered_norm(values):
    """Return a diverging normalization centered at zero, when possible."""
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        return None
    bound = float(np.max(np.abs(finite)))
    if bound <= 0.0:
        return None

    from matplotlib.colors import TwoSlopeNorm

    return TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound)


def scatter_with_colorbar(fig, ax, x, y, values, title, xlabel, ylabel, cmap, norm=None):
    """Draw one scatter panel and attach a compact colorbar."""
    marker_size = max(6.0, min(22.0, 50000.0 / max(1, values.size)))
    artist = ax.scatter(x, y, c=values, s=marker_size, cmap=cmap, norm=norm, linewidths=0.0)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, linestyle=":", linewidth=0.4)
    fig.colorbar(artist, ax=ax, fraction=0.046, pad=0.04)


def plot_tripole_output(output_nc, plot_path):
    """Plot horizontal and vertical tripole sections from output.nc."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fields = read_plot_fields(output_nc)
    x_km = fields["xCell"] / 1000.0
    y_km = fields["yCell"] / 1000.0
    z_km = fields["zMid"] / 1000.0
    mid = fields["mid_level"]
    section = fields["cross_section_cell_indices"]

    panels = (
        ("rho_charge", "rho_charge", "RdBu_r", centered_norm(fields["rho_charge"])),
        ("phi", "phi", "RdBu_r", centered_norm(fields["phi"])),
        ("E_mag", "|E|", "viridis", None),
    )

    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.0), constrained_layout=True)
    for col, (field_name, label, cmap, norm) in enumerate(panels):
        values = fields[field_name]
        scatter_with_colorbar(
            fig,
            axes[0, col],
            x_km,
            y_km,
            values[mid, :],
            f"{label}, z={float(np.nanmean(z_km[mid, :])):.1f} km",
            "x (km)",
            "y (km)",
            cmap,
            norm,
        )

        section_values = values[:, section]
        section_y = np.broadcast_to(y_km[section][None, :], section_values.shape)
        scatter_with_colorbar(
            fig,
            axes[1, col],
            section_y.ravel(),
            z_km[:, section].ravel(),
            section_values.ravel(),
            f"{label}, vertical centerline",
            "y (km)",
            "z (km)",
            cmap,
            norm,
        )

    plot_path = plot_path.expanduser()
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    return plot_path


def scalar_diagnostic(output_nc, name, default=math.nan):
    """Read one scalar diagnostic from output.nc if present."""
    with nc.Dataset(output_nc) as dataset:
        return tier_a1.scalar_diagnostic(dataset, name, default=default)


def analyze_output(run_dir, plot_path, residual_tol):
    """Plot output.nc and enforce a final residual gate when available."""
    output = run_dir / "output.nc"
    if not output.exists():
        raise FileNotFoundError(f"{output} was not created")

    plot = plot_tripole_output(output, plot_path)
    final_residual = float(scalar_diagnostic(output, "cg_residual_final"))
    if math.isfinite(final_residual) and final_residual > residual_tol:
        raise RuntimeError(
            f"final CG residual {final_residual:.6e} exceeds {residual_tol:.6e}"
        )

    print(f"Wrote {plot}")
    if math.isfinite(final_residual):
        print(f"Final CG residual: {final_residual:.6e}")
    return plot


def remove_stale_outputs(run_dir):
    """Remove stale output products before rerunning MPAS."""
    for pattern in STALE_PATTERNS:
        for path in run_dir.glob(pattern):
            if path.is_file() or path.is_symlink():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--template-run-dir",
        type=pathlib.Path,
        default=pathlib.Path("~/Data/MPAS/supercell"),
        help="Existing supercell run directory to seed inputs from.",
    )
    parser.add_argument(
        "--run-dir",
        type=pathlib.Path,
        default=pathlib.Path("~/Data/MPAS/poisson_tripole_supercell/run"),
        help="Isolated tripole run directory to create or reuse.",
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
        "--plot",
        type=pathlib.Path,
        default=None,
        help="Output PNG path. Defaults to run-dir/../results/tripole_supercell.png.",
    )
    parser.add_argument(
        "--analysis-only",
        action="store_true",
        help="Skip MPAS execution and plot an existing output.nc.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare the run directory but do not run or plot.",
    )
    parser.add_argument(
        "--poisson-tol",
        type=float,
        default=1.0e-10,
        help="PCG tolerance written to namelist.atmosphere.",
    )
    parser.add_argument(
        "--poisson-max-iter",
        type=int,
        default=5000,
        help="PCG iteration cap written to namelist.atmosphere.",
    )
    parser.add_argument(
        "--residual-tol",
        type=float,
        default=1.0e-8,
        help="Maximum accepted final PCG residual when the diagnostic is present.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    run_dir, partition = prepare_run_dir(
        args.template_run_dir,
        args.run_dir,
        args.model,
        args.ranks,
        poisson_tol=args.poisson_tol,
        poisson_max_iter=args.poisson_max_iter,
    )

    print(f"Prepared {run_dir}")
    print(f"Using {partition.name} with {args.ranks} ranks")
    if args.prepare_only:
        return 0

    if not args.analysis_only:
        remove_stale_outputs(run_dir)
        tier_a1.run_mpas(run_dir, args.ranks, args.mpiexec)

    plot_path = args.plot
    if plot_path is None:
        plot_path = run_dir.parent / "results" / "tripole_supercell.png"
    analyze_output(run_dir, plot_path, args.residual_tol)
    return 0


if __name__ == "__main__":
    sys.exit(main())
