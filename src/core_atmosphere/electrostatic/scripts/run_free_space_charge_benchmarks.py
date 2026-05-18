#!/usr/bin/env python3
"""
Free-space Gaussian charge benchmark runner for MPAS electrostatics.

This task-level runner prepares the zero-duration MPAS solve for smooth
Gaussian free-space charge sources, compares output against analytic fields,
and writes summary metrics and diagnostic plots.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import shutil
import subprocess
import sys

import netCDF4 as nc
import numpy as np


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import free_space_charge_analytics as analytic  # noqa: E402
import run_tier_A1_cartesian_mms as tier_a1  # noqa: E402
from run_tripole_supercell import ensure_output_stream_list, seed_run_dir_from_template  # noqa: E402


SOURCES = analytic.SOURCES

ELECTROSTATIC_CONFIG = (
    ("config_electrostatic_enable", ".true."),
    ("config_electrostatic_solve_at_init", ".true."),
    ("config_electrostatic_source", "'gaussian_dipole_y'"),
    ("config_poisson_preconditioner", "'jacobi'"),
    ("config_poisson_tol", "1.0e-10"),
    ("config_poisson_max_iter", "5000"),
    ("config_electrostatic_bc_ground", "0.0"),
    ("config_electrostatic_source_charge", "20.0"),
    ("config_electrostatic_source_sigma", "2000.0"),
    ("config_electrostatic_dipole_separation", "8000.0"),
)

STALE_PATTERNS = (
    "output.nc",
    "run.out",
    "log.atmosphere.*.out",
    "log.atmosphere.*.err",
)


def configure_namelist_text(
    text,
    *,
    source,
    charge,
    sigma,
    separation,
    poisson_tol,
    poisson_max_iter,
):
    """Return namelist text configured for a zero-duration Gaussian solve."""
    if source not in SOURCES:
        raise ValueError(f"unsupported Gaussian source: {source}")

    text = tier_a1.set_namelist_value(
        text,
        "nhyd_model",
        "config_run_duration",
        "'00_00:00:00'",
    )
    entries = tuple(
        (key, f"'{source}'")
        if key == "config_electrostatic_source"
        else (key, str(poisson_tol))
        if key == "config_poisson_tol"
        else (key, str(poisson_max_iter))
        if key == "config_poisson_max_iter"
        else (key, str(charge))
        if key == "config_electrostatic_source_charge"
        else (key, str(sigma))
        if key == "config_electrostatic_source_sigma"
        else (key, str(separation))
        if key == "config_electrostatic_dipole_separation"
        else (key, value)
        for key, value in ELECTROSTATIC_CONFIG
    )
    text = tier_a1.replace_namelist_block(text, "electrostatic", entries)
    text = re.sub(r"\n{2,}(&electrostatic)", r"\n\1", text)
    return text.rstrip() + "\n"


def configure_namelist(path, **kwargs):
    """Configure namelist.atmosphere in place."""
    path.write_text(configure_namelist_text(path.read_text(), **kwargs))


def prepare_run_dir(
    template_run_dir,
    run_dir,
    model,
    ranks,
    source,
    charge,
    sigma,
    separation,
    poisson_tol,
    poisson_max_iter,
):
    """Create and configure the isolated free-space benchmark run directory."""
    run_dir = run_dir.expanduser().resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    seed_run_dir_from_template(template_run_dir, run_dir)
    tier_a1.symlink_model(model.expanduser(), run_dir)

    namelist = run_dir / "namelist.atmosphere"
    streams = run_dir / "streams.atmosphere"
    if not namelist.exists() or not streams.exists():
        raise FileNotFoundError(
            f"{run_dir} needs namelist.atmosphere and streams.atmosphere; "
            f"seed it from a complete supercell run directory"
        )

    configure_namelist(
        namelist,
        source=source,
        charge=charge,
        sigma=sigma,
        separation=separation,
        poisson_tol=poisson_tol,
        poisson_max_iter=poisson_max_iter,
    )
    tier_a1.ensure_streams_netcdf(streams)
    ensure_output_stream_list(run_dir / "stream_list.atmosphere.output")

    partition = tier_a1.find_partition_file(run_dir, ranks, namelist.read_text())
    return run_dir, partition


def remove_stale_outputs(run_dir):
    """Remove stale output products before rerunning MPAS."""
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


def read_e_vector(dataset):
    """Read E_vector as an array with shape (R3, nVertLevels, nCells)."""
    data, dims = tier_a1.read_without_time(dataset, "E_vector")
    expected = ("R3", "nVertLevels", "nCells")
    if sorted(dims) != sorted(expected):
        raise ValueError(f"E_vector has unsupported dimensions {dims}")
    axes = [dims.index(dim) for dim in expected]
    return np.transpose(data, axes)


def read_fields(output_nc):
    """Read MPAS output fields needed for free-space Gaussian diagnostics."""
    with nc.Dataset(output_nc) as dataset:
        x = tier_a1.require_cell(dataset, "xCell")
        y = tier_a1.require_cell(dataset, "yCell")
        area = tier_a1.require_cell(dataset, "areaCell")
        zgrid = tier_a1.require_level_cell(dataset, "zgrid", "nVertLevelsP1")
        rho = tier_a1.require_level_cell(dataset, "rho_charge", "nVertLevels")
        phi = tier_a1.require_level_cell(dataset, "phi", "nVertLevels")
        e_vector = read_e_vector(dataset)
        residual = tier_a1.scalar_diagnostic(dataset, "cg_residual_final")

    zmid = 0.5 * (zgrid[:-1, :] + zgrid[1:, :])
    volume = area[np.newaxis, :] * np.diff(zgrid, axis=0)
    return {
        "x": x,
        "y": y,
        "zgrid": zgrid,
        "zmid": zmid,
        "area": area,
        "volume": volume,
        "rho": rho,
        "phi": phi,
        "ex": e_vector[0, :, :],
        "ey": e_vector[1, :, :],
        "ez": e_vector[2, :, :],
        "residual": residual,
    }


def analytic_center(fields):
    """Infer the analytic source center from the output domain bounds."""
    return (
        0.5 * (float(np.min(fields["x"])) + float(np.max(fields["x"]))),
        0.5 * (float(np.min(fields["y"])) + float(np.max(fields["y"]))),
        0.5 * (float(np.min(fields["zgrid"])) + float(np.max(fields["zgrid"]))),
    )


def _field_bounds(fields):
    return (
        (float(np.min(fields["x"])), float(np.max(fields["x"]))),
        (float(np.min(fields["y"])), float(np.max(fields["y"]))),
        (float(np.min(fields["zgrid"])), float(np.max(fields["zgrid"]))),
    )


def _broadcast_cell_coordinates(fields):
    x = np.broadcast_to(fields["x"], fields["zmid"].shape)
    y = np.broadcast_to(fields["y"], fields["zmid"].shape)
    return x, y


def _rho_relative_error(fields, source, charge, sigma, separation, center):
    x, y = _broadcast_cell_coordinates(fields)
    candidate = analytic.evaluate_gaussian_source(
        source,
        x,
        y,
        fields["zmid"],
        center=center,
        charge=charge,
        sigma=sigma,
        separation=separation,
    )
    diff = fields["rho"] - candidate.rho
    numerator = float(np.sum(fields["volume"] * diff * diff))
    denominator = float(np.sum(fields["volume"] * candidate.rho * candidate.rho))
    return math.sqrt(numerator / denominator) if denominator > 0.0 else math.inf


def _peak_inferred_center(fields, source, separation):
    x, y = _broadcast_cell_coordinates(fields)
    index = np.unravel_index(int(np.argmax(np.abs(fields["rho"]))), fields["rho"].shape)
    center = [float(x[index]), float(y[index]), float(fields["zmid"][index])]
    if source == "gaussian_dipole_y":
        center[1] += math.copysign(0.5 * separation, fields["rho"][index])
    elif source == "gaussian_dipole_z":
        center[2] += math.copysign(0.5 * separation, fields["rho"][index])
    return tuple(center)


def infer_source_center(fields, source, charge, sigma, separation):
    """Use the domain center unless the emitted source clearly says otherwise."""
    center = analytic_center(fields)
    if _rho_relative_error(fields, source, charge, sigma, separation, center) < 1.0e-6:
        return center
    return _peak_inferred_center(fields, source, separation)


def analyze_output(
    output_nc,
    summary_path,
    plot_path,
    *,
    source,
    charge,
    sigma,
    separation,
    boundary_margin,
    core_radius,
    e_relative_floor,
):
    """Analyze MPAS output against the matching free-space analytic source."""
    fields = read_fields(pathlib.Path(output_nc).expanduser())
    x, y = _broadcast_cell_coordinates(fields)
    center = infer_source_center(fields, source, charge, sigma, separation)
    exact = analytic.evaluate_gaussian_source(
        source,
        x,
        y,
        fields["zmid"],
        center=center,
        charge=charge,
        sigma=sigma,
        separation=separation,
    )
    bounds = _field_bounds(fields)
    mask = analytic.interior_comparison_mask(
        x,
        y,
        fields["zmid"],
        bounds=bounds,
        centers=exact.centers,
        sigma=sigma,
        boundary_margin=boundary_margin,
        core_radius=core_radius,
    )
    phi_aligned, phi_offset = analytic.align_potential_gauge(
        fields["phi"],
        exact.phi,
        fields["volume"],
        mask,
    )
    phi_norms = analytic.weighted_error_norms(
        phi_aligned,
        exact.phi,
        fields["volume"],
        mask,
    )

    e_mpas = np.sqrt(fields["ex"] ** 2 + fields["ey"] ** 2 + fields["ez"] ** 2)
    e_exact = np.sqrt(exact.ex**2 + exact.ey**2 + exact.ez**2)
    e_norms = analytic.weighted_error_norms(
        e_mpas,
        e_exact,
        fields["volume"],
        mask,
        relative_floor=e_relative_floor,
    )

    result = {
        "source": source,
        "charge": float(charge),
        "sigma": float(sigma),
        "separation": float(separation),
        "boundary_margin": float(boundary_margin),
        "core_radius": None if core_radius is None else float(core_radius),
        "e_relative_floor": float(e_relative_floor),
        "center": [float(value) for value in center],
        "bounds": [[float(lo), float(hi)] for lo, hi in bounds],
        "mask_count": int(np.count_nonzero(mask)),
        "cell_count": int(mask.size),
        "phi_gauge_offset": float(phi_offset),
        "phi_l2_absolute": phi_norms["l2_absolute"],
        "phi_l2_relative": phi_norms["l2_relative"],
        "phi_linf_absolute": phi_norms["linf_absolute"],
        "e_l2_absolute": e_norms["l2_absolute"],
        "e_l2_relative": e_norms["l2_relative"],
        "e_linf_absolute": e_norms["linf_absolute"],
        "cg_residual_final": float(fields["residual"]),
    }

    summary_path = pathlib.Path(summary_path).expanduser()
    plot_path = pathlib.Path(plot_path).expanduser()
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    plot_diagnostics(fields, exact, phi_aligned, mask, plot_path, source)
    return result


def _flatten_xy(fields):
    x, y = _broadcast_cell_coordinates(fields)
    return x.ravel(), y.ravel()


def _normalized_quiver_components(ex, ey):
    scale = np.sqrt(ex * ex + ey * ey)
    safe = scale > 0.0
    u = np.zeros_like(ex, dtype=float)
    v = np.zeros_like(ey, dtype=float)
    u[safe] = ex[safe] / scale[safe]
    v[safe] = ey[safe] / scale[safe]
    return u, v


def plot_diagnostics(fields, exact, phi_aligned, mask, plot_path, source):
    """Write a 2x4 scatter diagnostic figure for the free-space benchmark."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x, y = _flatten_xy(fields)
    e_mpas = np.sqrt(fields["ex"] ** 2 + fields["ey"] ** 2 + fields["ez"] ** 2)
    e_exact = np.sqrt(exact.ex**2 + exact.ey**2 + exact.ez**2)
    panels = (
        ("MPAS rho", fields["rho"]),
        ("MPAS Phi", phi_aligned),
        ("MPAS |E|", e_mpas),
        ("Phi error", phi_aligned - exact.phi),
        ("Analytic rho", exact.rho),
        ("Analytic Phi", exact.phi),
        ("Analytic |E|", e_exact),
        ("Comparison mask", mask.astype(float)),
    )

    fig, axes = plt.subplots(2, 4, figsize=(16, 8), constrained_layout=True)
    for ax, (title, values) in zip(axes.flat, panels):
        scatter = ax.scatter(x, y, c=np.asarray(values).ravel(), s=18, cmap="viridis")
        ax.set_title(title)
        ax.set_xlabel("xCell (m)")
        ax.set_ylabel("yCell (m)")
        ax.set_aspect("equal", adjustable="box")
        fig.colorbar(scatter, ax=ax, shrink=0.78)

    mpas_u, mpas_v = _normalized_quiver_components(fields["ex"], fields["ey"])
    exact_u, exact_v = _normalized_quiver_components(exact.ex, exact.ey)
    axes.flat[1].quiver(x, y, mpas_u.ravel(), mpas_v.ravel(), color="black", scale=24.0, width=0.004)
    axes.flat[5].quiver(x, y, exact_u.ravel(), exact_v.ravel(), color="black", scale=24.0, width=0.004)
    fig.suptitle(f"Free-space charge diagnostics: {source}")
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        choices=SOURCES,
        default="gaussian_dipole_y",
        help="Gaussian charge source mode to write to config_electrostatic_source.",
    )
    parser.add_argument(
        "--template-run-dir",
        type=pathlib.Path,
        default=pathlib.Path("~/Data/MPAS/supercell"),
        help="Existing supercell run directory to seed inputs from.",
    )
    parser.add_argument(
        "--run-dir",
        type=pathlib.Path,
        default=None,
        help="Isolated run directory to create or reuse.",
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
        "--charge",
        type=float,
        default=20.0,
        help="Gaussian lobe charge written to namelist.atmosphere.",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=2000.0,
        help="Gaussian width in meters written to namelist.atmosphere.",
    )
    parser.add_argument(
        "--separation",
        type=float,
        default=8000.0,
        help="Dipole lobe separation in meters written to namelist.atmosphere.",
    )
    parser.add_argument(
        "--boundary-margin",
        type=float,
        default=12000.0,
        help="Interior comparison margin in meters for the follow-on analysis step.",
    )
    parser.add_argument(
        "--core-radius",
        type=float,
        default=None,
        help="Optional source-core exclusion radius in meters for follow-on analysis.",
    )
    parser.add_argument(
        "--e-relative-floor",
        type=float,
        default=0.0,
        help="Relative-error floor for electric-field analysis.",
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
        "--summary",
        type=pathlib.Path,
        default=None,
        help="Output JSON summary path for the follow-on analysis step.",
    )
    parser.add_argument(
        "--plot",
        type=pathlib.Path,
        default=None,
        help="Output PNG diagnostics path for the follow-on analysis step.",
    )
    parser.add_argument(
        "--analysis-only",
        action="store_true",
        help="Skip MPAS execution and analyze an existing output.nc.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare the run directory but do not run or analyze output.",
    )
    args = parser.parse_args(argv)

    root = pathlib.Path("~/Data/MPAS/poisson_free_space_charge")
    if args.run_dir is None:
        args.run_dir = root / args.source / "run"
    if args.summary is None:
        args.summary = root / "results" / f"{args.source}_summary.json"
    if args.plot is None:
        args.plot = root / "results" / f"{args.source}_diagnostics.png"
    return args


def main(argv=None):
    args = parse_args(argv)
    run_dir, partition = prepare_run_dir(
        args.template_run_dir,
        args.run_dir,
        args.model,
        args.ranks,
        args.source,
        args.charge,
        args.sigma,
        args.separation,
        args.poisson_tol,
        args.poisson_max_iter,
    )
    print(f"Prepared run directory: {run_dir}")
    print(f"Using partition file: {partition}")
    if args.prepare_only:
        return 0

    if not args.analysis_only:
        run_mpas(run_dir, args.ranks, args.mpiexec)

    output_nc = run_dir / "output.nc"
    result = analyze_output(
        output_nc,
        args.summary,
        args.plot,
        source=args.source,
        charge=args.charge,
        sigma=args.sigma,
        separation=args.separation,
        boundary_margin=args.boundary_margin,
        core_radius=args.core_radius,
        e_relative_floor=args.e_relative_floor,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Wrote summary: {pathlib.Path(args.summary).expanduser()}")
    print(f"Wrote plot: {pathlib.Path(args.plot).expanduser()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
