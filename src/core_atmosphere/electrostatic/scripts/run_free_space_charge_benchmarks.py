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
KM_PER_M = 1.0e-3
SCIENTIFIC_POWER_LIMITS = (-2, 3)
MAX_QUIVER_VECTORS = 441

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


def comparison_mask(
    x,
    y,
    z,
    *,
    bounds,
    centers,
    sigma,
    horizontal_boundary_margin,
    vertical_boundary_margin,
    core_radius=None,
):
    """Build the interior comparison mask with separate horizontal/vertical margins."""
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = bounds
    mask = (
        (x >= xmin + horizontal_boundary_margin)
        & (x <= xmax - horizontal_boundary_margin)
        & (y >= ymin + horizontal_boundary_margin)
        & (y <= ymax - horizontal_boundary_margin)
        & (z >= zmin + vertical_boundary_margin)
        & (z <= zmax - vertical_boundary_margin)
    )
    exclusion_radius = core_radius if core_radius is not None else 3.0 * sigma
    for cx, cy, cz in centers:
        radius = np.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2)
        mask = mask & (radius >= exclusion_radius)
    return mask


def representative_level(fields, center):
    """Return the vertical level nearest the analytic center height."""
    level_z = np.mean(fields["zmid"], axis=1)
    level = int(np.argmin(np.abs(level_z - center[2])))
    return level, float(level_z[level])


def vector_error_norms(actual_components, exact_components, weights, mask, *, relative_floor=0.0):
    """Compute weighted norms for vector component errors."""
    if relative_floor < 0.0:
        raise ValueError("relative_floor must be nonnegative")
    mask = np.asarray(mask, dtype=bool)
    weights = np.asarray(weights, dtype=float)
    if not np.any(mask):
        raise ValueError("comparison mask is empty")
    masked_weights = weights[mask]
    if not np.all(np.isfinite(masked_weights)) or np.any(masked_weights < 0.0):
        raise ValueError("masked weights must be finite and nonnegative")
    masked_weight_sum = float(np.sum(masked_weights))
    if masked_weight_sum == 0.0:
        raise ValueError("masked weight sum must be nonzero")

    error_squared = np.zeros_like(weights, dtype=float)
    exact_squared = np.zeros_like(weights, dtype=float)
    for actual, exact in zip(actual_components, exact_components):
        actual = np.asarray(actual, dtype=float)
        exact = np.asarray(exact, dtype=float)
        error_squared = error_squared + (actual - exact) ** 2
        exact_squared = exact_squared + exact**2

    numerator = float(np.sum(masked_weights * error_squared[mask]))
    denominator = float(np.sum(masked_weights * exact_squared[mask]))
    if relative_floor > 0.0:
        denominator = max(denominator, masked_weight_sum * relative_floor**2)
    return {
        "l2_absolute": float(np.sqrt(numerator)),
        "l2_relative": float(np.sqrt(numerator / denominator)) if denominator > 0.0 else float("nan"),
        "linf_absolute": float(np.max(np.sqrt(error_squared[mask]))),
        "n": int(np.count_nonzero(mask)),
    }


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
    vertical_boundary_margin=0.0,
    core_radius,
    e_relative_floor,
):
    """Analyze MPAS output against the matching free-space analytic source."""
    fields = read_fields(pathlib.Path(output_nc).expanduser())
    x, y = _broadcast_cell_coordinates(fields)
    center = analytic_center(fields)
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
    mask = comparison_mask(
        x,
        y,
        fields["zmid"],
        bounds=bounds,
        centers=exact.centers,
        sigma=sigma,
        horizontal_boundary_margin=boundary_margin,
        vertical_boundary_margin=vertical_boundary_margin,
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

    e_norms = vector_error_norms(
        (fields["ex"], fields["ey"], fields["ez"]),
        (exact.ex, exact.ey, exact.ez),
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
        "vertical_boundary_margin": float(vertical_boundary_margin),
        "core_radius": None if core_radius is None else float(core_radius),
        "e_relative_floor": float(e_relative_floor),
        "center": [float(value) for value in center],
        "bounds": [[float(lo), float(hi)] for lo, hi in bounds],
        "mask_count": int(np.count_nonzero(mask)),
        "comparison_points": int(np.count_nonzero(mask)),
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
    result["plot_level"], result["plot_z"] = representative_level(fields, center)

    summary_path = pathlib.Path(summary_path).expanduser()
    plot_path = pathlib.Path(plot_path).expanduser()
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    plot_diagnostics(
        fields,
        exact,
        phi_aligned,
        mask,
        plot_path,
        source,
        result["plot_level"],
        result["plot_z"],
    )
    return result


def _normalized_quiver_components(ex, ey):
    scale = np.sqrt(ex * ex + ey * ey)
    safe = scale > 0.0
    u = np.zeros_like(ex, dtype=float)
    v = np.zeros_like(ey, dtype=float)
    u[safe] = ex[safe] / scale[safe]
    v[safe] = ey[safe] / scale[safe]
    return u, v


def _quiver_sample_indices(x, y, u, v, max_vectors=MAX_QUIVER_VECTORS):
    """Return spatially distributed vector indices for readable quiver overlays."""
    if max_vectors < 1:
        raise ValueError("max_vectors must be positive")

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    if not (x.shape == y.shape == u.shape == v.shape):
        raise ValueError("x, y, u, and v must have matching shapes")

    magnitude = np.sqrt(u * u + v * v)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(u) & np.isfinite(v) & (magnitude > 0.0)
    valid_indices = np.flatnonzero(valid)
    if valid_indices.size <= max_vectors:
        return valid_indices

    bins = max(1, int(np.floor(np.sqrt(max_vectors))))
    xv = x[valid_indices]
    yv = y[valid_indices]
    xmin = float(np.min(xv))
    xmax = float(np.max(xv))
    ymin = float(np.min(yv))
    ymax = float(np.max(yv))
    x_span = xmax - xmin
    y_span = ymax - ymin
    if x_span == 0.0 or y_span == 0.0:
        stride = int(np.ceil(valid_indices.size / max_vectors))
        return valid_indices[::stride][:max_vectors]

    xbin = np.floor((xv - xmin) / x_span * bins).astype(int)
    ybin = np.floor((yv - ymin) / y_span * bins).astype(int)
    xbin = np.clip(xbin, 0, bins - 1)
    ybin = np.clip(ybin, 0, bins - 1)
    bin_id = xbin * bins + ybin

    selected = []
    for current_bin in np.unique(bin_id):
        candidates = np.flatnonzero(bin_id == current_bin)
        bx = current_bin // bins
        by = current_bin % bins
        cx = xmin + (float(bx) + 0.5) * x_span / bins
        cy = ymin + (float(by) + 0.5) * y_span / bins
        distance2 = (xv[candidates] - cx) ** 2 + (yv[candidates] - cy) ** 2
        selected.append(valid_indices[candidates[int(np.argmin(distance2))]])
    return np.asarray(selected, dtype=int)


def _scientific_scalar_formatter():
    """Return the plot formatter used for compact scientific colorbar labels."""
    import matplotlib.ticker as mticker

    formatter = mticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits(SCIENTIFIC_POWER_LIMITS)
    formatter.set_useOffset(False)
    return formatter


def _shared_scalar_norm(*arrays, signed=False):
    """Return a shared color normalization for one comparable field group."""
    import matplotlib.colors as mcolors

    finite_parts = []
    for values in arrays:
        values = np.asarray(values, dtype=float)
        finite_parts.append(values[np.isfinite(values)])
    finite_parts = [values for values in finite_parts if values.size > 0]
    if not finite_parts:
        return None
    finite = np.concatenate(finite_parts)
    if finite.size == 0:
        return None
    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    if vmin == vmax:
        return None
    if signed:
        limit = float(np.max(np.abs(finite)))
        if limit == 0.0:
            return None
        return mcolors.TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    return mcolors.Normalize(vmin=vmin, vmax=vmax)


def _scatter_panel(fig, ax, x, y, values, *, title, cmap, norm=None):
    scatter = ax.scatter(x, y, c=np.asarray(values), s=18, cmap=cmap, norm=norm)
    ax.set_title(title)
    ax.set_xlabel("xCell (km)")
    ax.set_ylabel("yCell (km)")
    ax.set_aspect("equal", adjustable="box")
    fig.colorbar(scatter, ax=ax, shrink=0.78, format=_scientific_scalar_formatter())
    return scatter


def _overlay_quiver(ax, x, y, u, v):
    indices = _quiver_sample_indices(x, y, u, v)
    if indices.size == 0:
        return None
    return ax.quiver(
        x[indices],
        y[indices],
        u[indices],
        v[indices],
        color="black",
        alpha=0.55,
        angles="xy",
        scale_units="width",
        scale=32.0,
        width=0.002,
        pivot="middle",
    )


def plot_diagnostics(fields, exact, phi_aligned, mask, plot_path, source, plot_level, plot_z):
    """Write a 2x4 scatter diagnostic figure for the free-space benchmark."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = KM_PER_M * fields["x"]
    y = KM_PER_M * fields["y"]
    e_mpas = np.sqrt(fields["ex"] ** 2 + fields["ey"] ** 2 + fields["ez"] ** 2)
    e_exact = np.sqrt(exact.ex**2 + exact.ey**2 + exact.ez**2)
    mpas_rho = fields["rho"][plot_level, :]
    analytic_rho = exact.rho[plot_level, :]
    mpas_phi = phi_aligned[plot_level, :]
    analytic_phi = exact.phi[plot_level, :]
    mpas_e = e_mpas[plot_level, :]
    analytic_e = e_exact[plot_level, :]
    phi_error = (phi_aligned - exact.phi)[plot_level, :]
    mask_values = mask[plot_level, :].astype(float)
    rho_norm = _shared_scalar_norm(mpas_rho, analytic_rho, signed=True)
    phi_norm = _shared_scalar_norm(mpas_phi, analytic_phi, signed=True)
    e_norm = _shared_scalar_norm(mpas_e, analytic_e)
    panels = (
        ("MPAS rho", mpas_rho, "RdBu_r", rho_norm),
        ("MPAS Phi", mpas_phi, "RdBu_r", phi_norm),
        ("MPAS |E|", mpas_e, "viridis", e_norm),
        ("Phi error", phi_error, "RdBu_r", _shared_scalar_norm(phi_error, signed=True)),
        ("Analytic rho", analytic_rho, "RdBu_r", rho_norm),
        ("Analytic Phi", analytic_phi, "RdBu_r", phi_norm),
        ("Analytic |E|", analytic_e, "viridis", e_norm),
        ("Comparison mask", mask_values, "viridis", _shared_scalar_norm(mask_values)),
    )

    fig, axes = plt.subplots(2, 4, figsize=(16, 8), constrained_layout=True)
    for ax, (title, values, cmap, norm) in zip(axes.flat, panels):
        _scatter_panel(
            fig,
            ax,
            x,
            y,
            values,
            title=title,
            cmap=cmap,
            norm=norm,
        )

    mpas_u, mpas_v = _normalized_quiver_components(fields["ex"][plot_level, :], fields["ey"][plot_level, :])
    exact_u, exact_v = _normalized_quiver_components(exact.ex[plot_level, :], exact.ey[plot_level, :])
    _overlay_quiver(axes.flat[1], x, y, mpas_u, mpas_v)
    _overlay_quiver(axes.flat[5], x, y, exact_u, exact_v)
    fig.suptitle(f"Free-space charge diagnostics: {source}, level {plot_level}, z={KM_PER_M * plot_z:g} km")
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
        help="Horizontal interior comparison margin in meters for the analysis step.",
    )
    parser.add_argument(
        "--vertical-boundary-margin",
        type=float,
        default=0.0,
        help="Vertical interior comparison margin in meters for the analysis step.",
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
    if args.analysis_only:
        run_dir = args.run_dir.expanduser().resolve()
        result = analyze_output(
            run_dir / "output.nc",
            args.summary,
            args.plot,
            source=args.source,
            charge=args.charge,
            sigma=args.sigma,
            separation=args.separation,
            boundary_margin=args.boundary_margin,
            vertical_boundary_margin=args.vertical_boundary_margin,
            core_radius=args.core_radius,
            e_relative_floor=args.e_relative_floor,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        print(f"Wrote summary: {pathlib.Path(args.summary).expanduser()}")
        print(f"Wrote plot: {pathlib.Path(args.plot).expanduser()}")
        return 0

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
        vertical_boundary_margin=args.vertical_boundary_margin,
        core_radius=args.core_radius,
        e_relative_floor=args.e_relative_floor,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Wrote summary: {pathlib.Path(args.summary).expanduser()}")
    print(f"Wrote plot: {pathlib.Path(args.plot).expanduser()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
