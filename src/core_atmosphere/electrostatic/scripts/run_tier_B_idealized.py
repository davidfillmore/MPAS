#!/usr/bin/env python3
"""
Generate Tier B idealized electrostatic paper figures.

The script analyzes existing zero-duration MPAS electrostatic outputs:

* gaussian_monopole: regularized point-charge radial field figure
* tripole: thundercloud-like vertical cross-section figure
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys

import netCDF4 as nc
import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import free_space_charge_analytics as analytic  # noqa: E402
import plot_charge_coupled_supercell as charge_style  # noqa: E402
import run_tier_A1_cartesian_mms as tier_a1  # noqa: E402
from run_tripole_supercell import centerline_indices, overlay_field_lines  # noqa: E402


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
PAPER_DIR = REPO_ROOT.parent / "MPAS-Papers/papers/mpas-poisson-electrostatics"
FIGURE_DIR = PAPER_DIR / "figures"

KM_PER_M = 1.0e-3
KVM_PER_VM = 1.0e-3
MV_PER_V = 1.0e-6
NC_PER_C = 1.0e9
DEFAULT_POINT_CHARGE = 20.0
DEFAULT_POINT_SIGMA = 2000.0
DEFAULT_POINT_MAX_RADIUS = 25000.0
TRIPOLE_PHI_COLORBAR_LABEL = r"$\phi$ (MV)"
TRIPOLE_EMAG_COLORBAR_LABEL = r"$|E|$ (kV m$^{-1}$)"
TRIPOLE_NEGATIVE_CHARGE_COLORS = charge_style.NEGATIVE_CHARGE_CONTOUR_COLORS
TRIPOLE_POSITIVE_CHARGE_COLORS = charge_style.POSITIVE_CHARGE_CONTOUR_COLORS
TRIPOLE_CHARGE_LINESTYLE = charge_style.CHARGE_CONTOUR_LINESTYLE


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate Tier B point-charge and tripole paper figures."
    )
    parser.add_argument(
        "--point-output",
        type=pathlib.Path,
        default=pathlib.Path("~/Data/MPAS/poisson_free_space_charge/gaussian_monopole/run/output.nc"),
        help="MPAS output.nc for the gaussian_monopole case.",
    )
    parser.add_argument(
        "--tripole-output",
        type=pathlib.Path,
        default=pathlib.Path("~/Data/MPAS/poisson_tripole_supercell/run/output.nc"),
        help="MPAS output.nc for the tripole case.",
    )
    parser.add_argument(
        "--point-figure",
        type=pathlib.Path,
        default=FIGURE_DIR / "tier_B1_point_charge.pdf",
        help="Output PDF for the regularized point-charge figure.",
    )
    parser.add_argument(
        "--tripole-figure",
        type=pathlib.Path,
        default=FIGURE_DIR / "tier_B2_tripole.pdf",
        help="Output PDF for the tripole cross-section figure.",
    )
    parser.add_argument(
        "--summary",
        type=pathlib.Path,
        default=pathlib.Path("~/Data/MPAS/poisson_tier_B/results/tier_B_summary.json"),
        help="JSON summary path for paper metrics.",
    )
    parser.add_argument(
        "--point-charge",
        type=float,
        default=DEFAULT_POINT_CHARGE,
        help="Gaussian monopole charge magnitude in C.",
    )
    parser.add_argument(
        "--point-sigma",
        type=float,
        default=DEFAULT_POINT_SIGMA,
        help="Gaussian monopole sigma in m.",
    )
    parser.add_argument(
        "--point-min-radius",
        type=float,
        default=None,
        help="Minimum radius in m for the radial profile fit. Defaults to 3*sigma.",
    )
    parser.add_argument(
        "--point-max-radius",
        type=float,
        default=DEFAULT_POINT_MAX_RADIUS,
        help="Maximum radius in m for the radial profile fit.",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=28,
        help="Number of logarithmic radial bins for the point-charge profile.",
    )
    return parser.parse_args(argv)


def expanded(path):
    return pathlib.Path(path).expanduser()


def read_e_vector(dataset):
    """Read E_vector as (R3, nVertLevels, nCells)."""
    data, dims = tier_a1.read_without_time(dataset, "E_vector")
    required = ("R3", "nVertLevels", "nCells")
    if sorted(dims) != sorted(required):
        raise ValueError(f"E_vector has unsupported dimensions {dims}")
    axes = [dims.index(dim) for dim in required]
    return np.transpose(data, axes)


def read_fields(output_nc):
    """Read fields needed by the Tier B figures."""
    output_nc = expanded(output_nc)
    with nc.Dataset(output_nc) as dataset:
        x = tier_a1.require_cell(dataset, "xCell")
        y = tier_a1.require_cell(dataset, "yCell")
        zgrid = tier_a1.require_level_cell(dataset, "zgrid", "nVertLevelsP1")
        rho = tier_a1.require_level_cell(dataset, "rho_charge", "nVertLevels")
        phi = tier_a1.require_level_cell(dataset, "phi", "nVertLevels")
        e_vector = read_e_vector(dataset)
        residual = tier_a1.scalar_diagnostic(dataset, "cg_residual_final")

    zmid = 0.5 * (zgrid[:-1, :] + zgrid[1:, :])
    e_x = e_vector[0, :, :]
    e_y = e_vector[1, :, :]
    e_z = e_vector[2, :, :]
    e_mag = np.sqrt(e_x * e_x + e_y * e_y + e_z * e_z)
    return {
        "x": x,
        "y": y,
        "zgrid": zgrid,
        "zmid": zmid,
        "rho": rho,
        "phi": phi,
        "E_x": e_x,
        "E_y": e_y,
        "E_z": e_z,
        "E_mag": e_mag,
        "residual": residual,
    }


def domain_center(fields):
    """Return the source center used by current Tier B zero-duration runs."""
    return (
        0.5 * (float(np.min(fields["x"])) + float(np.max(fields["x"]))),
        0.5 * (float(np.min(fields["y"])) + float(np.max(fields["y"]))),
        0.5 * (float(np.min(fields["zgrid"])) + float(np.max(fields["zgrid"]))),
    )


def coulomb_field_strength(radius, charge):
    radius = np.asarray(radius, dtype=float)
    return charge / (4.0 * math.pi * analytic.EPSILON0 * radius * radius)


def gaussian_field_strength(radius, charge, sigma):
    radius = np.asarray(radius, dtype=float)
    _, ex, _, _ = analytic.gaussian_lobe_phi_e(
        radius,
        np.zeros_like(radius),
        np.zeros_like(radius),
        center=(0.0, 0.0, 0.0),
        charge=charge,
        sigma=sigma,
    )
    return np.abs(ex)


def _broadcast_cell_coordinates(fields):
    zmid = np.asarray(fields["zmid"], dtype=float)
    x = np.broadcast_to(np.asarray(fields["x"], dtype=float), zmid.shape)
    y = np.broadcast_to(np.asarray(fields["y"], dtype=float), zmid.shape)
    return x, y, zmid


def point_charge_radial_profile(
    fields,
    *,
    center,
    charge,
    sigma,
    min_radius=None,
    max_radius=None,
    bins=28,
):
    """Return binned radial E-field metrics for the Gaussian monopole."""
    x, y, z = _broadcast_cell_coordinates(fields)
    cx, cy, cz = center
    radius = np.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2).ravel()
    e_mag = np.asarray(fields["E_mag"], dtype=float).ravel()

    r_min = 3.0 * sigma if min_radius is None else float(min_radius)
    r_max = float(np.nanmax(radius)) if max_radius is None else float(max_radius)
    mask = (
        np.isfinite(radius)
        & np.isfinite(e_mag)
        & (radius >= r_min)
        & (radius <= r_max)
        & (e_mag > 0.0)
    )
    if np.count_nonzero(mask) < 3:
        raise ValueError("not enough finite point-charge samples for a radial profile")

    edges = np.geomspace(max(r_min, 1.0), r_max, max(2, int(bins)) + 1)
    binned_radius = []
    binned_e = []
    counts = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        bin_mask = mask & (radius >= lo) & (radius < hi)
        if not np.any(bin_mask):
            continue
        binned_radius.append(float(np.median(radius[bin_mask])))
        binned_e.append(float(np.median(e_mag[bin_mask])))
        counts.append(int(np.count_nonzero(bin_mask)))

    binned_radius = np.asarray(binned_radius, dtype=float)
    binned_e = np.asarray(binned_e, dtype=float)
    if binned_radius.size < 3:
        raise ValueError("not enough populated radial bins for a slope fit")

    slope, intercept = np.polyfit(np.log(binned_radius), np.log(binned_e), deg=1)
    e_coulomb = coulomb_field_strength(binned_radius, charge)
    e_gaussian = gaussian_field_strength(binned_radius, charge, sigma)
    ratio = binned_e / e_coulomb

    return {
        "radius": binned_radius,
        "e_median": binned_e,
        "e_coulomb": e_coulomb,
        "e_gaussian": e_gaussian,
        "counts": np.asarray(counts, dtype=int),
        "far_field_slope": float(slope),
        "fit_intercept": float(intercept),
        "median_coulomb_ratio": float(np.median(ratio)),
        "min_radius": float(r_min),
        "max_radius": float(r_max),
    }


def plot_point_charge(output_nc, figure_path, *, charge, sigma, min_radius, max_radius, bins):
    fields = read_fields(output_nc)
    center = domain_center(fields)
    profile = point_charge_radial_profile(
        fields,
        center=center,
        charge=charge,
        sigma=sigma,
        min_radius=min_radius,
        max_radius=max_radius,
        bins=bins,
    )

    figure_path = expanded(figure_path)
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.3, 3.8))
    ax.loglog(
        profile["radius"] * KM_PER_M,
        profile["e_median"] * KVM_PER_VM,
        "o",
        label="MPAS median bins",
        markersize=4.0,
    )
    ax.loglog(
        profile["radius"] * KM_PER_M,
        profile["e_gaussian"] * KVM_PER_VM,
        "-",
        label="analytic Gaussian",
        linewidth=1.5,
    )
    ax.loglog(
        profile["radius"] * KM_PER_M,
        profile["e_coulomb"] * KVM_PER_VM,
        "--",
        label=r"$Q/(4\pi\varepsilon_0 r^2)$",
        linewidth=1.2,
    )
    ax.axvline(3.0 * sigma * KM_PER_M, color="0.5", linestyle=":", linewidth=1.0)
    ax.set_xlabel("radius from charge center (km)")
    ax.set_ylabel(r"$|E|$ (kV m$^{-1}$)")
    ax.grid(True, which="both", linestyle=":", linewidth=0.4)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Tier B.1 regularized point charge")
    fig.tight_layout()
    fig.savefig(figure_path)
    plt.close(fig)

    return {
        "figure": str(figure_path),
        "charge": float(charge),
        "sigma": float(sigma),
        "center": tuple(float(v) for v in center),
        "cg_residual_final": float(fields["residual"]),
        "far_field_slope": profile["far_field_slope"],
        "median_coulomb_ratio": profile["median_coulomb_ratio"],
        "profile_min_radius": profile["min_radius"],
        "profile_max_radius": profile["max_radius"],
        "profile_bins": int(profile["radius"].size),
    }


def tripole_metrics(fields):
    e_mag = np.asarray(fields["E_mag"], dtype=float)
    flat_index = int(np.nanargmax(e_mag))
    k, i_cell = np.unravel_index(flat_index, e_mag.shape)
    x_key = "x" if "x" in fields else "xCell"
    y_key = "y" if "y" in fields else "yCell"
    z_key = "zmid" if "zmid" in fields else "zMid"
    return {
        "peak_e": float(e_mag[k, i_cell]),
        "peak_e_kv_m": float(e_mag[k, i_cell] * KVM_PER_VM),
        "peak_x": float(fields[x_key][i_cell]),
        "peak_y": float(fields[y_key][i_cell]),
        "peak_z": float(fields[z_key][k, i_cell]),
    }


def _section_indices(fields):
    return centerline_indices(np.asarray(fields["x"], dtype=float), np.asarray(fields["y"], dtype=float))


def tripole_plot_quantities(fields, section):
    """Return Tier B.2 vertical-section quantities in paper plotting units."""
    y = np.broadcast_to(fields["y"][section], fields["zmid"][:, section].shape)
    z = fields["zmid"][:, section]
    return {
        "y_km": y * KM_PER_M,
        "z_km": z * KM_PER_M,
        "phi_mv": fields["phi"][:, section] * MV_PER_V,
        "e_mag_kv_m": fields["E_mag"][:, section] * KVM_PER_VM,
        "rho_nc_m3": fields["rho"][:, section] * NC_PER_C,
        "e_y": fields["E_y"][:, section],
        "e_z": fields["E_z"][:, section],
    }


def tripole_phi_colormap():
    """Return the shared signed potential colormap."""
    return charge_style.phi_colormap()


def tripole_positive_colormap():
    """Return the shared positive-field colormap with white at zero."""
    return charge_style.lwc_colormap()


def _style_tripole_axis(ax):
    ax.set_xlabel("y (km)")
    ax.set_ylabel("z (km)")
    ax.grid(True, linestyle=":", linewidth=0.35, alpha=0.55)


def _add_charge_contour_legend(ax):
    handles = [
        Line2D(
            [0],
            [0],
            color=TRIPOLE_NEGATIVE_CHARGE_COLORS[-1],
            lw=0.9,
            ls=TRIPOLE_CHARGE_LINESTYLE,
            label=r"$\rho_q < 0$",
        ),
        Line2D(
            [0],
            [0],
            color=TRIPOLE_POSITIVE_CHARGE_COLORS[-1],
            lw=0.9,
            ls=TRIPOLE_CHARGE_LINESTYLE,
            label=r"$\rho_q > 0$",
        ),
    ]
    ax.legend(
        handles=handles,
        title=r"$\rho_q$ (nC m$^{-3}$)",
        loc="upper right",
        fontsize=8,
        title_fontsize=8,
        frameon=True,
        framealpha=0.86,
        borderpad=0.35,
        handlelength=1.8,
    )


def _plot_section(ax, y, z, values, title, cmap, contour_values=None, contour_color="k"):
    mesh = ax.contourf(y * KM_PER_M, z * KM_PER_M, values, levels=31, cmap=cmap)
    if contour_values is not None:
        finite = np.asarray(contour_values)[np.isfinite(contour_values)]
        if finite.size and float(np.max(np.abs(finite))) > 0.0:
            bound = float(np.max(np.abs(finite)))
            levels = np.linspace(-bound, bound, 7)
            levels = levels[np.abs(levels) > 0.0]
            if levels.size:
                ax.contour(
                    y * KM_PER_M,
                    z * KM_PER_M,
                    contour_values,
                    levels=levels,
                    colors=contour_color,
                    linewidths=0.45,
                    alpha=0.65,
                )
    ax.set_title(title)
    ax.set_xlabel("y (km)")
    ax.set_ylabel("z (km)")
    return mesh


def plot_tripole(output_nc, figure_path):
    fields = read_fields(output_nc)
    metrics = tripole_metrics(fields)
    section = _section_indices(fields)
    quantities = tripole_plot_quantities(fields, section)
    y_km = quantities["y_km"]
    z_km = quantities["z_km"]
    phi_mv = quantities["phi_mv"]
    e_mag_kv_m = quantities["e_mag_kv_m"]
    rho_nc_m3 = quantities["rho_nc_m3"]
    e_y = quantities["e_y"]
    e_z = quantities["e_z"]

    figure_path = expanded(figure_path)
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.35), sharey=True, constrained_layout=True)

    phi_levels, phi_ticks, phi_norm = charge_style.signed_filled_levels(phi_mv)
    charge_style.contourf_from_samples(
        fig,
        axes[0],
        y_km,
        z_km,
        phi_mv,
        levels=phi_levels,
        cmap=tripole_phi_colormap(),
        norm=phi_norm,
        colorbar_label=TRIPOLE_PHI_COLORBAR_LABEL,
        ticks=phi_ticks,
    )
    charge_style.overlay_signed_contours(axes[0], y_km, z_km, rho_nc_m3)
    overlay_field_lines(
        axes[0],
        y_km.ravel(),
        z_km.ravel(),
        e_y.ravel(),
        e_z.ravel(),
        nx=75,
        ny=55,
    )
    axes[0].set_title("Potential and Electric Field", pad=8)
    _style_tripole_axis(axes[0])
    _add_charge_contour_legend(axes[0])

    emag_levels, emag_ticks = charge_style.positive_filled_levels(e_mag_kv_m)
    charge_style.contourf_from_samples(
        fig,
        axes[1],
        y_km,
        z_km,
        e_mag_kv_m,
        levels=emag_levels,
        cmap=tripole_positive_colormap(),
        norm=None,
        colorbar_label=TRIPOLE_EMAG_COLORBAR_LABEL,
        ticks=emag_ticks,
    )
    charge_style.overlay_signed_contours(axes[1], y_km, z_km, rho_nc_m3)
    axes[1].plot(
        metrics["peak_y"] * KM_PER_M,
        metrics["peak_z"] * KM_PER_M,
        "o",
        color="#00a6d6",
        markeredgecolor="white",
        markeredgewidth=0.7,
        markersize=4.5,
    )
    axes[1].set_title("Electric Field Magnitude", pad=8)
    _style_tripole_axis(axes[1])

    for ax in axes:
        ax.set_xlim(float(np.nanmin(y_km)), float(np.nanmax(y_km)))
        ax.set_ylim(float(np.nanmin(z_km)), float(np.nanmax(z_km)))

    fig.savefig(figure_path, bbox_inches="tight")
    plt.close(fig)

    return {
        "figure": str(figure_path),
        "cg_residual_final": float(fields["residual"]),
        **metrics,
    }


def _json_safe(value):
    if isinstance(value, dict):
        return {str(key): _json_safe(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_summary(path, *, point, tripole):
    path = expanded(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "point": _json_safe(point),
        "tripole": _json_safe(tripole),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main(argv=None):
    args = parse_args(argv)
    point = plot_point_charge(
        args.point_output,
        args.point_figure,
        charge=args.point_charge,
        sigma=args.point_sigma,
        min_radius=args.point_min_radius,
        max_radius=args.point_max_radius,
        bins=args.bins,
    )
    tripole = plot_tripole(args.tripole_output, args.tripole_figure)
    payload = write_summary(args.summary, point=point, tripole=tripole)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
