#!/usr/bin/env python3
"""Plot terrain-aware charge-coupled supercell electrostatic figures."""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import plot_charge_coupled_supercell as base  # noqa: E402
import plot_style  # noqa: E402
import run_tripole_supercell as tripole  # noqa: E402

KM_PER_M = base.KM_PER_M
TERRAIN_LINE_COLOR = plot_style.NCAR_COLORS["gray"]
TERRAIN_FILL_COLOR = plot_style.NCAR_COLORS["light_gray"]


def terrain_slice_section(fields):
    indices = base.select_cross_section(fields)
    ordered = indices[np.argsort(np.asarray(fields["yCell"], dtype=float)[indices])]
    return {
        "terrain_y": np.asarray(fields["yCell"], dtype=float)[ordered],
        "terrain_z": np.asarray(fields["zgrid"], dtype=float)[0, ordered],
    }


def terrain_xmean_section(fields, core_half_width_m):
    x = np.asarray(fields["xCell"], dtype=float)
    y = np.asarray(fields["yCell"], dtype=float)
    area = np.asarray(fields.get("areaCell", np.ones(x.size)), dtype=float)
    core_x = base.storm_core_x(fields)
    half_width = max(float(core_half_width_m), 0.0)
    indices = np.where(np.abs(x - core_x) <= half_width)[0]
    if indices.size == 0:
        indices = np.array([int(np.nanargmin(np.abs(x - core_x)))])

    groups = base.y_row_groups(y, indices)
    terrain = np.asarray(fields["zgrid"], dtype=float)[0, :]
    terrain_y = np.asarray(
        [np.average(y[group], weights=base.clean_weights(area, group)) for group in groups],
        dtype=float,
    )
    terrain_z = np.asarray(
        [np.average(terrain[group], weights=base.clean_weights(area, group)) for group in groups],
        dtype=float,
    )
    return {"terrain_y": terrain_y, "terrain_z": terrain_z}


def terrain_section(fields, section_mode, core_half_width_km):
    if section_mode == "slice":
        return terrain_slice_section(fields)
    if section_mode == "xmean":
        return terrain_xmean_section(fields, core_half_width_km / KM_PER_M)
    raise ValueError(f"Unsupported section mode {section_mode!r}")


def overlay_terrain(ax, terrain):
    y_km = np.asarray(terrain["terrain_y"], dtype=float) * KM_PER_M
    z_km = np.asarray(terrain["terrain_z"], dtype=float) * KM_PER_M
    order = np.argsort(y_km)
    y_km = y_km[order]
    z_km = z_km[order]
    ax.fill_between(y_km, 0.0, z_km, color=TERRAIN_FILL_COLOR, alpha=0.35, linewidth=0.0)
    (line,) = ax.plot(y_km, z_km, color=TERRAIN_LINE_COLOR, linewidth=1.4, label="Terrain")
    return line


def plot_terrain_charge_coupled_output(
    output_nc,
    plot_path,
    time_index=-1,
    *,
    section_mode="slice",
    core_half_width_km=10.0,
):
    """Plot terrain-aware water/charge and potential/field vertical sections."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fields = base.read_figure_fields(output_nc, time_index=time_index)
    section = base.select_plot_section(fields, section_mode, core_half_width_km)
    terrain = terrain_section(fields, section_mode, core_half_width_km)
    y_km = section["y"] * KM_PER_M
    z_km = section["z"] * KM_PER_M

    liquid = section["liquid"] * base.G_PER_KG
    rho = section["rho"] * base.NC_PER_C
    phi = section["phi"] * base.MV_PER_V
    e_y = section["e_y"]
    e_z = section["e_z"]

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.35), sharey=True, constrained_layout=True)
    liquid_levels, liquid_ticks = base.positive_filled_levels(liquid)
    phi_levels, phi_ticks, phi_norm = base.signed_filled_levels(phi)
    title_parts = []
    if fields["simulation_time_label"]:
        title_parts.append(fields["simulation_time_label"])
    if section["label"]:
        title_parts.append(section["label"])
    if title_parts:
        fig.suptitle("; ".join(title_parts), fontsize=10, y=1.02)

    base.contourf_from_samples(
        fig,
        axes[0],
        y_km,
        z_km,
        liquid,
        levels=liquid_levels,
        cmap=base.lwc_colormap(),
        norm=None,
        colorbar_label=base.LWC_COLORBAR_LABEL,
        ticks=liquid_ticks,
    )
    base.overlay_signed_contours(axes[0], y_km, z_km, rho)
    axes[0].set_title(base.LEFT_PANEL_TITLE, pad=8)
    base.style_section_axis(axes[0])
    legend_handles = [
        Line2D(
            [0],
            [0],
            color=base.NEGATIVE_CHARGE_CONTOUR_COLORS[-1],
            lw=0.9,
            ls=base.CHARGE_CONTOUR_LINESTYLE,
            label=r"$\rho_q < 0$",
        ),
        Line2D(
            [0],
            [0],
            color=base.POSITIVE_CHARGE_CONTOUR_COLORS[-1],
            lw=0.9,
            ls=base.CHARGE_CONTOUR_LINESTYLE,
            label=r"$\rho_q > 0$",
        ),
    ]
    axes[0].legend(
        handles=legend_handles,
        title=base.RHO_CONTOUR_LEGEND_TITLE,
        loc="upper right",
        fontsize=8,
        title_fontsize=8,
        frameon=True,
        framealpha=0.86,
        borderpad=0.35,
        handlelength=1.8,
    )

    base.contourf_from_samples(
        fig,
        axes[1],
        y_km,
        z_km,
        phi,
        levels=phi_levels,
        cmap=base.phi_colormap(),
        norm=phi_norm,
        colorbar_label=base.PHI_COLORBAR_LABEL,
        ticks=phi_ticks,
    )
    tripole.overlay_field_lines(
        axes[1],
        y_km.ravel(),
        z_km.ravel(),
        e_y.ravel(),
        e_z.ravel(),
        nx=75,
        ny=55,
    )
    axes[1].set_title(base.RIGHT_PANEL_TITLE, pad=8)
    base.style_section_axis(axes[1])

    for ax in axes:
        overlay_terrain(ax, terrain)

    zmin = min(float(np.nanmin(z_km)), float(np.nanmin(terrain["terrain_z"]) * KM_PER_M), 0.0)
    zmax = float(np.nanmax(z_km))
    for ax in axes:
        ax.set_xlim(float(np.nanmin(y_km)), float(np.nanmax(y_km)))
        ax.set_ylim(zmin, zmax)

    plot_path = pathlib.Path(plot_path).expanduser()
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, bbox_inches="tight", dpi=220)
    plt.close(fig)
    return plot_path


def default_plot_path(output_nc, *, section_mode="slice"):
    output_nc = pathlib.Path(output_nc).expanduser()
    name = f"terrain_charge_coupled_{section_mode}.png"
    if output_nc.parent.name == "run":
        return output_nc.parent.parent / "results" / name
    return output_nc.with_name(name)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-nc",
        type=pathlib.Path,
        default=pathlib.Path(
            "~/Data/MPAS/poisson_charge_coupled_supercell_terrain_h1000_v2/run/output.nc"
        ),
        help="Terrain charge-coupled supercell output.nc file.",
    )
    parser.add_argument(
        "--plot",
        type=pathlib.Path,
        default=None,
        help="Figure path. Defaults to a results directory next to the run.",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time record to plot, using Python indexing.",
    )
    parser.add_argument(
        "--section-mode",
        choices=("slice", "xmean"),
        default="slice",
        help="Vertical section mode: strongest-column x slice or core-window x mean.",
    )
    parser.add_argument(
        "--core-half-width-km",
        type=float,
        default=10.0,
        help="Half-width of the x-mean storm-core window in km.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    output_nc = args.output_nc.expanduser()
    plot_path = args.plot.expanduser() if args.plot else default_plot_path(
        output_nc, section_mode=args.section_mode
    )
    plot = plot_terrain_charge_coupled_output(
        output_nc,
        plot_path,
        time_index=args.time_index,
        section_mode=args.section_mode,
        core_half_width_km=args.core_half_width_km,
    )
    print(f"Wrote {plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
