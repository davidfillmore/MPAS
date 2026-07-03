#!/usr/bin/env python3
"""Plot a two-panel charge-coupled supercell electrostatic figure."""

from __future__ import annotations

import argparse
import pathlib
import re
import sys

import netCDF4 as nc
import numpy as np


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_tripole_supercell as tripole  # noqa: E402


KM_PER_M = 1.0e-3
G_PER_KG = 1.0e3
NC_PER_C = 1.0e9
MV_PER_V = 1.0e-6
NICE_MANTISSAS = (1.0, 2.0, 2.5, 5.0, 10.0)

LWC_COLORBAR_LABEL = "LWC (g m⁻³)"
RHO_CONTOUR_LEGEND_TITLE = r"$\rho_q$ (nC m⁻³)"
PHI_COLORBAR_LABEL = r"$\phi$ (MV)"
LEFT_PANEL_TITLE = "Liquid Water Content and Charge Source"
RIGHT_PANEL_TITLE = "Potential and Electric Field"
NEGATIVE_CHARGE_CONTOUR_COLORS = ("#fcae91", "#fb6a4a", "#de2d26", "#a50f15", "#67000d")
POSITIVE_CHARGE_CONTOUR_COLORS = ("#bdbdbd", "#969696", "#636363", "#252525", "#000000")
CHARGE_CONTOUR_LINESTYLE = "solid"


def read_xtime(dataset, time_index=-1):
    """Read one MPAS xtime record as a stripped ASCII string."""
    if "xtime" not in dataset.variables:
        return None

    variable = dataset.variables["xtime"]
    data = variable[:]
    dims = list(variable.dimensions)
    if "Time" in dims:
        axis = dims.index("Time")
        data = np.take(data, time_index, axis=axis)

    array = np.asarray(data)
    if array.shape == ():
        value = array.item()
        if isinstance(value, bytes):
            return value.decode("ascii", errors="ignore").replace("\x00", "").strip()
        return str(value).replace("\x00", "").strip()

    return array.astype("S1").tobytes().decode("ascii", errors="ignore").replace("\x00", "").strip()


def xtime_to_seconds(xtime):
    """Convert an MPAS xtime string to elapsed seconds from Jan. 1."""
    match = re.match(r"^(\d+)-(\d+)-(\d+)_(\d+):(\d+):(\d+)$", xtime.strip())
    if match is None:
        return None

    year, month, day, hour, minute, second = (int(value) for value in match.groups())
    month_lengths = [31, 29 if year % 4 == 0 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    days = sum(month_lengths[: max(0, month - 1)]) + max(0, day - 1)
    return days * 86400 + hour * 3600 + minute * 60 + second


def xtime_to_simulation_time_label(xtime):
    """Return a compact title-case simulation-time label from MPAS xtime."""
    seconds = xtime_to_seconds(xtime)
    if seconds is None:
        return f"Simulation Time: {xtime}"

    days, remainder = divmod(seconds, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)

    parts = []
    if days:
        parts.append(f"{days} d")
    if hours:
        parts.append(f"{hours} h")
    if minutes:
        parts.append(f"{minutes} min")
    if seconds or not parts:
        parts.append(f"{seconds} s" if seconds else "0 min")
    return "Simulation Time: " + " ".join(parts)


def read_without_time(dataset, name, time_index=-1):
    """Read a variable, selecting one Time record if present."""
    if name not in dataset.variables:
        raise KeyError(f"output.nc does not contain required variable {name}")

    variable = dataset.variables[name]
    data = variable[:]
    dims = list(variable.dimensions)
    if "Time" in dims:
        axis = dims.index("Time")
        data = np.take(data, time_index, axis=axis)
        dims.pop(axis)
    return np.asarray(data), dims


def read_cell(dataset, name):
    """Read an nCells variable."""
    data, dims = read_without_time(dataset, name)
    if dims != ["nCells"]:
        raise ValueError(f"{name} has unsupported dimensions {dims}")
    return data


def read_optional_cell(dataset, name, size):
    """Read an optional nCells variable or return unit weights."""
    if name not in dataset.variables:
        return np.ones(size, dtype=float)
    return read_cell(dataset, name)


def read_level_cell(dataset, name, level_dim="nVertLevels", time_index=-1):
    """Read a level-by-cell variable as (level, cell)."""
    data, dims = read_without_time(dataset, name, time_index=time_index)
    if dims == [level_dim, "nCells"]:
        return data
    if dims == ["nCells", level_dim]:
        return data.T
    raise ValueError(f"{name} has unsupported dimensions {dims}")


def read_e_vector(dataset, time_index=-1):
    """Read E_vector as (R3, nVertLevels, nCells)."""
    data, dims = read_without_time(dataset, "E_vector", time_index=time_index)
    required = ("R3", "nVertLevels", "nCells")
    if sorted(dims) != sorted(required):
        raise ValueError(f"E_vector has unsupported dimensions {dims}")
    order = [dims.index(dim) for dim in required]
    return np.transpose(data, order)


def read_optional_level_cell(dataset, name, shape, time_index=-1):
    """Read an optional level-cell variable or return zeros with shape."""
    if name not in dataset.variables:
        return np.zeros(shape, dtype=float)
    return read_level_cell(dataset, name, time_index=time_index)


def read_figure_fields(output_nc, time_index=-1):
    """Read fields needed for the charge-coupled supercell figure."""
    with nc.Dataset(output_nc) as dataset:
        xtime = read_xtime(dataset, time_index=time_index)
        x_cell = read_cell(dataset, "xCell")
        y_cell = read_cell(dataset, "yCell")
        area_cell = read_optional_cell(dataset, "areaCell", x_cell.size)
        zgrid = read_level_cell(dataset, "zgrid", level_dim="nVertLevelsP1")
        air_density = read_level_cell(dataset, "rho", time_index=time_index)
        qc = read_optional_level_cell(dataset, "qc", (zgrid.shape[0] - 1, x_cell.size), time_index)
        qr = read_optional_level_cell(dataset, "qr", qc.shape, time_index)
        rho_charge = read_level_cell(dataset, "rho_charge", time_index=time_index)
        phi = read_level_cell(dataset, "phi", time_index=time_index)
        e_vector = read_e_vector(dataset, time_index=time_index)

    z_mid = 0.5 * (zgrid[:-1, :] + zgrid[1:, :])
    liquid_water_mixing_ratio = qc + qr
    liquid_water_content = air_density * liquid_water_mixing_ratio
    e_mag = np.sqrt(np.sum(e_vector * e_vector, axis=0))

    return {
        "xCell": x_cell,
        "yCell": y_cell,
        "areaCell": area_cell,
        "zgrid": zgrid,
        "zMid": z_mid,
        "air_density": air_density,
        "liquid_water_mixing_ratio": liquid_water_mixing_ratio,
        "liquid_water_content": liquid_water_content,
        "rho_charge": rho_charge,
        "phi": phi,
        "E_x": e_vector[0, :, :],
        "E_y": e_vector[1, :, :],
        "E_z": e_vector[2, :, :],
        "E_mag": e_mag,
        "time_index": time_index,
        "simulation_time_label": xtime_to_simulation_time_label(xtime) if xtime else None,
    }


def select_cross_section(fields):
    """Return cell indices for an x-slice through the strongest liquid column."""
    x = np.asarray(fields["xCell"], dtype=float)
    y = np.asarray(fields["yCell"], dtype=float)
    x_slice = storm_core_x(fields)
    atol = max(1.0e-8, 1.0e-8 * float(np.ptp(x)))
    indices = np.where(np.isclose(x, x_slice, rtol=1.0e-8, atol=atol))[0]
    if indices.size == 0:
        indices = np.array([int(np.nanargmin(np.abs(x - x_slice)))])
    return indices[np.argsort(y[indices])]


def storm_core_x(fields):
    """Return the x coordinate of the strongest liquid or charge column."""
    x = np.asarray(fields["xCell"], dtype=float)
    liquid = np.asarray(fields["liquid_water_content"], dtype=float)
    rho = np.asarray(fields["rho_charge"], dtype=float)

    liquid_score = np.nansum(np.maximum(liquid, 0.0), axis=0)
    if np.any(np.isfinite(liquid_score) & (liquid_score > 0.0)):
        score = liquid_score
    else:
        score = np.nanmax(np.abs(rho), axis=0)

    peak_cell = int(np.nanargmax(score))
    return float(x[peak_cell])


def section_coordinates(fields, indices):
    """Return y,z sample coordinates for a vertical section."""
    y_section = np.broadcast_to(
        np.asarray(fields["yCell"], dtype=float)[indices][None, :],
        np.asarray(fields["zMid"])[:, indices].shape,
    )
    z_section = np.asarray(fields["zMid"], dtype=float)[:, indices]
    return y_section, z_section


def y_row_groups(y, indices):
    """Group selected cell indices by nearly equal y coordinate."""
    y = np.asarray(y, dtype=float)
    indices = np.asarray(indices, dtype=int)
    if indices.size == 0:
        return []

    ordered = indices[np.argsort(y[indices])]
    tolerance = max(1.0e-8, 1.0e-8 * float(np.ptp(y[indices])))
    groups = []
    current = [int(ordered[0])]
    current_y = float(y[ordered[0]])
    for index in ordered[1:]:
        if abs(float(y[index]) - current_y) <= tolerance:
            current.append(int(index))
        else:
            groups.append(np.asarray(current, dtype=int))
            current = [int(index)]
            current_y = float(y[index])
    groups.append(np.asarray(current, dtype=int))
    return groups


def clean_weights(area, group):
    """Return positive finite area weights for one y-row group."""
    weights = np.asarray(area, dtype=float)[group]
    if weights.size == 0 or not np.all(np.isfinite(weights)) or float(np.sum(weights)) <= 0.0:
        return np.ones(group.size, dtype=float)
    return weights


def weighted_row_mean(values, area, groups):
    """Average level-cell values across x for each y-row group."""
    values = np.asarray(values, dtype=float)
    columns = []
    for group in groups:
        weights = clean_weights(area, group)
        columns.append(np.average(values[:, group], axis=1, weights=weights))
    if not columns:
        return np.empty((values.shape[0], 0), dtype=float)
    return np.column_stack(columns)


def core_x_mean_section(fields, core_half_width_m):
    """Return an area-weighted core-window x-mean y-z section."""
    x = np.asarray(fields["xCell"], dtype=float)
    y = np.asarray(fields["yCell"], dtype=float)
    area = np.asarray(fields.get("areaCell", np.ones(x.size)), dtype=float)
    core_x = storm_core_x(fields)
    half_width = max(float(core_half_width_m), 0.0)
    indices = np.where(np.abs(x - core_x) <= half_width)[0]
    if indices.size == 0:
        indices = np.array([int(np.nanargmin(np.abs(x - core_x)))])

    groups = y_row_groups(y, indices)
    y_values = np.asarray(
        [np.average(y[group], weights=clean_weights(area, group)) for group in groups],
        dtype=float,
    )
    z_values = weighted_row_mean(fields["zMid"], area, groups)
    y_values = np.broadcast_to(y_values[None, :], z_values.shape)
    half_width_km = half_width * KM_PER_M

    return {
        "label": f"Core x Mean, Half-Width: {half_width_km:g} km",
        "y": y_values,
        "z": z_values,
        "liquid": weighted_row_mean(fields["liquid_water_content"], area, groups),
        "rho": weighted_row_mean(fields["rho_charge"], area, groups),
        "phi": weighted_row_mean(fields["phi"], area, groups),
        "e_y": weighted_row_mean(fields["E_y"], area, groups),
        "e_z": weighted_row_mean(fields["E_z"], area, groups),
    }


def slice_section(fields):
    """Return the existing strongest-column x-slice section."""
    section = select_cross_section(fields)
    y_section, z_section = section_coordinates(fields, section)
    return {
        "label": None,
        "y": y_section,
        "z": z_section,
        "liquid": fields["liquid_water_content"][:, section],
        "rho": fields["rho_charge"][:, section],
        "phi": fields["phi"][:, section],
        "e_y": fields["E_y"][:, section],
        "e_z": fields["E_z"][:, section],
    }


def select_plot_section(fields, section_mode, core_half_width_km):
    """Return plot-ready section data for the requested section mode."""
    if section_mode == "slice":
        return slice_section(fields)
    if section_mode == "xmean":
        return core_x_mean_section(fields, core_half_width_km / KM_PER_M)
    raise ValueError(f"Unsupported section mode {section_mode!r}")


def centered_norm(values):
    """Return a zero-centered normalization for signed fields."""
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        return None
    bound = float(np.max(np.abs(finite)))
    if bound <= 0.0:
        return None

    from matplotlib.colors import TwoSlopeNorm

    return TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound)


def nice_step(value, *, round_up=True):
    """Return a 1/2/2.5/5/10-based step near value."""
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        return 1.0

    exponent = np.floor(np.log10(value))
    scale = 10.0**exponent
    scaled = value / scale
    if round_up:
        for mantissa in NICE_MANTISSAS:
            if scaled <= mantissa:
                return mantissa * scale
        return 10.0 * scale

    mantissa = min(NICE_MANTISSAS, key=lambda candidate: abs(np.log10(scaled / candidate)))
    return mantissa * scale


def positive_filled_levels(values, tick_intervals=6, subdivisions=4):
    """Return round filled-contour levels and colorbar ticks for positive fields."""
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        return np.linspace(0.0, 1.0, tick_intervals * subdivisions + 1), np.linspace(
            0.0, 1.0, tick_intervals + 1
        )
    vmax = float(np.max(finite))
    if vmax <= 0.0:
        return np.linspace(0.0, 1.0, tick_intervals * subdivisions + 1), np.linspace(
            0.0, 1.0, tick_intervals + 1
        )

    tick_step = nice_step(vmax / tick_intervals, round_up=True)
    upper = np.ceil(vmax / tick_step) * tick_step
    tick_values = np.arange(0.0, upper + 0.5 * tick_step, tick_step)
    level_step = tick_step / subdivisions
    levels = np.arange(0.0, upper + 0.5 * level_step, level_step)
    return levels, tick_values


def lwc_colormap():
    """Return a liquid-water colormap with white at zero."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "lwc_white_ygnbu",
        ("#ffffff", "#e5f5b9", "#99d8c9", "#2ca25f", "#006d2c", "#08306b"),
    )


def phi_colormap():
    """Return a signed potential colormap with white fixed at zero."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "phi_zero_white_rdbu",
        (
            (0.0, "#2166ac"),
            (0.25, "#67a9cf"),
            (0.5, "#ffffff"),
            (0.75, "#ef8a62"),
            (1.0, "#b2182b"),
        ),
        N=257,
    )


def signed_line_levels(values, negative_count=6, positive_count=4):
    """Return negative and positive contour levels for a signed field."""
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        return np.array([]), np.array([])
    negative_values = np.abs(finite[finite < 0.0])
    positive_values = finite[finite > 0.0]

    negative = np.array([])
    positive = np.array([])
    if negative_values.size:
        negative = -log_nice_levels(float(np.max(negative_values)), minimum=1.0e-2)[::-1]
    if positive_values.size:
        positive = log_nice_levels(float(np.max(positive_values)), minimum=1.0e-3)
    return negative, positive


def log_nice_levels(high, minimum):
    """Return 1/2/5 log-spaced nice levels up to high."""
    if not np.isfinite(high) or high <= 0.0:
        return np.array([])

    lower = max(float(minimum), 10.0 ** np.floor(np.log10(high)) * 1.0e-3)
    start_exp = int(np.floor(np.log10(lower)))
    stop_exp = int(np.ceil(np.log10(high)))
    levels = []
    for exponent in range(start_exp, stop_exp + 1):
        scale = 10.0**exponent
        for mantissa in (1.0, 2.0, 5.0):
            value = mantissa * scale
            if lower <= value <= high:
                levels.append(value)
    return np.asarray(levels, dtype=float)


def contour_colors(levels, colors):
    """Map contour levels from weakest to strongest color."""
    levels = np.asarray(levels, dtype=float)
    if levels.size == 0:
        return []
    if levels.size == 1:
        return [colors[-1]]

    indices = np.linspace(0, len(colors) - 1, levels.size)
    return [colors[int(round(index))] for index in indices]


def signed_norm(lower, upper):
    """Return a norm that keeps zero at the neutral color."""
    from matplotlib.colors import Normalize, TwoSlopeNorm

    if lower < 0.0 and upper > 0.0:
        return TwoSlopeNorm(vmin=lower, vcenter=0.0, vmax=upper)
    if lower < 0.0 and upper == 0.0:
        return TwoSlopeNorm(vmin=lower, vcenter=0.0, vmax=-lower)
    if lower == 0.0 and upper > 0.0:
        return TwoSlopeNorm(vmin=-upper, vcenter=0.0, vmax=upper)
    if upper > lower:
        return Normalize(vmin=lower, vmax=upper)
    return None


def signed_filled_levels(values, tick_intervals=4, subdivisions=4):
    """Return round filled-contour levels, ticks, and norm for signed fields."""
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        levels = np.linspace(-1.0, 1.0, tick_intervals * subdivisions + 1)
        ticks = np.linspace(-1.0, 1.0, tick_intervals + 1)
        return levels, ticks, signed_norm(-1.0, 1.0)

    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    if vmin < 0.0 and vmax > 0.0:
        bound = max(abs(vmin), abs(vmax))
        tick_step = nice_step(bound / tick_intervals, round_up=True)
        upper = np.ceil(bound / tick_step) * tick_step
        lower = -upper
    elif vmin < 0.0:
        tick_step = nice_step(abs(vmin) / tick_intervals, round_up=True)
        lower = -np.ceil(abs(vmin) / tick_step) * tick_step
        upper = 0.0
    elif vmax > 0.0:
        tick_step = nice_step(vmax / tick_intervals, round_up=True)
        lower = 0.0
        upper = np.ceil(vmax / tick_step) * tick_step
    else:
        tick_step = 1.0
        lower = -1.0
        upper = 1.0

    ticks = np.arange(lower, upper + 0.5 * tick_step, tick_step)
    level_step = tick_step / subdivisions
    levels = np.arange(lower, upper + 0.5 * level_step, level_step)
    return levels, ticks, signed_norm(lower, upper)


def contourf_from_samples(fig, ax, y, z, values, *, levels, cmap, norm, colorbar_label, ticks=None):
    """Draw filled contours from unstructured y-z samples."""
    grid = tripole.scalar_grid_from_samples(y.ravel(), z.ravel(), values.ravel(), nx=95, ny=64)
    if grid is None:
        artist = ax.scatter(
            y.ravel(),
            z.ravel(),
            c=values.ravel(),
            s=10.0,
            cmap=cmap,
            norm=norm,
            linewidths=0.0,
        )
    else:
        grid_y, grid_z, grid_values = grid
        artist = ax.contourf(
            grid_y,
            grid_z,
            grid_values,
            levels=levels,
            cmap=cmap,
            norm=norm,
            extend="both",
        )

    cbar = fig.colorbar(artist, ax=ax, fraction=0.052, pad=0.025, ticks=ticks)
    cbar.set_label(colorbar_label)
    return artist


def overlay_signed_contours(ax, y, z, values):
    """Overlay negative and positive charge-density contours."""
    grid = tripole.scalar_grid_from_samples(y.ravel(), z.ravel(), values.ravel(), nx=95, ny=64)
    if grid is None:
        return []

    grid_y, grid_z, grid_values = grid
    vmin = float(np.nanmin(grid_values))
    vmax = float(np.nanmax(grid_values))
    if not vmax > vmin:
        return []

    def usable_levels(levels):
        levels = np.unique(np.asarray(levels, dtype=float))
        return levels[(levels > vmin) & (levels < vmax)]

    negative, positive = signed_line_levels(grid_values)
    artists = []
    negative = usable_levels(negative)
    positive = usable_levels(positive)
    if negative.size:
        contour = ax.contour(
            grid_y,
            grid_z,
            grid_values,
            levels=negative,
            colors=contour_colors(np.abs(negative), NEGATIVE_CHARGE_CONTOUR_COLORS),
            linewidths=0.75,
            linestyles=CHARGE_CONTOUR_LINESTYLE,
        )
        artists.append(contour)
    if positive.size:
        contour = ax.contour(
            grid_y,
            grid_z,
            grid_values,
            levels=positive,
            colors=contour_colors(positive, POSITIVE_CHARGE_CONTOUR_COLORS),
            linewidths=0.75,
            linestyles=CHARGE_CONTOUR_LINESTYLE,
        )
        artists.append(contour)
    return artists


def style_section_axis(ax):
    """Apply common y-z section axis styling."""
    ax.set_xlabel("y (km)")
    ax.set_ylabel("z (km)")
    ax.grid(True, linestyle=":", linewidth=0.35, alpha=0.55)


def plot_charge_coupled_output(
    output_nc,
    plot_path,
    time_index=-1,
    *,
    section_mode="slice",
    core_half_width_km=10.0,
):
    """Plot a two-panel water/charge and potential/field vertical section."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fields = read_figure_fields(output_nc, time_index=time_index)
    section = select_plot_section(fields, section_mode, core_half_width_km)
    y_km = section["y"] * KM_PER_M
    z_km = section["z"] * KM_PER_M

    liquid = section["liquid"] * G_PER_KG
    rho = section["rho"] * NC_PER_C
    phi = section["phi"] * MV_PER_V
    e_y = section["e_y"]
    e_z = section["e_z"]

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.35), sharey=True, constrained_layout=True)
    liquid_levels, liquid_ticks = positive_filled_levels(liquid)
    phi_levels, phi_ticks, phi_norm = signed_filled_levels(phi)
    title_parts = []
    if fields["simulation_time_label"]:
        title_parts.append(fields["simulation_time_label"])
    if section["label"]:
        title_parts.append(section["label"])
    if title_parts:
        fig.suptitle("; ".join(title_parts), fontsize=10, y=1.02)

    contourf_from_samples(
        fig,
        axes[0],
        y_km,
        z_km,
        liquid,
        levels=liquid_levels,
        cmap=lwc_colormap(),
        norm=None,
        colorbar_label=LWC_COLORBAR_LABEL,
        ticks=liquid_ticks,
    )
    overlay_signed_contours(axes[0], y_km, z_km, rho)
    axes[0].set_title(LEFT_PANEL_TITLE, pad=8)
    style_section_axis(axes[0])
    legend_handles = [
        Line2D(
            [0],
            [0],
            color=NEGATIVE_CHARGE_CONTOUR_COLORS[-1],
            lw=0.9,
            ls=CHARGE_CONTOUR_LINESTYLE,
            label=r"$\rho_q < 0$",
        ),
        Line2D(
            [0],
            [0],
            color=POSITIVE_CHARGE_CONTOUR_COLORS[-1],
            lw=0.9,
            ls=CHARGE_CONTOUR_LINESTYLE,
            label=r"$\rho_q > 0$",
        ),
    ]
    axes[0].legend(
        handles=legend_handles,
        title=RHO_CONTOUR_LEGEND_TITLE,
        loc="upper right",
        fontsize=8,
        title_fontsize=8,
        frameon=True,
        framealpha=0.86,
        borderpad=0.35,
        handlelength=1.8,
    )

    contourf_from_samples(
        fig,
        axes[1],
        y_km,
        z_km,
        phi,
        levels=phi_levels,
        cmap=phi_colormap(),
        norm=phi_norm,
        colorbar_label=PHI_COLORBAR_LABEL,
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
    axes[1].set_title(RIGHT_PANEL_TITLE, pad=8)
    style_section_axis(axes[1])

    for ax in axes:
        ax.set_xlim(float(np.nanmin(y_km)), float(np.nanmax(y_km)))
        ax.set_ylim(float(np.nanmin(z_km)), float(np.nanmax(z_km)))

    plot_path = pathlib.Path(plot_path).expanduser()
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return plot_path


def default_plot_path(output_nc):
    """Return the default result path adjacent to a run directory."""
    output_nc = pathlib.Path(output_nc).expanduser()
    if output_nc.parent.name == "run":
        return output_nc.parent.parent / "results" / "charge_coupled_two_panel.png"
    return output_nc.with_name("charge_coupled_two_panel.png")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-nc",
        type=pathlib.Path,
        default=pathlib.Path(
            "~/Data/MPAS/poisson_charge_coupled_supercell/calibrated_1e-6/run/output.nc"
        ),
        help="Charge-coupled supercell output.nc file.",
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
    plot_path = args.plot.expanduser() if args.plot else default_plot_path(output_nc)
    plot = plot_charge_coupled_output(
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
