#!/usr/bin/env python3
"""Free-space analytic fields for smooth Gaussian electrostatic charges."""

from dataclasses import dataclass
import math

import numpy as np


EPSILON0 = 8.8541878128e-12
SOURCES = ("gaussian_monopole", "gaussian_dipole_y", "gaussian_dipole_z")


@dataclass(frozen=True)
class GaussianField:
    rho: np.ndarray
    phi: np.ndarray
    ex: np.ndarray
    ey: np.ndarray
    ez: np.ndarray
    centers: tuple[tuple[float, float, float], ...]
    charges: tuple[float, ...]


def _as_float_array(values):
    return np.asarray(values, dtype=float)


def _erf_array(values):
    return np.vectorize(math.erf, otypes=[float])(values)


def gaussian_lobe_density(x, y, z, *, center, charge, sigma):
    x = _as_float_array(x)
    y = _as_float_array(y)
    z = _as_float_array(z)
    cx, cy, cz = center
    r2 = (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2
    norm = (2.0 * math.pi) ** 1.5 * sigma**3
    return charge / norm * np.exp(-r2 / (2.0 * sigma**2))


def gaussian_lobe_phi_e(x, y, z, *, center, charge, sigma):
    x = _as_float_array(x)
    y = _as_float_array(y)
    z = _as_float_array(z)
    cx, cy, cz = center
    dx = x - cx
    dy = y - cy
    dz = z - cz
    r = np.sqrt(dx * dx + dy * dy + dz * dz)
    prefactor = charge / (4.0 * math.pi * EPSILON0)
    a = 1.0 / (math.sqrt(2.0) * sigma)

    phi = np.empty_like(r, dtype=float)
    ex = np.zeros_like(r, dtype=float)
    ey = np.zeros_like(r, dtype=float)
    ez = np.zeros_like(r, dtype=float)

    nonzero = r > 0.0
    erf_values = _erf_array(a * r[nonzero])
    phi[nonzero] = prefactor * erf_values / r[nonzero]
    phi[~nonzero] = prefactor * math.sqrt(2.0 / math.pi) / sigma

    bracket = np.zeros_like(r, dtype=float)
    bracket[nonzero] = (
        erf_values / r[nonzero] ** 2
        - math.sqrt(2.0 / math.pi)
        * np.exp(-(r[nonzero] ** 2) / (2.0 * sigma**2))
        / (sigma * r[nonzero])
    )
    factor = np.zeros_like(r, dtype=float)
    factor[nonzero] = prefactor * bracket[nonzero] / r[nonzero]
    ex[nonzero] = factor[nonzero] * dx[nonzero]
    ey[nonzero] = factor[nonzero] * dy[nonzero]
    ez[nonzero] = factor[nonzero] * dz[nonzero]
    return phi, ex, ey, ez


def source_lobes(source, *, center, charge, separation):
    cx, cy, cz = center
    if source == "gaussian_monopole":
        return ((cx, cy, cz),), (charge,)
    if source == "gaussian_dipole_y":
        half = 0.5 * separation
        return ((cx, cy - half, cz), (cx, cy + half, cz)), (charge, -charge)
    if source == "gaussian_dipole_z":
        half = 0.5 * separation
        return ((cx, cy, cz - half), (cx, cy, cz + half)), (charge, -charge)
    raise ValueError(f"unsupported Gaussian source: {source}")


def evaluate_gaussian_source(source, x, y, z, *, center, charge, sigma, separation):
    centers, charges = source_lobes(source, center=center, charge=charge, separation=separation)
    x = _as_float_array(x)
    y = _as_float_array(y)
    z = _as_float_array(z)
    rho = np.zeros_like(x, dtype=float)
    phi = np.zeros_like(x, dtype=float)
    ex = np.zeros_like(x, dtype=float)
    ey = np.zeros_like(x, dtype=float)
    ez = np.zeros_like(x, dtype=float)

    for lobe_center, lobe_charge in zip(centers, charges):
        rho = rho + gaussian_lobe_density(x, y, z, center=lobe_center, charge=lobe_charge, sigma=sigma)
        lobe_phi, lobe_ex, lobe_ey, lobe_ez = gaussian_lobe_phi_e(
            x, y, z, center=lobe_center, charge=lobe_charge, sigma=sigma
        )
        phi = phi + lobe_phi
        ex = ex + lobe_ex
        ey = ey + lobe_ey
        ez = ez + lobe_ez

    return GaussianField(rho=rho, phi=phi, ex=ex, ey=ey, ez=ez, centers=centers, charges=charges)


def interior_comparison_mask(x, y, z, *, bounds, centers, sigma, boundary_margin, core_radius=None):
    x = _as_float_array(x)
    y = _as_float_array(y)
    z = _as_float_array(z)
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = bounds
    mask = (
        (x >= xmin + boundary_margin)
        & (x <= xmax - boundary_margin)
        & (y >= ymin + boundary_margin)
        & (y <= ymax - boundary_margin)
        & (z >= zmin + boundary_margin)
        & (z <= zmax - boundary_margin)
    )
    exclusion_radius = core_radius if core_radius is not None else 3.0 * sigma
    for cx, cy, cz in centers:
        r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2)
        mask = mask & (r >= exclusion_radius)
    return mask


def align_potential_gauge(phi_mpas, phi_exact, weights, mask):
    phi_mpas = _as_float_array(phi_mpas)
    phi_exact = _as_float_array(phi_exact)
    weights = _as_float_array(weights)
    mask = np.asarray(mask, dtype=bool)
    if not np.any(mask):
        raise ValueError("comparison mask is empty")
    offset = float(np.sum(weights[mask] * (phi_mpas[mask] - phi_exact[mask])) / np.sum(weights[mask]))
    return phi_mpas - offset, offset


def weighted_error_norms(actual, exact, weights, mask, *, relative_floor=0.0):
    actual = _as_float_array(actual)
    exact = _as_float_array(exact)
    weights = _as_float_array(weights)
    mask = np.asarray(mask, dtype=bool)
    if not np.any(mask):
        raise ValueError("comparison mask is empty")
    diff = actual[mask] - exact[mask]
    numerator = float(np.sum(weights[mask] * diff * diff))
    denominator = float(np.sum(weights[mask] * exact[mask] * exact[mask]))
    if relative_floor > 0.0:
        denominator = max(denominator, float(np.sum(weights[mask])) * relative_floor**2)
    l2_absolute = math.sqrt(numerator)
    l2_relative = math.sqrt(numerator / denominator) if denominator > 0.0 else math.nan
    return {
        "l2_absolute": l2_absolute,
        "l2_relative": l2_relative,
        "linf_absolute": float(np.max(np.abs(diff))),
        "n": int(np.count_nonzero(mask)),
    }
