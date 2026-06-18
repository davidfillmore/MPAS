#!/usr/bin/env python3
"""
Prototype Tier A.2 spherical MMS defect-neighborhood corrections.

This script is intentionally diagnostic. It does not modify MPAS inputs or
Fortran source. It reads MPAS spherical meshes and compares the baseline
two-point operator-only residual for the current Y_4^2 target against two
experimental corrections:

- positive shared edge factors near non-hex cells, preserving the two-point
  symmetric operator structure;
- local tangent-plane quadratic least-squares replacement at selected cells,
  which is useful as a non-production fallback test of whether the defect
  error is locally correctable.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import sys
import warnings
from collections import deque
from types import SimpleNamespace

warnings.filterwarnings(
    "ignore",
    message="numpy.ndarray size changed.*",
    category=RuntimeWarning,
)

import netCDF4 as nc
import numpy as np


EARTH_RADIUS_M = 6371229.0
UNREACHED = 10_000


def read_var(dataset, name):
    """Return a NetCDF variable as a NumPy array."""
    return np.asarray(dataset.variables[name][:])


def cell_major(array, n_cells):
    """Return a connectivity array with shape (nCells, maxEdges)."""
    result = np.asarray(array)
    if result.shape[0] != n_cells:
        result = result.T
    return result


def graph_distance_from_defects(n_edges_on_cell, cells_on_cell, max_ring):
    """Return graph distance from cells whose edge count differs from six."""
    n_edges_on_cell = np.asarray(n_edges_on_cell, dtype=int)
    cells_on_cell = np.asarray(cells_on_cell, dtype=int)
    distances = np.full(n_edges_on_cell.size, UNREACHED, dtype=int)
    queue = deque()

    for cell in np.where(n_edges_on_cell != 6)[0]:
        distances[cell] = 0
        queue.append(int(cell))

    while queue:
        cell = queue.popleft()
        if distances[cell] >= max_ring:
            continue
        for edge_index in range(min(n_edges_on_cell[cell], cells_on_cell.shape[1])):
            neighbor = int(cells_on_cell[cell, edge_index])
            if neighbor < 0 or neighbor >= n_edges_on_cell.size:
                continue
            if distances[neighbor] > distances[cell] + 1:
                distances[neighbor] = distances[cell] + 1
                queue.append(neighbor)

    return distances


def apply_edge_factors_to_laplacian(
    phi,
    area,
    edge_weight,
    n_edges_on_cell,
    cells_on_cell,
    edges_on_cell,
    edge_factors=None,
):
    """Apply a cell-centered two-point Laplacian with optional shared edge factors."""
    phi = np.asarray(phi, dtype=float)
    area = np.asarray(area, dtype=float)
    edge_weight = np.asarray(edge_weight, dtype=float)
    n_edges_on_cell = np.asarray(n_edges_on_cell, dtype=int)
    cells_on_cell = np.asarray(cells_on_cell, dtype=int)
    edges_on_cell = np.asarray(edges_on_cell, dtype=int)

    if edge_factors is None:
        edge_factors = np.ones_like(edge_weight)
    else:
        edge_factors = np.asarray(edge_factors, dtype=float)

    laplacian = np.zeros_like(phi, dtype=float)
    max_edges = cells_on_cell.shape[1]
    for edge_index_on_cell in range(max_edges):
        cell_mask = edge_index_on_cell < n_edges_on_cell
        cells = np.where(cell_mask)[0]
        if cells.size == 0:
            continue
        edges = edges_on_cell[cells, edge_index_on_cell]
        neighbors = cells_on_cell[cells, edge_index_on_cell]
        valid = (
            (edges >= 0)
            & (edges < edge_weight.size)
            & (neighbors >= 0)
            & (neighbors < phi.size)
        )
        if not np.any(valid):
            continue
        cells = cells[valid]
        edges = edges[valid]
        neighbors = neighbors[valid]
        laplacian[cells] += (
            edge_factors[edges]
            * edge_weight[edges]
            * (phi[neighbors] - phi[cells])
            / area[cells]
        )

    return laplacian


def fit_symmetric_edge_factors(
    samples,
    targets,
    area,
    edge_weight,
    n_edges_on_cell,
    cells_on_cell,
    edges_on_cell,
    active_edges,
    fit_cells,
    min_factor,
    max_factor,
):
    """Fit bounded shared edge factors for a set of sample fields."""
    samples = np.asarray(samples, dtype=float)
    targets = np.asarray(targets, dtype=float)
    active_edges = np.asarray(active_edges, dtype=int)
    fit_cells = np.asarray(fit_cells, dtype=int)

    factors = np.ones_like(edge_weight, dtype=float)
    if active_edges.size == 0 or fit_cells.size == 0:
        return SimpleNamespace(
            factors=factors,
            active_edges=active_edges,
            residual_before=np.nan,
            residual_after=np.nan,
            cost=np.nan,
            status=0,
        )

    edge_to_column = {int(edge): column for column, edge in enumerate(active_edges)}
    rows = []
    rhs = []

    for sample, target in zip(samples, targets):
        baseline = apply_edge_factors_to_laplacian(
            sample,
            area,
            edge_weight,
            n_edges_on_cell,
            cells_on_cell,
            edges_on_cell,
        )
        for cell in fit_cells:
            row = np.zeros(active_edges.size, dtype=float)
            for edge_index_on_cell in range(n_edges_on_cell[cell]):
                edge = int(edges_on_cell[cell, edge_index_on_cell])
                neighbor = int(cells_on_cell[cell, edge_index_on_cell])
                column = edge_to_column.get(edge)
                if column is None or neighbor < 0:
                    continue
                row[column] = (
                    edge_weight[edge] * (sample[neighbor] - sample[cell]) / area[cell]
                )
            if np.any(row):
                rows.append(row)
                rhs.append(target[cell] - baseline[cell])

    if not rows:
        return SimpleNamespace(
            factors=factors,
            active_edges=active_edges,
            residual_before=np.nan,
            residual_after=np.nan,
            cost=np.nan,
            status=0,
        )

    matrix = np.vstack(rows)
    rhs = np.asarray(rhs)
    residual_before = float(np.linalg.norm(rhs))
    lower = np.full(active_edges.size, min_factor - 1.0)
    upper = np.full(active_edges.size, max_factor - 1.0)

    try:
        from scipy.optimize import lsq_linear

        result = lsq_linear(matrix, rhs, bounds=(lower, upper), lsmr_tol="auto")
        delta = result.x
        cost = float(result.cost)
        status = int(result.status)
    except Exception:
        delta, *_ = np.linalg.lstsq(matrix, rhs, rcond=None)
        delta = np.clip(delta, lower, upper)
        cost = 0.5 * float(np.linalg.norm(matrix @ delta - rhs) ** 2)
        status = 0

    factors[active_edges] = 1.0 + delta
    residual_after = float(np.linalg.norm(matrix @ delta - rhs))
    return SimpleNamespace(
        factors=factors,
        active_edges=active_edges,
        residual_before=residual_before,
        residual_after=residual_after,
        cost=cost,
        status=status,
    )


def quadratic_laplacian_from_offsets(offsets, center_value, neighbor_values):
    """Estimate the tangent-plane Laplacian from a local quadratic fit."""
    offsets = np.asarray(offsets, dtype=float)
    neighbor_values = np.asarray(neighbor_values, dtype=float)
    if offsets.shape[0] < 5:
        return np.nan

    x = offsets[:, 0]
    y = offsets[:, 1]
    matrix = np.column_stack((x, y, 0.5 * x * x, x * y, 0.5 * y * y))
    rhs = neighbor_values - float(center_value)
    coeffs, *_ = np.linalg.lstsq(matrix, rhs, rcond=None)
    return float(coeffs[2] + coeffs[4])


def tangent_basis(normal):
    """Return two orthonormal tangent vectors at a unit normal."""
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)
    reference = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(normal, reference))) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    east = np.cross(reference, normal)
    east = east / np.linalg.norm(east)
    north = np.cross(normal, east)
    return east, north


def neighbor_offsets(mesh, cell):
    """Return geodesic tangent-plane offsets and neighbor indices for one cell."""
    center = mesh.xyz_cell[cell]
    basis_x, basis_y = tangent_basis(center)
    offsets = []
    neighbors = []
    for edge_index_on_cell in range(mesh.n_edges_on_cell[cell]):
        neighbor = int(mesh.cells_on_cell[cell, edge_index_on_cell])
        if neighbor < 0 or neighbor >= mesh.xyz_cell.shape[0]:
            continue
        neighbor_xyz = mesh.xyz_cell[neighbor]
        dot = float(np.clip(np.dot(center, neighbor_xyz), -1.0, 1.0))
        distance = mesh.sphere_radius * math.acos(dot)
        tangent = neighbor_xyz - dot * center
        tangent_norm = float(np.linalg.norm(tangent))
        if tangent_norm <= 0.0:
            continue
        direction = tangent / tangent_norm
        offsets.append(
            [distance * np.dot(direction, basis_x), distance * np.dot(direction, basis_y)]
        )
        neighbors.append(neighbor)
    return np.asarray(offsets), np.asarray(neighbors, dtype=int)


def neighbor_cells_within_rings(n_edges_on_cell, cells_on_cell, cell, rings):
    """Return sorted neighbors within a graph distance, excluding the center cell."""
    visited = {int(cell)}
    frontier = {int(cell)}
    for _ in range(rings):
        next_frontier = set()
        for current in frontier:
            for edge_index_on_cell in range(
                min(int(n_edges_on_cell[current]), cells_on_cell.shape[1])
            ):
                neighbor = int(cells_on_cell[current, edge_index_on_cell])
                if neighbor < 0 or neighbor >= n_edges_on_cell.size or neighbor in visited:
                    continue
                visited.add(neighbor)
                next_frontier.add(neighbor)
        frontier = next_frontier
        if not frontier:
            break
    visited.remove(int(cell))
    return np.asarray(sorted(visited), dtype=int)


def offsets_to_neighbors(mesh, cell, neighbors):
    """Return geodesic tangent-plane offsets from one cell to selected neighbors."""
    center = mesh.xyz_cell[cell]
    basis_x, basis_y = tangent_basis(center)
    offsets = []
    valid_neighbors = []
    for neighbor in neighbors:
        neighbor_xyz = mesh.xyz_cell[neighbor]
        dot = float(np.clip(np.dot(center, neighbor_xyz), -1.0, 1.0))
        distance = mesh.sphere_radius * math.acos(dot)
        tangent = neighbor_xyz - dot * center
        tangent_norm = float(np.linalg.norm(tangent))
        if tangent_norm <= 0.0:
            continue
        direction = tangent / tangent_norm
        offsets.append(
            [distance * np.dot(direction, basis_x), distance * np.dot(direction, basis_y)]
        )
        valid_neighbors.append(neighbor)
    return np.asarray(offsets), np.asarray(valid_neighbors, dtype=int)


def apply_local_lsq_laplacian(mesh, phi, correction_cells, stencil_rings=1):
    """Apply baseline Laplacian, replacing selected cells with local quadratic fits."""
    laplacian = apply_edge_factors_to_laplacian(
        phi,
        mesh.area,
        mesh.edge_weight,
        mesh.n_edges_on_cell,
        mesh.cells_on_cell,
        mesh.edges_on_cell,
    )
    for cell in correction_cells:
        if stencil_rings <= 1:
            offsets, neighbors = neighbor_offsets(mesh, int(cell))
        else:
            neighbors = neighbor_cells_within_rings(
                mesh.n_edges_on_cell,
                mesh.cells_on_cell,
                int(cell),
                stencil_rings,
            )
            offsets, neighbors = offsets_to_neighbors(mesh, int(cell), neighbors)
        if neighbors.size < 5:
            continue
        local_laplacian = quadratic_laplacian_from_offsets(
            offsets,
            phi[cell],
            phi[neighbors],
        )
        if math.isfinite(local_laplacian):
            laplacian[cell] = local_laplacian
    return laplacian


def _triangle_plane_coords(triangle_xyz):
    """Return planar (3, 2) coordinates for a triangle's vertices in its own plane."""
    triangle_xyz = np.asarray(triangle_xyz, dtype=float)
    edge_a = triangle_xyz[1] - triangle_xyz[0]
    edge_b = triangle_xyz[2] - triangle_xyz[0]
    basis_u = edge_a / np.linalg.norm(edge_a)
    perp = edge_b - np.dot(edge_b, basis_u) * basis_u
    basis_v = perp / np.linalg.norm(perp)
    return np.array(
        [
            [0.0, 0.0],
            [np.dot(edge_a, basis_u), 0.0],
            [np.dot(edge_b, basis_u), np.dot(edge_b, basis_v)],
        ]
    )


def whitney_triangle_mass(triangle_xyz):
    """Return the 3x3 lowest-order Whitney 1-form mass for one triangle.

    Rows/columns follow the local edge order [(0, 1), (1, 2), (0, 2)]. The edge
    element for local edge (i, j) is w = lambda_i grad(lambda_j) -
    lambda_j grad(lambda_i); the entries are the closed-form integrals
    ``M[e, e'] = integral_T w_e . w_e'`` using ``integral lambda_p lambda_q =
    area (1 + delta_pq) / 12`` and constant barycentric gradients.
    """
    plane_xy = _triangle_plane_coords(triangle_xyz)
    vandermonde = np.column_stack((np.ones(3), plane_xy[:, 0], plane_xy[:, 1]))
    inverse = np.linalg.inv(vandermonde)
    gradients = [np.array([inverse[1, i], inverse[2, i]]) for i in range(3)]
    area = 0.5 * abs(
        (plane_xy[1, 0] - plane_xy[0, 0]) * (plane_xy[2, 1] - plane_xy[0, 1])
        - (plane_xy[2, 0] - plane_xy[0, 0]) * (plane_xy[1, 1] - plane_xy[0, 1])
    )

    local_pairs = [(0, 1), (1, 2), (0, 2)]

    def product_integral(p, q):
        return area / 12.0 * (1.0 + (1.0 if p == q else 0.0))

    mass = np.zeros((3, 3))
    for row, (i_a, j_a) in enumerate(local_pairs):
        for column, (i_b, j_b) in enumerate(local_pairs):
            mass[row, column] = (
                product_integral(i_a, i_b) * np.dot(gradients[j_a], gradients[j_b])
                - product_integral(i_a, j_b) * np.dot(gradients[j_a], gradients[i_b])
                - product_integral(j_a, i_b) * np.dot(gradients[i_a], gradients[j_b])
                + product_integral(j_a, j_b) * np.dot(gradients[i_a], gradients[i_b])
            )
    return mass, local_pairs


def _reconstruct_local_triangles(cell_xyz, vertices_xyz):
    """Return (triangle cell triples, generating vertex index) for a neighborhood.

    Each dual vertex (Voronoi corner) is the circumcenter of one Delaunay
    triangle, so its three nearest cell centers recover that triangle.
    """
    triangles = []
    for vertex_index, vertex_xyz in enumerate(vertices_xyz):
        distances = np.linalg.norm(cell_xyz - vertex_xyz, axis=1)
        triple = tuple(sorted(int(c) for c in np.argsort(distances)[:3]))
        triangles.append((triple, vertex_index))
    return triangles


def whitney_hodge_block(cell_xyz, edge_list, vertices_xyz):
    """Return the un-lumped (Whitney 1-form) Hodge block for a local edge set.

    ``cell_xyz`` holds the cell-center (primal vertex) coordinates, ``edge_list``
    holds the primal edges as ``(cell_a, cell_b)`` index pairs, and
    ``vertices_xyz`` holds the dual (Voronoi) corner coordinates. The block is
    the lowest-order Whitney edge-mass matrix on the local Delaunay
    triangulation, symmetrically calibrated so its diagonal reproduces the
    DEC-consistent Hodge ``l_e / d_e`` (the discrete Whitney edge-mass diagonal
    is a fixed constant times ``l_e / d_e`` in the regular limit; the
    calibration removes that constant while preserving the off-diagonal Whitney
    coupling that distinguishes irregular neighborhoods). The result is
    symmetric, SPD by construction, and reduces to the diagonal lumped Hodge
    when the local edges decouple (the regular-hexagon limit).
    """
    cell_xyz = np.asarray(cell_xyz, dtype=float)
    vertices_xyz = np.asarray(vertices_xyz, dtype=float)
    edge_list = np.asarray(edge_list, dtype=int)
    n_edge = edge_list.shape[0]

    edge_index = {}
    for column, (cell_a, cell_b) in enumerate(edge_list):
        edge_index[(min(int(cell_a), int(cell_b)), max(int(cell_a), int(cell_b)))] = column

    triangles = _reconstruct_local_triangles(cell_xyz, vertices_xyz)
    mass = np.zeros((n_edge, n_edge))
    incident_vertices = {column: [] for column in range(n_edge)}

    for triple, vertex_index in triangles:
        local_mass, local_pairs = whitney_triangle_mass(cell_xyz[list(triple)])
        local_columns = []
        for (local_i, local_j) in local_pairs:
            cell_a, cell_b = triple[local_i], triple[local_j]
            local_columns.append(
                edge_index.get((min(cell_a, cell_b), max(cell_a, cell_b)))
            )
        for row in range(3):
            if local_columns[row] is None:
                continue
            incident_vertices[local_columns[row]].append(vertex_index)
            for column in range(3):
                if local_columns[column] is None:
                    continue
                mass[local_columns[row], local_columns[column]] += local_mass[row, column]

    scale = np.ones(n_edge)
    for column, (cell_a, cell_b) in enumerate(edge_list):
        dual_distance = float(np.linalg.norm(cell_xyz[cell_a] - cell_xyz[cell_b]))
        flanking = list(dict.fromkeys(incident_vertices[column]))
        if len(flanking) >= 2:
            primal_length = float(
                np.linalg.norm(vertices_xyz[flanking[0]] - vertices_xyz[flanking[1]])
            )
        elif len(flanking) == 1:
            primal_length = dual_distance / math.sqrt(3.0)
        else:
            primal_length = dual_distance
        target = primal_length / dual_distance if dual_distance > 0.0 else 0.0
        if mass[column, column] > 0.0 and target > 0.0:
            scale[column] = math.sqrt(target / mass[column, column])

    hodge = (scale[:, None] * mass) * scale[None, :]
    hodge = 0.5 * (hodge + hodge.T)
    return hodge


def regular_hex_patch():
    """Return a synthetic regular-hexagon edge neighborhood for Whitney tests.

    The patch is a set of congruent equilateral diamonds (two equilateral
    triangles sharing a primal edge), one per listed edge, placed far enough
    apart that the listed edges decouple - the regular-hexagon limit in which
    the Whitney Hodge reduces to the diagonal lumped Hodge. Exposes
    ``cell_xyz``, ``edge_list``, ``vertices_xyz``, ``edge_len`` (dual/Voronoi
    edge length l_e) and ``edge_dc`` (primal edge length d_e).
    """
    apex_height = math.sqrt(3.0) / 2.0
    circumcenter_offset = math.sqrt(3.0) / 6.0

    cells = []
    edges = []
    vertices = []
    edge_len = []
    edge_dc = []
    for diamond in range(6):
        origin = 10.0 * diamond
        cell_a = len(cells)
        cells.append([origin - 0.5, 0.0, 0.0])
        cell_b = len(cells)
        cells.append([origin + 0.5, 0.0, 0.0])
        cells.append([origin, apex_height, 0.0])
        cells.append([origin, -apex_height, 0.0])
        vertices.append([origin, circumcenter_offset, 0.0])
        vertices.append([origin, -circumcenter_offset, 0.0])
        edges.append([cell_a, cell_b])
        edge_dc.append(1.0)
        edge_len.append(2.0 * circumcenter_offset)

    return SimpleNamespace(
        cell_xyz=np.asarray(cells, dtype=float),
        edge_list=np.asarray(edges, dtype=int),
        vertices_xyz=np.asarray(vertices, dtype=float),
        edge_len=np.asarray(edge_len, dtype=float),
        edge_dc=np.asarray(edge_dc, dtype=float),
    )


def infer_sphere_radius(grid_dataset, init_dataset):
    """Infer the physical sphere radius from mesh/init metadata."""
    for dataset in (init_dataset, grid_dataset):
        if dataset is None:
            continue
        value = getattr(dataset, "sphere_radius", None)
        if value is not None and float(value) > 1000.0:
            return float(value)
    return EARTH_RADIUS_M


def load_mesh(bundle_dir):
    """Load the mesh fields needed by the prototype."""
    bundle_dir = pathlib.Path(bundle_dir)
    grid_path = bundle_dir / "grid.nc"
    init_path = bundle_dir / "init.nc"
    init_dataset = nc.Dataset(init_path) if init_path.exists() else None
    try:
        with nc.Dataset(grid_path) as grid:
            sphere_radius = infer_sphere_radius(grid, init_dataset)
            lat = read_var(grid, "latCell")
            lon = read_var(grid, "lonCell")
            x_cell = read_var(grid, "xCell")
            y_cell = read_var(grid, "yCell")
            z_cell = read_var(grid, "zCell")
            area_cell = read_var(grid, "areaCell")
            dc_edge = read_var(grid, "dcEdge")
            dv_edge = read_var(grid, "dvEdge")
            n_edges_on_cell = read_var(grid, "nEdgesOnCell").astype(int)
            n_cells = n_edges_on_cell.size
            cells_on_cell = cell_major(read_var(grid, "cellsOnCell").astype(int), n_cells) - 1
            edges_on_cell = cell_major(read_var(grid, "edgesOnCell").astype(int), n_cells) - 1

        metric_scale = (
            sphere_radius if sphere_radius > 1000.0 and np.max(np.abs(dc_edge)) < 1000.0 else 1.0
        )
        area = area_cell * metric_scale * metric_scale
        edge_weight = (dv_edge * metric_scale) / (dc_edge * metric_scale)
        xyz_cell = np.column_stack((x_cell, y_cell, z_cell))
        xyz_cell = xyz_cell / np.linalg.norm(xyz_cell, axis=1)[:, None]
        return SimpleNamespace(
            lat=lat,
            lon=lon,
            xyz_cell=xyz_cell,
            area=area,
            edge_weight=edge_weight,
            n_edges_on_cell=n_edges_on_cell,
            cells_on_cell=cells_on_cell,
            edges_on_cell=edges_on_cell,
            sphere_radius=sphere_radius,
            metric_scale=metric_scale,
            n_edges=edge_weight.size,
        )
    finally:
        if init_dataset is not None:
            init_dataset.close()


def y42(lat, lon):
    """Return the current Tier A.2 unnormalized Y_4^2 horizontal target."""
    return np.cos(lat) ** 2 * (7.0 * np.sin(lat) ** 2 - 1.0) * np.cos(2.0 * lon)


def real_spherical_sample(l_degree, m_order, phase, lat, lon):
    """Return an unnormalized real spherical-harmonic-like sample."""
    if l_degree == 4 and m_order == 2 and phase == "cos":
        return y42(lat, lon)

    from scipy.special import lpmv

    sample = lpmv(m_order, l_degree, np.sin(lat))
    if m_order > 0:
        if phase == "sin":
            sample = sample * np.sin(float(m_order) * lon)
        else:
            sample = sample * np.cos(float(m_order) * lon)
    return sample


def fitting_modes(lmax):
    """Return real spherical modes used for the least-squares fit."""
    modes = []
    for l_degree in range(1, lmax + 1):
        modes.append((l_degree, 0, "cos"))
        for m_order in range(1, l_degree + 1):
            modes.append((l_degree, m_order, "cos"))
            modes.append((l_degree, m_order, "sin"))
    return modes


def sample_and_target(mesh, modes):
    """Build normalized sample fields and continuous Laplacian targets."""
    samples = []
    targets = []
    for l_degree, m_order, phase in modes:
        sample = real_spherical_sample(l_degree, m_order, phase, mesh.lat, mesh.lon)
        scale = np.max(np.abs(sample))
        if scale <= 0.0:
            continue
        sample = sample / scale
        samples.append(sample)
        targets.append(-float(l_degree * (l_degree + 1)) / mesh.sphere_radius**2 * sample)
    return np.asarray(samples), np.asarray(targets)


def active_edges_from_distances(mesh, distances, active_rings):
    """Return unique edges incident to cells inside active_rings."""
    active_cells = np.where(distances <= active_rings)[0]
    active_edges = set()
    for cell in active_cells:
        for edge_index_on_cell in range(mesh.n_edges_on_cell[cell]):
            edge = int(mesh.edges_on_cell[cell, edge_index_on_cell])
            if 0 <= edge < mesh.n_edges:
                active_edges.add(edge)
    return np.asarray(sorted(active_edges), dtype=int)


def residual_metrics(laplacian, target, area, mask):
    """Return relative L2 and relative Linf operator residual metrics."""
    denominator = float(np.sum(area[mask] * target[mask] ** 2))
    numerator = float(np.sum(area[mask] * (laplacian[mask] - target[mask]) ** 2))
    l2 = math.sqrt(numerator / denominator) if denominator > 0.0 else math.sqrt(numerator)
    target_linf = float(np.max(np.abs(target[mask]))) if np.any(mask) else 0.0
    error_linf = float(np.max(np.abs(laplacian[mask] - target[mask]))) if np.any(mask) else 0.0
    linf = error_linf / target_linf if target_linf > 0.0 else error_linf
    return l2, linf


def convergence_slope(rows, key):
    """Return the finest-pair convergence slope for a metric key."""
    usable = [row for row in rows if row["h_m"] > 0.0 and row[key] > 0.0]
    usable = sorted(usable, key=lambda row: row["h_m"], reverse=True)
    if len(usable) < 2:
        return math.nan
    coarse = usable[-2]
    fine = usable[-1]
    return math.log(coarse[key] / fine[key]) / math.log(coarse["h_m"] / fine["h_m"])


def mesh_spacing_m(mesh_name):
    """Infer mesh spacing in meters from labels such as 120km."""
    name = mesh_name.lower()
    if name.endswith("km"):
        return float(name[:-2]) * 1000.0
    if name.endswith("m"):
        return float(name[:-1])
    return math.nan


def evaluate_mesh(args, mesh_name):
    """Fit and evaluate the correction on one mesh."""
    mesh = load_mesh(args.run_root / "meshes" / mesh_name)
    distances = graph_distance_from_defects(
        mesh.n_edges_on_cell,
        mesh.cells_on_cell,
        max(args.fit_rings, args.active_rings, args.report_exclusion_rings),
    )
    fit_cells = np.where(distances <= args.fit_rings)[0]
    active_edges = active_edges_from_distances(mesh, distances, args.active_rings)
    correction_cells = np.where(distances <= args.active_rings)[0]
    if args.method == "edge-factors":
        modes = fitting_modes(args.fit_lmax)
        samples, targets = sample_and_target(mesh, modes)
        fit = fit_symmetric_edge_factors(
            samples,
            targets,
            mesh.area,
            mesh.edge_weight,
            mesh.n_edges_on_cell,
            mesh.cells_on_cell,
            mesh.edges_on_cell,
            active_edges,
            fit_cells,
            args.min_factor,
            args.max_factor,
        )
    else:
        fit = SimpleNamespace(
            factors=np.ones_like(mesh.edge_weight),
            active_edges=active_edges,
            residual_before=np.nan,
            residual_after=np.nan,
        )

    target_sample = y42(mesh.lat, mesh.lon)
    target_sample = target_sample / np.max(np.abs(target_sample))
    target_laplacian = -20.0 / mesh.sphere_radius**2 * target_sample
    baseline_laplacian = apply_edge_factors_to_laplacian(
        target_sample,
        mesh.area,
        mesh.edge_weight,
        mesh.n_edges_on_cell,
        mesh.cells_on_cell,
        mesh.edges_on_cell,
    )
    if args.method == "edge-factors":
        corrected_laplacian = apply_edge_factors_to_laplacian(
            target_sample,
            mesh.area,
            mesh.edge_weight,
            mesh.n_edges_on_cell,
            mesh.cells_on_cell,
            mesh.edges_on_cell,
            fit.factors,
        )
    else:
        corrected_laplacian = apply_local_lsq_laplacian(
            mesh,
            target_sample,
            correction_cells,
            stencil_rings=args.lsq_stencil_rings,
        )

    rows = []
    for exclusion in [-1] + list(range(args.report_exclusion_rings + 1)):
        if exclusion < 0:
            mask = np.ones(mesh.n_edges_on_cell.size, dtype=bool)
            label = "global"
        else:
            mask = distances > exclusion
            label = f"exclude_through_{exclusion}"
        if not np.any(mask):
            continue
        baseline_l2, baseline_linf = residual_metrics(
            baseline_laplacian, target_laplacian, mesh.area, mask
        )
        corrected_l2, corrected_linf = residual_metrics(
            corrected_laplacian, target_laplacian, mesh.area, mask
        )
        rows.append(
            {
                "mesh": mesh_name,
                "h_m": mesh_spacing_m(mesh_name),
                "exclusion": label,
                "baseline_l2": baseline_l2,
                "corrected_l2": corrected_l2,
                "baseline_linf": baseline_linf,
                "corrected_linf": corrected_linf,
                "fit_cells": fit_cells.size,
                "correction_cells": correction_cells.size,
                "active_edges": active_edges.size,
                "method": args.method,
                "factor_min": float(np.min(fit.factors[active_edges]))
                if active_edges.size
                else 1.0,
                "factor_max": float(np.max(fit.factors[active_edges]))
                if active_edges.size
                else 1.0,
                "fit_residual_before": fit.residual_before,
                "fit_residual_after": fit.residual_after,
            }
        )
    return rows


def write_csv(rows, output_path):
    """Write prototype result rows to CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = (
        "mesh",
        "h_m",
        "exclusion",
        "baseline_l2",
        "corrected_l2",
        "baseline_linf",
        "corrected_linf",
        "fit_cells",
        "correction_cells",
        "active_edges",
        "method",
        "factor_min",
        "factor_max",
        "fit_residual_before",
        "fit_residual_after",
    )
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def print_summary(rows):
    """Print global baseline/corrected metrics and finest-pair slopes."""
    global_rows = [row for row in rows if row["exclusion"] == "global"]
    for row in global_rows:
        print(
            f"{row['mesh']}: baseline L2={row['baseline_l2']:.6e}, "
            f"corrected L2={row['corrected_l2']:.6e}, "
            f"baseline Linf={row['baseline_linf']:.6e}, "
            f"corrected Linf={row['corrected_linf']:.6e}, "
            f"correction_cells={row['correction_cells']}, "
            f"active_edges={row['active_edges']}, "
            f"factor_range=[{row['factor_min']:.3f}, {row['factor_max']:.3f}], "
            f"fit_residual={row['fit_residual_before']:.3e}->{row['fit_residual_after']:.3e}"
        )
    if len(global_rows) >= 2:
        print(
            "Global finest-pair slopes: "
            f"baseline L2={convergence_slope(global_rows, 'baseline_l2'):.3f}, "
            f"corrected L2={convergence_slope(global_rows, 'corrected_l2'):.3f}, "
            f"baseline Linf={convergence_slope(global_rows, 'baseline_linf'):.3f}, "
            f"corrected Linf={convergence_slope(global_rows, 'corrected_linf'):.3f}"
        )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-root",
        type=pathlib.Path,
        required=True,
        help="Tier A.2 run root containing meshes/<mesh>/grid.nc.",
    )
    parser.add_argument("--mesh-list", nargs="+", required=True, help="Mesh labels to analyze.")
    parser.add_argument("--fit-rings", type=int, default=1, help="Defect rings used as fit rows.")
    parser.add_argument(
        "--active-rings",
        type=int,
        default=1,
        help="Defect rings whose incident edges get fitted factors.",
    )
    parser.add_argument("--fit-lmax", type=int, default=4, help="Maximum spherical degree to fit.")
    parser.add_argument(
        "--lsq-stencil-rings",
        type=int,
        default=1,
        help="Neighbor graph rings used by the local-lsq quadratic fit.",
    )
    parser.add_argument(
        "--method",
        choices=("edge-factors", "local-lsq"),
        default="edge-factors",
        help="Correction prototype to apply.",
    )
    parser.add_argument("--min-factor", type=float, default=0.05, help="Lower factor bound.")
    parser.add_argument("--max-factor", type=float, default=4.0, help="Upper factor bound.")
    parser.add_argument(
        "--report-exclusion-rings",
        type=int,
        default=8,
        help="Maximum defect-ring exclusion reported in the CSV.",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=None,
        help="CSV output path. Defaults to run-root/results/tier_A2_defect_correction_prototype.csv.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.run_root = args.run_root.expanduser().resolve()
    if args.output is None:
        args.output = args.run_root / "results" / "tier_A2_defect_correction_prototype.csv"
    else:
        args.output = args.output.expanduser().resolve()

    rows = []
    for mesh_name in args.mesh_list:
        rows.extend(evaluate_mesh(args, mesh_name))

    write_csv(rows, args.output)
    print_summary(rows)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
