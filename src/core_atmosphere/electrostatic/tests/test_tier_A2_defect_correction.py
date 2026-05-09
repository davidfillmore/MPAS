#!/usr/bin/env python3
"""Unit tests for the Tier A.2 defect-correction prototype helpers."""

import importlib.util
import pathlib
import unittest

import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/prototype_tier_A2_defect_correction.py"
)


def load_script():
    spec = importlib.util.spec_from_file_location(
        "prototype_tier_A2_defect_correction", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class DefectCorrectionPrototypeTests(unittest.TestCase):
    def test_graph_distance_from_defects_marks_neighbor_rings(self):
        script = load_script()

        n_edges_on_cell = np.array([5, 6, 6, 6], dtype=int)
        cells_on_cell = np.array(
            [
                [1, -1, -1],
                [0, 2, -1],
                [1, 3, -1],
                [2, -1, -1],
            ],
            dtype=int,
        )

        distances = script.graph_distance_from_defects(
            n_edges_on_cell,
            cells_on_cell,
            max_ring=3,
        )

        np.testing.assert_array_equal(distances, np.array([0, 1, 2, 3]))

    def test_apply_edge_factors_to_laplacian_uses_shared_symmetric_edge_weight(self):
        script = load_script()

        phi = np.array([1.0, 3.0])
        area = np.array([2.0, 4.0])
        edge_weight = np.array([2.0])
        n_edges_on_cell = np.array([1, 1], dtype=int)
        cells_on_cell = np.array([[1], [0]], dtype=int)
        edges_on_cell = np.array([[0], [0]], dtype=int)
        edge_factors = np.array([1.5])

        laplacian = script.apply_edge_factors_to_laplacian(
            phi,
            area,
            edge_weight,
            n_edges_on_cell,
            cells_on_cell,
            edges_on_cell,
            edge_factors,
        )

        np.testing.assert_allclose(laplacian, np.array([3.0, -1.5]))

    def test_fit_symmetric_edge_factors_reduces_constrained_residual(self):
        script = load_script()

        area = np.ones(3)
        edge_weight = np.ones(2)
        n_edges_on_cell = np.array([1, 2, 1], dtype=int)
        cells_on_cell = np.array(
            [
                [1, -1],
                [0, 2],
                [1, -1],
            ],
            dtype=int,
        )
        edges_on_cell = np.array(
            [
                [0, -1],
                [0, 1],
                [1, -1],
            ],
            dtype=int,
        )
        active_edges = np.array([0, 1], dtype=int)
        fit_cells = np.array([0, 1, 2], dtype=int)
        samples = np.array([[0.0, 1.0, 0.0]])
        targets = np.array([[0.5, -1.0, 0.5]])

        baseline = script.apply_edge_factors_to_laplacian(
            samples[0],
            area,
            edge_weight,
            n_edges_on_cell,
            cells_on_cell,
            edges_on_cell,
        )
        baseline_norm = np.linalg.norm(baseline[fit_cells] - targets[0, fit_cells])

        result = script.fit_symmetric_edge_factors(
            samples,
            targets,
            area,
            edge_weight,
            n_edges_on_cell,
            cells_on_cell,
            edges_on_cell,
            active_edges,
            fit_cells,
            min_factor=0.05,
            max_factor=4.0,
        )
        corrected = script.apply_edge_factors_to_laplacian(
            samples[0],
            area,
            edge_weight,
            n_edges_on_cell,
            cells_on_cell,
            edges_on_cell,
            result.factors,
        )
        corrected_norm = np.linalg.norm(corrected[fit_cells] - targets[0, fit_cells])

        self.assertGreaterEqual(np.min(result.factors), 0.05)
        self.assertLessEqual(np.max(result.factors), 4.0)
        self.assertLess(corrected_norm, 0.5 * baseline_norm)

    def test_quadratic_laplacian_from_offsets_recovers_planar_polynomial(self):
        script = load_script()

        offsets = np.array(
            [
                [1.0, 0.0],
                [-1.0, 0.0],
                [0.0, 1.0],
                [0.0, -1.0],
                [1.0, 1.0],
                [-1.0, -1.0],
            ]
        )
        center_value = 4.0
        x = offsets[:, 0]
        y = offsets[:, 1]
        neighbor_values = center_value + 3.0 * x - 2.0 * y + x * x + 2.0 * y * y

        laplacian = script.quadratic_laplacian_from_offsets(
            offsets,
            center_value,
            neighbor_values,
        )

        self.assertAlmostEqual(laplacian, 6.0, places=12)

    def test_neighbor_cells_within_rings_expands_from_center_cell(self):
        script = load_script()

        n_edges_on_cell = np.array([1, 2, 2, 1], dtype=int)
        cells_on_cell = np.array(
            [
                [1, -1],
                [0, 2],
                [1, 3],
                [2, -1],
            ],
            dtype=int,
        )

        ring_one = script.neighbor_cells_within_rings(
            n_edges_on_cell,
            cells_on_cell,
            cell=1,
            rings=1,
        )
        ring_two = script.neighbor_cells_within_rings(
            n_edges_on_cell,
            cells_on_cell,
            cell=1,
            rings=2,
        )

        np.testing.assert_array_equal(ring_one, np.array([0, 2]))
        np.testing.assert_array_equal(ring_two, np.array([0, 2, 3]))


if __name__ == "__main__":
    unittest.main()
