#!/usr/bin/env python3
"""Unit tests for the Tier A.2 defect-correction prototype helpers."""

import importlib.util
import math
import pathlib
import unittest

import numpy as np
import scipy.sparse.linalg as spla


SCVT_MESH_ROOT = pathlib.Path(
    "~/Data/MPAS/poisson_tier_A2_scvt/meshes"
).expanduser()


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

    def test_whitney_hodge_block_is_spd(self):
        script = load_script()

        geom = script.regular_hex_patch()
        H = script.whitney_hodge_block(
            geom.cell_xyz, geom.edge_list, geom.vertices_xyz
        )
        self.assertTrue(np.allclose(H, H.T, atol=1e-12))
        self.assertTrue(np.all(np.linalg.eigvalsh(H) > 0.0))

    def test_whitney_hodge_diagonal_is_five_sixths_lumped_on_connected_hex(self):
        # The lowest-order Whitney edge-mass diagonal is NOT the lumped Hodge: for
        # an interior edge in the equilateral (regular-hex) limit it equals exactly
        # 5/6 * (l_e / d_e). Exercise a connected patch whose spokes share fan
        # triangles, so the block has genuinely non-zero off-diagonal coupling.
        script = load_script()

        geom = script.connected_hex_patch()
        H = script.whitney_hodge_block(
            geom.cell_xyz, geom.edge_list, geom.vertices_xyz
        )

        self.assertTrue(np.allclose(H, H.T, atol=1e-12))
        self.assertTrue(np.all(np.linalg.eigvalsh(H) > 0.0))

        off_diagonal = H - np.diag(np.diag(H))
        self.assertGreater(np.max(np.abs(off_diagonal)), 1e-3 * np.max(np.diag(H)))

        expected_diagonal = (5.0 / 6.0) * (geom.edge_len / geom.edge_dc)
        self.assertTrue(np.allclose(np.diag(H), expected_diagonal, rtol=1e-9))

    def test_whitney_hodge_block_requires_geometry(self):
        script = load_script()

        with self.assertRaises(ValueError):
            script.whitney_hodge_block(
                np.zeros((3, 3)), np.array([[0, 1]], dtype=int)
            )

    @unittest.skipUnless(
        (SCVT_MESH_ROOT / "480km" / "grid.nc").exists(),
        "480km SCVT mesh bundle not available",
    )
    def test_enriched_operator_is_spd_on_scvt(self):
        script = load_script()

        mesh = script.load_mesh(SCVT_MESH_ROOT / "480km")
        A = script.assemble_enriched_laplacian(mesh, k_rings=1)

        # Symmetric to machine precision.
        self.assertEqual((abs(A - A.T) > 1e-9).nnz, 0)

        # Positive-definite: smallest algebraic eigenvalue strictly positive.
        lam_min = spla.eigsh(A, k=1, which="SA", return_eigenvectors=False)[0]
        self.assertGreater(lam_min, 0.0)

    @unittest.skipUnless(
        (SCVT_MESH_ROOT / "60km" / "grid.nc").exists(),
        "60km SCVT mesh bundle not available",
    )
    @unittest.expectedFailure  # NO-GO 2026-06-18; see notes/2026-06-18-poisson-defect-operator-no-go.md
    def test_enriched_restores_second_order_Y42(self):
        """Phase-A gate: pure Whitney-Hodge enrichment does NOT restore second
        order (achieved L2 ~0.002 vs the >=1.9 target). The un-lumped defect block
        is inconsistent with the lumped bulk Hodge at the defect/bulk interface,
        giving an O(1/h) ring-1 error. Defect-aware operator track shelved;
        retained as an expected failure documenting the gate verdict."""
        script = load_script()

        slopes = script.convergence_Y42(
            meshes=["480km", "240km", "120km", "60km"], enrich=True, k_rings=1
        )
        self.assertGreaterEqual(slopes["l2"], 1.9)  # NO-GO: enriched L2 ~0.002

    @unittest.skipUnless(
        (SCVT_MESH_ROOT / "480km" / "grid.nc").exists(),
        "480km SCVT mesh bundle not available",
    )
    def test_cotangent_operator_is_spd_on_scvt(self):
        """The global cotangent (P1 nodal stiffness) operator A = d0^T diag(*1) d0
        is symmetric and SPD once grounded (constant is its only null vector)."""
        script = load_script()

        mesh = script.load_mesh(SCVT_MESH_ROOT / "480km")
        A = script.assemble_cotangent_laplacian_ungrounded(mesh)
        n = mesh.n_edges_on_cell.size
        grounded = A[1:n][:, 1:n].tocsr()

        # Symmetric to machine precision.
        self.assertEqual((abs(grounded - grounded.T) > 1e-9).nnz, 0)

        # Positive-definite: smallest algebraic eigenvalue strictly positive.
        lam_min = spla.eigsh(grounded, k=1, which="SA", return_eigenvectors=False)[0]
        self.assertGreater(lam_min, 0.0)

    @unittest.skipUnless(
        (SCVT_MESH_ROOT / "480km" / "grid.nc").exists(),
        "480km SCVT mesh bundle not available",
    )
    def test_cotangent_matches_baseline_in_bulk(self):
        """In the hexagonal bulk (well-centred Delaunay-dual edges) the accumulated
        cotangent Hodge weight 1/2(cot a + cot b) equals the baseline lumped Hodge
        edge_weight = dvEdge/dcEdge to a few percent. The two differ only near the
        pentagons (where circumcentric l_e/d_e is not the true cotangent)."""
        script = load_script()

        mesh = script.load_mesh(SCVT_MESH_ROOT / "480km")
        cot_weight = script.cotangent_edge_weight(mesh)
        baseline = np.asarray(mesh.edge_weight, dtype=float)

        # Bulk edges: both incident cells are hexagons and every cell of every
        # Delaunay triangle touching the edge is a hexagon (away from pentagons).
        bulk = script.bulk_edge_mask(mesh)
        self.assertGreater(int(np.count_nonzero(bulk)), 0)

        relative = np.abs(cot_weight[bulk] - baseline[bulk]) / np.abs(baseline[bulk])
        self.assertLess(float(np.max(relative)), 0.05)

    @unittest.skipUnless(
        (SCVT_MESH_ROOT / "60km" / "grid.nc").exists(),
        "60km SCVT mesh bundle not available",
    )
    @unittest.expectedFailure  # NO-GO 2026-06-18; see notes/2026-06-18-poisson-mfd-bakeoff.md
    def test_cotangent_restores_second_order_Y42(self):
        """THE SPHERE GATE (variant=cotangent). NO-GO: the global, interface-free
        diagonal cotangent Hodge (P1 nodal stiffness on the Delaunay triangulation)
        does NOT restore second order (achieved L2 ~1.22 vs the >=1.9 target).

        Root cause: the icosahedral SCVT mesh is essentially well-centred EVERYWHERE
        - even at the 12 pentagons the dual edge sits at the Delaunay circumcentre,
        so the cotangent identity 1/2(cot a + cot b) = l_e/d_e holds to ~6e-4 on
        every edge. The cotangent and baseline edge weights therefore agree to ~6e-4
        relative on every edge, so the cotangent residual differs from baseline by
        only ~6e-4 and inherits its ~1.2 slope (the raw max|(cot-base) @ Y42| / area
        ~ 3e-16 is a unit artifact of the sphere-area normalization, not operator
        coincidence). The residual is the cell-centred scheme's intrinsic finite-volume
        truncation error on the optimized SCVT, not a pentagon Hodge defect the
        cotangent weight could correct. Retained as an expected failure documenting
        the gate verdict. See notes/2026-06-18-poisson-mfd-bakeoff.md."""
        script = load_script()

        slopes = script.convergence_Y42(
            meshes=["480km", "240km", "120km", "60km"], variant="cotangent"
        )
        self.assertGreaterEqual(slopes["l2"], 1.9)  # NO-GO: cotangent L2 ~1.22

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


class ConvergenceY42SolutionErrorTests(unittest.TestCase):
    @unittest.skipUnless(
        (SCVT_MESH_ROOT / "480km" / "grid.nc").exists(),
        "480km SCVT mesh bundle not available",
    )
    def test_convergence_Y42_solution_error_cotangent_slope_is_positive(self):
        """Solution-error slope for the cotangent sphere operator is finite and positive.
        Diagnostic — does NOT gate on >= 1.9. Records the actual slope."""
        script = load_script()
        result = script.convergence_Y42_solution_error(
            meshes=["480km", "240km", "120km", "60km"], variant="cotangent"
        )
        slope = result["l2_slope"]
        self.assertTrue(math.isfinite(slope), f"cotangent solution-error slope is NaN/inf: {slope}")
        self.assertGreater(slope, 0.0, f"cotangent solution-error slope non-positive: {slope}")

    @unittest.skipUnless(
        (SCVT_MESH_ROOT / "480km" / "grid.nc").exists(),
        "480km SCVT mesh bundle not available",
    )
    def test_convergence_Y42_solution_error_baseline_slope_is_positive(self):
        """Solution-error slope for the baseline sphere operator is finite and positive."""
        script = load_script()
        result = script.convergence_Y42_solution_error(
            meshes=["480km", "240km", "120km", "60km"], variant="baseline"
        )
        slope = result["l2_slope"]
        self.assertTrue(math.isfinite(slope), f"baseline solution-error slope is NaN/inf: {slope}")
        self.assertGreater(slope, 0.0, f"baseline solution-error slope non-positive: {slope}")


if __name__ == "__main__":
    unittest.main()
