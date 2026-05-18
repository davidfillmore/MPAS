#!/usr/bin/env python3
"""Tests for free-space Gaussian charge analytic utilities."""

import importlib.util
import math
import pathlib
import unittest

import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
MODULE_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/free_space_charge_analytics.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("free_space_charge_analytics", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class GaussianAnalyticTests(unittest.TestCase):
    def test_monopole_center_limit_is_finite_and_field_zero(self):
        analytic = load_module()
        phi, ex, ey, ez = analytic.gaussian_lobe_phi_e(
            np.array([0.0]),
            np.array([0.0]),
            np.array([0.0]),
            center=(0.0, 0.0, 0.0),
            charge=2.0,
            sigma=3.0,
        )

        expected_phi = 2.0 / (4.0 * math.pi * analytic.EPSILON0) * math.sqrt(2.0 / math.pi) / 3.0
        self.assertAlmostEqual(float(phi[0]), expected_phi, delta=abs(expected_phi) * 1.0e-14)
        self.assertEqual(float(ex[0]), 0.0)
        self.assertEqual(float(ey[0]), 0.0)
        self.assertEqual(float(ez[0]), 0.0)

    def test_monopole_far_field_approaches_coulomb_solution(self):
        analytic = load_module()
        radius = 100.0
        charge = 4.0
        phi, ex, ey, ez = analytic.gaussian_lobe_phi_e(
            np.array([radius]),
            np.array([0.0]),
            np.array([0.0]),
            center=(0.0, 0.0, 0.0),
            charge=charge,
            sigma=1.0,
        )

        prefactor = charge / (4.0 * math.pi * analytic.EPSILON0)
        self.assertAlmostEqual(float(phi[0]), prefactor / radius, delta=abs(prefactor / radius) * 1.0e-12)
        self.assertAlmostEqual(float(ex[0]), prefactor / radius**2, delta=abs(prefactor / radius**2) * 1.0e-12)
        self.assertAlmostEqual(float(ey[0]), 0.0, delta=1.0e-18)
        self.assertAlmostEqual(float(ez[0]), 0.0, delta=1.0e-18)

    def test_monopole_small_nonzero_radius_matches_linear_field_limit(self):
        analytic = load_module()
        radius = 1.0e-8
        phi, ex, ey, ez = analytic.gaussian_lobe_phi_e(
            np.array([radius]),
            np.array([0.0]),
            np.array([0.0]),
            center=(0.0, 0.0, 0.0),
            charge=1.0,
            sigma=1.0,
        )

        prefactor = 1.0 / (4.0 * math.pi * analytic.EPSILON0)
        expected_ex = prefactor * math.sqrt(2.0 / math.pi) / 3.0 * radius
        self.assertTrue(math.isfinite(float(phi[0])))
        self.assertGreater(float(ex[0]), 0.0)
        self.assertAlmostEqual(float(ex[0]), expected_ex, delta=abs(expected_ex) * 1.0e-8)
        self.assertEqual(float(ey[0]), 0.0)
        self.assertEqual(float(ez[0]), 0.0)

    def test_horizontal_dipole_is_linear_superposition(self):
        analytic = load_module()
        x = np.array([0.0, 0.0])
        y = np.array([-2.0, 2.0])
        z = np.array([0.0, 0.0])
        center = (0.0, 0.0, 0.0)

        dipole = analytic.evaluate_gaussian_source(
            "gaussian_dipole_y",
            x,
            y,
            z,
            center=center,
            charge=5.0,
            sigma=1.0,
            separation=4.0,
        )
        positive = analytic.gaussian_lobe_phi_e(x, y, z, center=(0.0, -2.0, 0.0), charge=5.0, sigma=1.0)
        negative = analytic.gaussian_lobe_phi_e(x, y, z, center=(0.0, 2.0, 0.0), charge=-5.0, sigma=1.0)

        self.assertTrue(np.allclose(dipole.phi, positive[0] + negative[0]))
        self.assertTrue(np.allclose(dipole.ey, positive[2] + negative[2]))
        self.assertGreater(dipole.ey[0], 0.0)
        self.assertGreater(dipole.ey[1], 0.0)

    def test_vertical_dipole_is_antisymmetric_in_potential_and_symmetric_in_ez(self):
        analytic = load_module()
        x = np.array([0.0, 0.0])
        y = np.array([0.0, 0.0])
        z = np.array([-2.0, 2.0])

        dipole = analytic.evaluate_gaussian_source(
            "gaussian_dipole_z",
            x,
            y,
            z,
            center=(0.0, 0.0, 0.0),
            charge=5.0,
            sigma=1.0,
            separation=4.0,
        )

        self.assertAlmostEqual(float(dipole.phi[0]), -float(dipole.phi[1]))
        self.assertAlmostEqual(float(dipole.ez[0]), float(dipole.ez[1]))
        self.assertGreater(float(dipole.ez[0]), 0.0)

    def test_mask_excludes_boundaries_and_charge_cores(self):
        analytic = load_module()
        x = np.array([0.0, 5.0, 10.0, 5.0, 8.0])
        y = np.array([5.0, 5.0, 5.0, 5.0, 8.0])
        z = np.array([5.0, 5.0, 5.0, 1.0, 8.0])
        mask = analytic.interior_comparison_mask(
            x,
            y,
            z,
            bounds=((0.0, 10.0), (0.0, 10.0), (0.0, 10.0)),
            centers=[(5.0, 5.0, 5.0)],
            sigma=1.0,
            boundary_margin=2.0,
            core_radius=2.0,
        )

        self.assertEqual(mask.tolist(), [False, False, False, False, True])

    def test_gauge_alignment_removes_weighted_mean_offset(self):
        analytic = load_module()
        phi_mpas = np.array([11.0, 12.0, 13.0])
        phi_exact = np.array([1.0, 2.0, 3.0])
        weights = np.array([1.0, 2.0, 1.0])
        mask = np.array([True, True, True])

        aligned, offset = analytic.align_potential_gauge(phi_mpas, phi_exact, weights, mask)

        self.assertAlmostEqual(offset, 10.0)
        self.assertTrue(np.allclose(aligned, phi_exact))

    def test_gauge_alignment_rejects_zero_masked_weight_sum(self):
        analytic = load_module()
        phi_mpas = np.array([11.0, 12.0])
        phi_exact = np.array([1.0, 2.0])
        weights = np.array([0.0, 0.0])
        mask = np.array([True, True])

        with self.assertRaises(ValueError):
            analytic.align_potential_gauge(phi_mpas, phi_exact, weights, mask)

    def test_error_norms_are_finite_and_relative(self):
        analytic = load_module()
        exact = np.array([2.0, 4.0])
        actual = np.array([3.0, 6.0])
        weights = np.array([1.0, 1.0])
        mask = np.array([True, True])

        norms = analytic.weighted_error_norms(actual, exact, weights, mask)

        self.assertAlmostEqual(norms["l2_relative"], 0.5)
        self.assertAlmostEqual(norms["linf_absolute"], 2.0)

    def test_error_norms_support_relative_floor(self):
        analytic = load_module()
        exact = np.array([0.0, 0.0])
        actual = np.array([1.0, 2.0])
        weights = np.array([1.0, 1.0])
        mask = np.array([True, True])

        norms = analytic.weighted_error_norms(actual, exact, weights, mask, relative_floor=10.0)

        self.assertAlmostEqual(norms["l2_relative"], math.sqrt(5.0 / 200.0))

    def test_error_norms_reject_invalid_weight_inputs(self):
        analytic = load_module()
        exact = np.array([1.0, 2.0])
        actual = np.array([1.0, 2.0])
        mask = np.array([True, True])

        for weights in (np.array([1.0, -1.0]), np.array([1.0, math.inf]), np.array([0.0, 0.0])):
            with self.subTest(weights=weights):
                with self.assertRaises(ValueError):
                    analytic.weighted_error_norms(actual, exact, weights, mask)

    def test_error_norms_reject_negative_relative_floor(self):
        analytic = load_module()
        exact = np.array([1.0, 2.0])
        actual = np.array([1.0, 2.0])
        weights = np.array([1.0, 1.0])
        mask = np.array([True, True])

        with self.assertRaises(ValueError):
            analytic.weighted_error_norms(actual, exact, weights, mask, relative_floor=-1.0)


if __name__ == "__main__":
    unittest.main()
