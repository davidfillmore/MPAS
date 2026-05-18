#!/usr/bin/env python3
"""Tests for the Tier B idealized electrostatic figure helper."""

import importlib.util
import pathlib
import tempfile
import unittest

import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/run_tier_B_idealized.py"
)


def load_script():
    spec = importlib.util.spec_from_file_location("run_tier_B_idealized", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TierBIdealizedFigureTests(unittest.TestCase):
    def test_default_paths_target_tier_b_outputs(self):
        script = load_script()

        args = script.parse_args([])

        self.assertEqual(
            args.point_output,
            pathlib.Path("~/Data/MPAS/poisson_free_space_charge/gaussian_monopole/run/output.nc"),
        )
        self.assertEqual(
            args.tripole_output,
            pathlib.Path("~/Data/MPAS/poisson_tripole_supercell/run/output.nc"),
        )
        self.assertEqual(args.point_figure.name, "tier_B1_point_charge.pdf")
        self.assertEqual(args.tripole_figure.name, "tier_B2_tripole.pdf")

    def test_coulomb_field_strength_is_inverse_square(self):
        script = load_script()

        e1 = script.coulomb_field_strength(np.array([1000.0]), charge=2.0)[0]
        e2 = script.coulomb_field_strength(np.array([2000.0]), charge=2.0)[0]

        self.assertAlmostEqual(e1 / e2, 4.0)

    def test_radial_profile_recovers_exact_far_field_slope(self):
        script = load_script()
        radius = np.array([6000.0, 9000.0, 14000.0, 22000.0, 35000.0, 55000.0])
        e_mag = script.coulomb_field_strength(radius, charge=3.0)
        fields = {
            "x": radius,
            "y": np.zeros_like(radius),
            "zmid": np.zeros((1, radius.size)),
            "E_mag": e_mag[np.newaxis, :],
        }

        profile = script.point_charge_radial_profile(
            fields,
            center=(0.0, 0.0, 0.0),
            charge=3.0,
            sigma=1000.0,
            min_radius=5000.0,
            max_radius=60000.0,
            bins=6,
        )

        self.assertGreaterEqual(profile["radius"].size, 4)
        self.assertAlmostEqual(profile["far_field_slope"], -2.0, places=6)
        self.assertTrue(np.allclose(profile["e_median"], profile["e_coulomb"]))

    def test_tripole_metrics_report_peak_location(self):
        script = load_script()
        fields = {
            "x": np.array([0.0, 1000.0]),
            "y": np.array([0.0, 0.0]),
            "zmid": np.array([[1000.0, 1000.0], [4000.0, 4000.0]]),
            "E_mag": np.array([[2.0, 3.0], [4.0, 12.0]]),
        }

        metrics = script.tripole_metrics(fields)

        self.assertEqual(metrics["peak_e"], 12.0)
        self.assertEqual(metrics["peak_z"], 4000.0)
        self.assertEqual(metrics["peak_x"], 1000.0)

    def test_write_summary_json_serializes_numpy_values(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            summary = pathlib.Path(tmp) / "summary.json"

            script.write_summary(
                summary,
                point={"far_field_slope": np.float64(-2.0)},
                tripole={"peak_e": np.float64(120000.0)},
            )

            text = summary.read_text()
            self.assertIn('"far_field_slope": -2.0', text)
            self.assertIn('"peak_e": 120000.0', text)


if __name__ == "__main__":
    unittest.main()
