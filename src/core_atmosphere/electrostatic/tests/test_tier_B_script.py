#!/usr/bin/env python3
"""Tests for the Tier B idealized electrostatic figure helper."""

import importlib.util
import pathlib
import tempfile
import unittest
from unittest import mock

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

    def test_tripole_style_matches_charge_coupled_theme(self):
        script = load_script()

        self.assertEqual(
            script.TRIPOLE_NEGATIVE_CHARGE_COLORS,
            script.charge_style.NEGATIVE_CHARGE_CONTOUR_COLORS,
        )
        self.assertEqual(
            script.TRIPOLE_POSITIVE_CHARGE_COLORS,
            script.charge_style.POSITIVE_CHARGE_CONTOUR_COLORS,
        )
        self.assertEqual(script.TRIPOLE_CHARGE_LINESTYLE, "solid")
        self.assertEqual(script.TRIPOLE_PHI_COLORBAR_LABEL, r"$\phi$ (MV)")
        self.assertEqual(script.TRIPOLE_EMAG_COLORBAR_LABEL, r"$|E|$ (kV m$^{-1}$)")
        self.assertTrue(np.allclose(script.tripole_phi_colormap()(0.5)[:3], (1.0, 1.0, 1.0)))
        self.assertTrue(
            np.allclose(script.tripole_positive_colormap()(0.0)[:3], (1.0, 1.0, 1.0))
        )

    def test_tripole_plot_quantities_use_paper_units(self):
        script = load_script()
        fields = {
            "y": np.array([0.0, 2000.0, 4000.0]),
            "zmid": np.array([[1000.0, 1000.0, 1000.0], [3000.0, 3000.0, 3000.0]]),
            "phi": np.array([[1.0e6, -2.0e6, 3.0e6], [4.0e6, -5.0e6, 6.0e6]]),
            "E_mag": np.array([[1000.0, 2000.0, 3000.0], [4000.0, 5000.0, 6000.0]]),
            "rho": np.array([[1.0e-11, -2.0e-11, 3.0e-11], [4.0e-11, -5.0e-11, 6.0e-11]]),
            "E_y": np.ones((2, 3)),
            "E_z": -np.ones((2, 3)),
        }

        quantities = script.tripole_plot_quantities(fields, np.array([0, 2]))

        self.assertTrue(np.allclose(quantities["y_km"], [[0.0, 4.0], [0.0, 4.0]]))
        self.assertTrue(np.allclose(quantities["z_km"], [[1.0, 1.0], [3.0, 3.0]]))
        self.assertTrue(np.allclose(quantities["phi_mv"], [[1.0, 3.0], [4.0, 6.0]]))
        self.assertTrue(np.allclose(quantities["e_mag_kv_m"], [[1.0, 3.0], [4.0, 6.0]]))
        self.assertTrue(np.allclose(quantities["rho_nc_m3"], [[0.01, 0.03], [0.04, 0.06]]))

    def test_tripole_plot_uses_field_line_overlay(self):
        script = load_script()
        fields = {
            "x": np.array([0.0, 1000.0, 2000.0]),
            "y": np.array([0.0, 1000.0, 2000.0]),
            "zmid": np.array([[1000.0, 1000.0, 1000.0], [3000.0, 3000.0, 3000.0]]),
            "phi": np.array([[1.0e6, -2.0e6, 1.0e6], [2.0e6, -3.0e6, 2.0e6]]),
            "E_mag": np.array([[1000.0, 2500.0, 1200.0], [1400.0, 2800.0, 1600.0]]),
            "rho": np.array([[1.0e-11, -2.0e-11, 1.0e-11], [2.0e-11, -3.0e-11, 2.0e-11]]),
            "E_y": np.ones((2, 3)),
            "E_z": -np.ones((2, 3)),
            "residual": 1.0e-10,
        }
        with tempfile.TemporaryDirectory() as tmp:
            output = pathlib.Path(tmp) / "tripole.pdf"
            with mock.patch.object(script, "read_fields", return_value=fields), mock.patch.object(
                script, "_section_indices", return_value=np.array([0, 1, 2])
            ), mock.patch.object(script, "overlay_field_lines", return_value=True, create=True) as overlay:
                script.plot_tripole(pathlib.Path("unused.nc"), output)

            self.assertTrue(output.exists())
            overlay.assert_called_once()
            self.assertEqual(overlay.call_args.kwargs["nx"], 75)
            self.assertEqual(overlay.call_args.kwargs["ny"], 55)

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
