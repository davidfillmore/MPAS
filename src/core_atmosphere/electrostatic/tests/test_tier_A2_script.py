#!/usr/bin/env python3
"""Unit tests for the Tier A.2 spherical MMS runner helpers."""

import importlib.util
import pathlib
import tempfile
import unittest

import netCDF4 as nc
import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/run_tier_A2_sphere_mms.py"
)
SETUP_SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/setup_tier_A2_sphere_meshes.py"
)


def load_script():
    spec = importlib.util.spec_from_file_location("run_tier_A2_sphere_mms", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_setup_script():
    spec = importlib.util.spec_from_file_location(
        "setup_tier_A2_sphere_meshes", SETUP_SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TierA2ScriptTests(unittest.TestCase):
    def test_spherical_mms_vertical_basis_matches_boundary_conditions(self):
        script = load_script()

        lat = np.array([0.25])
        lon = np.array([0.5])
        z_top = 20000.0

        self.assertLess(abs(script.phi_exact(lat, lon, np.array([[0.0]]), z_top)), 1.0e-14)

        h = 1.0e-3
        top_deriv = (
            script.phi_exact(lat, lon, np.array([[z_top + h]]), z_top)
            - script.phi_exact(lat, lon, np.array([[z_top - h]]), z_top)
        ) / (2.0 * h)
        self.assertLess(abs(top_deriv.item()), 1.0e-9)

    def test_compute_errors_is_exact_for_synthetic_spherical_output(self):
        script = load_script()

        with tempfile.TemporaryDirectory() as tmp:
            output = pathlib.Path(tmp) / "output.nc"
            lat = np.array([0.0, 0.25, -0.25])
            lon = np.array([0.0, 0.5, 1.0])
            zgrid = np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10000.0, 10000.0, 10000.0],
                    [20000.0, 20000.0, 20000.0],
                ]
            )
            zmid = 0.5 * (zgrid[:-1, :] + zgrid[1:, :])
            phi = script.phi_exact(lat[None, :], lon[None, :], zmid, 20000.0)

            with nc.Dataset(output, "w") as ds:
                ds.createDimension("Time", 1)
                ds.createDimension("nCells", 3)
                ds.createDimension("nVertLevels", 2)
                ds.createDimension("nVertLevelsP1", 3)
                ds.createVariable("phi", "f8", ("Time", "nVertLevels", "nCells"))[
                    0, :, :
                ] = phi
                ds.createVariable("latCell", "f8", ("nCells",))[:] = lat
                ds.createVariable("lonCell", "f8", ("nCells",))[:] = lon
                ds.createVariable("zgrid", "f8", ("nVertLevelsP1", "nCells"))[:] = zgrid
                ds.createVariable("areaCell", "f8", ("nCells",))[:] = 1.0
                ds.createVariable("cg_iter_count", "i4", ("Time",))[:] = 5
                ds.createVariable("cg_residual_initial", "f8", ("Time",))[:] = 1.0
                ds.createVariable("cg_residual_final", "f8", ("Time",))[:] = 1.0e-10

            result = script.compute_errors(output)

        self.assertLess(result["l2"], 1.0e-14)
        self.assertLess(result["linf"], 1.0e-14)
        self.assertEqual(result["cg_iter_count"], 5)
        self.assertEqual(result["nCells"], 3)

    def test_setup_script_mesh_labels_and_bundle_inputs(self):
        setup = load_setup_script()

        self.assertEqual(setup.mesh_spacing_km("240km"), 240.0)
        self.assertEqual(setup.mesh_spacing_km("3750m"), 3.75)
        self.assertIn("config_nvertlevels = 16", setup.init_namelist_text(16, 30000.0))
        self.assertIn("mms_sphere_horizontal", setup.atmosphere_namelist_text())
        self.assertIn("latCell", setup.output_stream_list_text().splitlines())

    def test_runner_accepts_strict_poisson_tolerance_options(self):
        script = load_script()

        args = script.parse_args(
            [
                "--run-root",
                "/tmp/tier_a2",
                "--mesh-list",
                "960km",
                "--ranks",
                "8",
                "--model",
                "/tmp/atmosphere_model",
                "--poisson-tol",
                "1e-13",
                "--poisson-max-iter",
                "9000",
            ]
        )

        self.assertEqual(args.poisson_tol, 1.0e-13)
        self.assertEqual(args.poisson_max_iter, 9000)


if __name__ == "__main__":
    unittest.main()
