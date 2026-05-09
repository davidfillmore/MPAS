#!/usr/bin/env python3
"""Unit tests for the Tier A.1 Cartesian MMS runner helpers."""

import importlib.util
import pathlib
import tempfile
import unittest

import netCDF4 as nc
import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/run_tier_A1_cartesian_mms.py"
)


def load_script():
    spec = importlib.util.spec_from_file_location("run_tier_A1_cartesian_mms", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TierA1ScriptTests(unittest.TestCase):
    def test_namelist_and_stream_updates_are_idempotent(self):
        script = load_script()

        namelist = """&nhyd_model
    config_run_duration = '00:03:00'
/
&electrostatic
    config_electrostatic_enable = .false.
/
"""
        updated = script.configure_namelist_text(namelist)
        updated_again = script.configure_namelist_text(updated)

        self.assertEqual(updated, updated_again)
        self.assertIn("config_run_duration = '00_00:00:00'", updated)
        self.assertIn("config_electrostatic_enable = .true.", updated)
        self.assertIn("config_poisson_max_iter = 2000", updated)

        streams = """<streams>
<stream name="output" type="output" filename_template="history.nc">
  <file name="stream_list.atmosphere.output"/>
</stream>
</streams>
"""
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "streams.atmosphere"
            path.write_text(streams)
            script.ensure_streams_netcdf(path)
            script.ensure_streams_netcdf(path)
            text = path.read_text()

        self.assertIn('io_type="netcdf"', text)
        self.assertIn('filename_template="output.nc"', text)
        self.assertIn('output_interval="initial_only"', text)

    def test_compute_errors_is_exact_for_synthetic_output(self):
        script = load_script()

        with tempfile.TemporaryDirectory() as tmp:
            output = pathlib.Path(tmp) / "output.nc"
            x = np.array([0.0, 0.5, 1.0])
            y = np.array([0.0, 0.5, 1.0])
            zgrid = np.array(
                [
                    [0.0, 0.0, 0.0],
                    [0.5, 0.5, 0.5],
                    [1.0, 1.0, 1.0],
                ]
            )
            zmid = 0.5 * (zgrid[:-1, :] + zgrid[1:, :])
            phi = script.phi_exact(x[None, :], y[None, :], zmid, 1.0, 1.0, 1.0)

            with nc.Dataset(output, "w") as ds:
                ds.createDimension("Time", 1)
                ds.createDimension("nCells", 3)
                ds.createDimension("nVertLevels", 2)
                ds.createDimension("nVertLevelsP1", 3)
                ds.createVariable("phi", "f8", ("Time", "nVertLevels", "nCells"))[
                    0, :, :
                ] = phi
                ds.createVariable("xCell", "f8", ("nCells",))[:] = x
                ds.createVariable("yCell", "f8", ("nCells",))[:] = y
                ds.createVariable("zgrid", "f8", ("nVertLevelsP1", "nCells"))[:] = zgrid
                ds.createVariable("areaCell", "f8", ("nCells",))[:] = 1.0
                ds.createDimension("nEdges", 1)
                ds.createDimension("R3", 3)
                ds.createVariable("rho_charge", "f8", ("Time", "nVertLevels", "nCells"))[
                    0, :, :
                ] = 0.0
                ds.createVariable("E_normal", "f8", ("Time", "nVertLevels", "nEdges"))[
                    0, :, :
                ] = 0.0
                ds.createVariable(
                    "E_vector", "f8", ("Time", "R3", "nVertLevels", "nCells")
                )[0, :, :, :] = 0.0
                ds.createVariable("cg_iter_count", "i4", ("Time",))[:] = 4
                ds.createVariable("cg_residual_initial", "f8", ("Time",))[:] = 1.0
                ds.createVariable("cg_residual_final", "f8", ("Time",))[:] = 1.0e-10

            result = script.compute_errors(output)

        self.assertLess(result["l2"], 1.0e-14)
        self.assertLess(result["linf"], 1.0e-14)
        self.assertEqual(result["cg_iter_count"], 4)
        self.assertEqual(result["nCells"], 3)

    def test_compute_errors_uses_periodic_mesh_attributes(self):
        script = load_script()

        with tempfile.TemporaryDirectory() as tmp:
            output = pathlib.Path(tmp) / "output.nc"
            x = np.array([0.0, 1.0, 2.0])
            y = np.array([0.0, 1.0, 2.0])
            zgrid = np.array(
                [
                    [0.0, 0.0, 0.0],
                    [0.5, 0.5, 0.5],
                    [1.0, 1.0, 1.0],
                ]
            )
            zmid = 0.5 * (zgrid[:-1, :] + zgrid[1:, :])
            phi = script.phi_exact(x[None, :], y[None, :], zmid, 4.0, 5.0, 1.0)

            with nc.Dataset(output, "w") as ds:
                ds.setncattr("is_periodic", "YES")
                ds.setncattr("x_period", 4.0)
                ds.setncattr("y_period", 5.0)
                ds.createDimension("Time", 1)
                ds.createDimension("nCells", 3)
                ds.createDimension("nVertLevels", 2)
                ds.createDimension("nVertLevelsP1", 3)
                ds.createDimension("nEdges", 1)
                ds.createDimension("R3", 3)
                ds.createVariable("phi", "f8", ("Time", "nVertLevels", "nCells"))[
                    0, :, :
                ] = phi
                ds.createVariable("xCell", "f8", ("nCells",))[:] = x
                ds.createVariable("yCell", "f8", ("nCells",))[:] = y
                ds.createVariable("zgrid", "f8", ("nVertLevelsP1", "nCells"))[:] = zgrid
                ds.createVariable("areaCell", "f8", ("nCells",))[:] = 1.0
                ds.createVariable("rho_charge", "f8", ("Time", "nVertLevels", "nCells"))[
                    0, :, :
                ] = 0.0
                ds.createVariable("E_normal", "f8", ("Time", "nVertLevels", "nEdges"))[
                    0, :, :
                ] = 0.0
                ds.createVariable(
                    "E_vector", "f8", ("Time", "R3", "nVertLevels", "nCells")
                )[0, :, :, :] = 0.0
                ds.createVariable("cg_iter_count", "i4", ("Time",))[:] = 4
                ds.createVariable("cg_residual_initial", "f8", ("Time",))[:] = 1.0
                ds.createVariable("cg_residual_final", "f8", ("Time",))[:] = 1.0e-10

            result = script.compute_errors(output)

        self.assertLess(result["l2"], 1.0e-14)
        self.assertLess(result["linf"], 1.0e-14)
        self.assertEqual(result["x_period"], 4.0)
        self.assertEqual(result["y_period"], 5.0)

    def test_namelist_source_can_select_horizontal_mms(self):
        script = load_script()

        namelist = """&nhyd_model
/
"""
        updated = script.configure_namelist_text(
            namelist,
            source="mms_cart_horizontal",
        )

        self.assertIn("config_electrostatic_source = 'mms_cart_horizontal'", updated)


if __name__ == "__main__":
    unittest.main()
