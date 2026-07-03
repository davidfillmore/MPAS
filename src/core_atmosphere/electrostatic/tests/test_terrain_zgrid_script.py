"""Unit tests for run_terrain_zgrid_mms helpers; no model execution."""

import importlib.util
import math
import pathlib
import tempfile
import unittest

import netCDF4 as nc
import numpy as np


SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestTerrainZgridScript(unittest.TestCase):
    def setUp(self):
        self.script = load("run_terrain_zgrid_mms")

    def test_phi_exact_boundary_conditions(self):
        script, h0, ztop, length = self.script, 1000.0, 20000.0, 90000.0
        x, y = 12345.0, 23456.0
        terrain = h0 * math.cos(2.0 * math.pi * x / length) * math.cos(
            2.0 * math.pi * y / length
        )

        self.assertAlmostEqual(
            script.phi_exact_terrain(x, y, terrain, length, length, ztop, h0),
            0.0,
            places=14,
        )

        dz = 0.01
        dphidz = (
            script.phi_exact_terrain(x, y, ztop, length, length, ztop, h0)
            - script.phi_exact_terrain(x, y, ztop - dz, length, length, ztop, h0)
        ) / dz
        self.assertLess(abs(dphidz), 1.0e-6)

    def test_gate_constants(self):
        self.assertEqual(self.script.GATE_MIN_SLOPE, 1.8)
        self.assertEqual(self.script.BUFFER_LAYERS, 2)
        self.assertEqual(
            self.script.parse_sequence(["15km:25", "7.5km:50"]),
            [("15km", 25), ("7.5km", 50)],
        )

    def test_rewrite_zgrid_terrain_columns(self):
        script, h0, ztop, length, nlevels = self.script, 500.0, 8000.0, 40000.0, 3

        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "init.nc"
            with nc.Dataset(path, "w") as dataset:
                dataset.createDimension("nCells", 4)
                dataset.createDimension("nVertLevelsP1", nlevels + 1)
                dataset.x_period = length
                dataset.y_period = length
                zgrid = dataset.createVariable(
                    "zgrid", "f8", ("nCells", "nVertLevelsP1")
                )
                zgrid[:] = 0.0
                xcell = dataset.createVariable("xCell", "f8", ("nCells",))
                ycell = dataset.createVariable("yCell", "f8", ("nCells",))
                xcell[:] = [0.0, 10000.0, 20000.0, 30000.0]
                ycell[:] = [0.0, 5000.0, 10000.0, 15000.0]

            script.rewrite_zgrid_terrain(path, hill_height=h0, ztop=ztop)
            with nc.Dataset(path) as dataset:
                zgrid = np.asarray(dataset.variables["zgrid"][:])
                xcell = np.asarray(dataset.variables["xCell"][:])
                ycell = np.asarray(dataset.variables["yCell"][:])

            terrain = h0 * np.cos(2.0 * np.pi * xcell / length) * np.cos(
                2.0 * np.pi * ycell / length
            )
            zeta = np.linspace(0.0, ztop, nlevels + 1)
            expected = zeta[None, :] + terrain[:, None] * (1.0 - zeta[None, :] / ztop)
            self.assertLess(float(np.abs(zgrid - expected).max()), 1.0e-10)

            script.rewrite_zgrid_terrain(path, hill_height=h0, ztop=ztop)
            with nc.Dataset(path) as dataset:
                zgrid2 = np.asarray(dataset.variables["zgrid"][:])
            self.assertLess(float(np.abs(zgrid2 - expected).max()), 1.0e-10)

    def test_check_residual_rejects_nan_and_excess(self):
        # Within tolerance: no exception.
        self.script.check_residual("ok", 9.9e-9, 1.0e-8)

        # Over tolerance: raises, with the value in the message.
        with self.assertRaisesRegex(RuntimeError, "2.000e-08"):
            self.script.check_residual("high", 2.0e-8, 1.0e-8)

        # NaN (diverged solve, or cg_residual_final absent so
        # scalar_diagnostic returned its NaN default): must also raise.
        with self.assertRaisesRegex(RuntimeError, "nan"):
            self.script.check_residual("diverged", math.nan, 1.0e-8)

    def test_finest_two_slope(self):
        rows = [
            {"mesh": "15km", "h_m": 15000.0, "l2_interior": 4.0e-3},
            {"mesh": "7.5km", "h_m": 7500.0, "l2_interior": 1.0e-3},
        ]
        # log(4e-3/1e-3) / log(15000/7500) = log 4 / log 2 = 2.0 exactly.
        self.assertAlmostEqual(
            self.script.finest_two_slope(rows, "l2_interior"), 2.0, places=12
        )
        # Fewer than two usable rows: NaN, not an exception.
        self.assertTrue(
            math.isnan(self.script.finest_two_slope(rows[:1], "l2_interior"))
        )

    def test_finest_two_slope_vertical_only_sequence_exits(self):
        rows = [
            {"mesh": "15km", "h_m": 15000.0, "l2_interior": 1.0e-3},
            {"mesh": "15km", "h_m": 15000.0, "l2_interior": 2.5e-4},
        ]
        with self.assertRaises(SystemExit):
            self.script.finest_two_slope(rows, "l2_interior")


if __name__ == "__main__":
    unittest.main()
