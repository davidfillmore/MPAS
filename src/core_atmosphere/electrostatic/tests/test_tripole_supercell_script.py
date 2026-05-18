#!/usr/bin/env python3
"""Unit tests for the synthetic tripole supercell runner helpers."""

import importlib.util
import pathlib
import tempfile
import unittest

import netCDF4 as nc
import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/run_tripole_supercell.py"
)


def load_script():
    spec = importlib.util.spec_from_file_location("run_tripole_supercell", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TripoleSupercellScriptTests(unittest.TestCase):
    def test_configure_namelist_selects_zero_duration_tripole(self):
        script = load_script()

        namelist = """&nhyd_model
    config_run_duration = '00:03:00'
/
"""
        updated = script.configure_namelist_text(namelist)
        updated_again = script.configure_namelist_text(updated)

        self.assertEqual(updated, updated_again)
        self.assertIn("config_run_duration = '00_00:00:00'", updated)
        self.assertIn("config_electrostatic_enable = .true.", updated)
        self.assertIn("config_electrostatic_solve_at_init = .true.", updated)
        self.assertIn("config_electrostatic_source = 'tripole'", updated)
        self.assertIn("config_poisson_tol = 1.0e-10", updated)

    def test_ensure_output_stream_list_adds_plot_fields(self):
        script = load_script()

        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "stream_list.atmosphere.output"
            path.write_text("theta\nphi\n")
            script.ensure_output_stream_list(path)
            script.ensure_output_stream_list(path)
            fields = [line.strip() for line in path.read_text().splitlines() if line.strip()]

        for field in ("xCell", "yCell", "zgrid", "rho_charge", "phi", "E_vector"):
            self.assertIn(field, fields)
            self.assertEqual(fields.count(field), 1)

    def test_read_plot_fields_handles_e_vector_dimension_order(self):
        script = load_script()

        with tempfile.TemporaryDirectory() as tmp:
            output = pathlib.Path(tmp) / "output.nc"
            write_synthetic_output(output)
            fields = script.read_plot_fields(output)

        self.assertEqual(fields["rho_charge"].shape, (3, 4))
        self.assertEqual(fields["phi"].shape, (3, 4))
        self.assertEqual(fields["E_mag"].shape, (3, 4))
        self.assertTrue(np.allclose(fields["E_mag"], 3.0))
        self.assertEqual(fields["mid_level"], 1)
        self.assertEqual(fields["cross_section_cell_indices"].tolist(), [0, 2])
        self.assertEqual(fields["E_x"].shape, (3, 4))
        self.assertEqual(fields["E_y"].shape, (3, 4))
        self.assertEqual(fields["E_z"].shape, (3, 4))
        self.assertTrue(np.allclose(fields["E_x"], 1.0))
        self.assertTrue(np.allclose(fields["E_y"], 2.0))
        self.assertTrue(np.allclose(fields["E_z"], 2.0))

    def test_stream_grid_from_samples_returns_regular_vectors(self):
        script = load_script()

        x = np.array([0.0, 1.0, 0.0, 1.0])
        y = np.array([0.0, 0.0, 1.0, 1.0])
        u = np.ones(4)
        v = np.full(4, 2.0)

        grid = script.stream_grid_from_samples(x, y, u, v, nx=5, ny=4)

        self.assertIsNotNone(grid)
        grid_x, grid_y, grid_u, grid_v = grid
        self.assertEqual(grid_x.shape, (4, 5))
        self.assertEqual(grid_y.shape, (4, 5))
        self.assertEqual(grid_u.shape, (4, 5))
        self.assertEqual(grid_v.shape, (4, 5))
        self.assertTrue(np.allclose(grid_u, 1.0))
        self.assertTrue(np.allclose(grid_v, 2.0))

    def test_scalar_grid_from_samples_returns_regular_values(self):
        script = load_script()

        x = np.array([0.0, 1.0, 0.0, 1.0])
        y = np.array([0.0, 0.0, 1.0, 1.0])
        values = np.full(4, 5.0)

        grid = script.scalar_grid_from_samples(x, y, values, nx=6, ny=5)

        self.assertIsNotNone(grid)
        grid_x, grid_y, grid_values = grid
        self.assertEqual(grid_x.shape, (5, 6))
        self.assertEqual(grid_y.shape, (5, 6))
        self.assertEqual(grid_values.shape, (5, 6))
        self.assertTrue(np.allclose(grid_values, 5.0))

    def test_plot_tripole_output_writes_png(self):
        script = load_script()

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output = tmp_path / "output.nc"
            plot = tmp_path / "tripole.png"
            write_synthetic_output(output)

            result = script.plot_tripole_output(output, plot)

            self.assertEqual(result, plot)
            self.assertTrue(plot.exists())
            self.assertGreater(plot.stat().st_size, 0)


def write_synthetic_output(path):
    x = np.array([0.0, 1.0, 0.0, 1.0])
    y = np.array([0.0, 0.0, 1.0, 1.0])
    zgrid = np.array(
        [
            [0.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0, 2.0],
            [3.0, 3.0, 3.0, 3.0],
        ]
    )
    rho = np.arange(12, dtype=float).reshape(3, 4)
    phi = rho + 10.0
    e_vector = np.zeros((3, 3, 4))
    e_vector[0, :, :] = 1.0
    e_vector[1, :, :] = 2.0
    e_vector[2, :, :] = 2.0

    with nc.Dataset(path, "w") as ds:
        ds.createDimension("Time", 1)
        ds.createDimension("nCells", 4)
        ds.createDimension("nVertLevels", 3)
        ds.createDimension("nVertLevelsP1", 4)
        ds.createDimension("R3", 3)
        ds.createVariable("xCell", "f8", ("nCells",))[:] = x
        ds.createVariable("yCell", "f8", ("nCells",))[:] = y
        ds.createVariable("zgrid", "f8", ("nVertLevelsP1", "nCells"))[:] = zgrid
        ds.createVariable("rho_charge", "f8", ("Time", "nVertLevels", "nCells"))[
            0, :, :
        ] = rho
        ds.createVariable("phi", "f8", ("Time", "nVertLevels", "nCells"))[0, :, :] = phi
        ds.createVariable("E_vector", "f8", ("Time", "R3", "nVertLevels", "nCells"))[
            0, :, :, :
        ] = e_vector


if __name__ == "__main__":
    unittest.main()
