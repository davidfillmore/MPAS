#!/usr/bin/env python3
"""Unit tests for terrain charge-coupled supercell figure helpers."""

import importlib.util
import pathlib
import tempfile
import unittest

import netCDF4 as nc
import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/plot_terrain_charge_coupled_supercell.py"
)


def load_script():
    spec = importlib.util.spec_from_file_location(
        "plot_terrain_charge_coupled_supercell", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TerrainChargeCoupledSupercellPlotTests(unittest.TestCase):
    def test_terrain_slice_section_uses_selected_slice_bottom_terrain(self):
        script = load_script()
        fields = synthetic_fields()

        section = script.terrain_slice_section(fields)

        self.assertTrue(np.allclose(section["terrain_y"], [0.0, 1.0]))
        self.assertTrue(np.allclose(section["terrain_z"], [1000.0, 1000.0]))

    def test_terrain_xmean_section_uses_area_weighted_row_terrain(self):
        script = load_script()
        fields = synthetic_fields()

        section = script.terrain_xmean_section(fields, core_half_width_m=2.0)

        self.assertTrue(np.allclose(section["terrain_y"], [0.0, 1.0]))
        self.assertTrue(np.allclose(section["terrain_z"], [750.0, 750.0]))

    def test_default_plot_path_targets_terrain_results_directory(self):
        script = load_script()

        output = pathlib.Path(
            "~/Data/MPAS/poisson_charge_coupled_supercell_terrain_h1000/run/output.nc"
        )
        expected = pathlib.Path(
            "~/Data/MPAS/poisson_charge_coupled_supercell_terrain_h1000/results/"
            "terrain_charge_coupled_slice.png"
        ).expanduser()

        self.assertEqual(script.default_plot_path(output, section_mode="slice"), expected)

    def test_plot_terrain_charge_coupled_output_writes_png(self):
        script = load_script()

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output = tmp_path / "output.nc"
            plot = tmp_path / "terrain_charge_coupled.png"
            write_synthetic_output(output)

            result = script.plot_terrain_charge_coupled_output(output, plot)

            self.assertEqual(result, plot)
            self.assertTrue(plot.exists())
            self.assertGreater(plot.stat().st_size, 0)


def synthetic_fields():
    return {
        "xCell": np.array([0.0, 2.0, 0.0, 2.0]),
        "yCell": np.array([0.0, 0.0, 1.0, 1.0]),
        "areaCell": np.array([1.0, 3.0, 1.0, 3.0]),
        "zgrid": np.array(
            [
                [0.0, 1000.0, 0.0, 1000.0],
                [10000.0, 10500.0, 10000.0, 10500.0],
                [20000.0, 20000.0, 20000.0, 20000.0],
            ]
        ),
        "zMid": np.array(
            [
                [5000.0, 5750.0, 5000.0, 5750.0],
                [15000.0, 15250.0, 15000.0, 15250.0],
            ]
        ),
        "liquid_water_content": np.array([[0.0, 4.0, 0.0, 1.0], [0.0, 5.0, 0.0, 1.0]]),
        "rho_charge": np.zeros((2, 4)),
        "phi": np.zeros((2, 4)),
        "E_y": np.zeros((2, 4)),
        "E_z": np.zeros((2, 4)),
    }


def write_synthetic_output(path):
    fields = synthetic_fields()
    rho = np.ones((1, 2, 4), dtype=float)
    qc = np.zeros((1, 2, 4), dtype=float)
    qr = np.zeros((1, 2, 4), dtype=float)
    qc[0, :, :] = fields["liquid_water_content"]
    rho_charge = fields["rho_charge"][None, :, :]
    phi = fields["phi"][None, :, :]
    e_vector = np.zeros((1, 3, 2, 4), dtype=float)
    e_vector[0, 1, :, :] = fields["E_y"]
    e_vector[0, 2, :, :] = fields["E_z"]

    with nc.Dataset(path, "w") as ds:
        ds.createDimension("Time", 1)
        ds.createDimension("StrLen", 64)
        ds.createDimension("nCells", 4)
        ds.createDimension("nVertLevels", 2)
        ds.createDimension("nVertLevelsP1", 3)
        ds.createDimension("R3", 3)

        xtime = np.full((1, 64), b"\x00", dtype="S1")
        text = "0000-01-01_00:10:00"
        xtime[0, : len(text)] = np.frombuffer(text.encode("ascii"), dtype="S1")
        ds.createVariable("xtime", "S1", ("Time", "StrLen"))[:] = xtime
        ds.createVariable("xCell", "f8", ("nCells",))[:] = fields["xCell"]
        ds.createVariable("yCell", "f8", ("nCells",))[:] = fields["yCell"]
        ds.createVariable("areaCell", "f8", ("nCells",))[:] = fields["areaCell"]
        ds.createVariable("zgrid", "f8", ("nVertLevelsP1", "nCells"))[:] = fields["zgrid"]
        ds.createVariable("qc", "f8", ("Time", "nVertLevels", "nCells"))[:] = qc
        ds.createVariable("qr", "f8", ("Time", "nVertLevels", "nCells"))[:] = qr
        ds.createVariable("rho", "f8", ("Time", "nVertLevels", "nCells"))[:] = rho
        ds.createVariable("rho_charge", "f8", ("Time", "nVertLevels", "nCells"))[:] = rho_charge
        ds.createVariable("phi", "f8", ("Time", "nVertLevels", "nCells"))[:] = phi
        ds.createVariable("E_vector", "f8", ("Time", "R3", "nVertLevels", "nCells"))[:] = (
            e_vector
        )


if __name__ == "__main__":
    unittest.main()
