#!/usr/bin/env python3
"""Tests for the terrain charge-coupled supercell runner."""

import importlib.util
import pathlib
import tempfile
import unittest

import netCDF4 as nc
import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/run_terrain_charge_coupled_supercell.py"
)


def load_script():
    spec = importlib.util.spec_from_file_location(
        "run_terrain_charge_coupled_supercell", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TerrainChargeCoupledSupercellScriptTests(unittest.TestCase):
    def test_default_paths_target_new_terrain_root(self):
        script = load_script()

        args = script.parse_args([])

        self.assertEqual(
            args.run_dir,
            pathlib.Path(
                "~/Data/MPAS/poisson_charge_coupled_supercell_terrain_h1000/run"
            ),
        )
        self.assertEqual(args.hill_height, 1000.0)
        self.assertEqual(args.hill_half_width_km, 20.0)

    def test_cosine_bell_terrain_is_centered_nonnegative_and_compact(self):
        script = load_script()
        x = np.array([0.0, 10.0, 20.0, 10.0])
        y = np.array([10.0, 25.0, 10.0, 40.0])

        terrain = script.cosine_bell_terrain(
            x, y, hill_height=1000.0, half_width_m=12.0
        )

        self.assertAlmostEqual(float(np.max(terrain)), 1000.0)
        self.assertTrue(np.all(terrain >= 0.0))
        self.assertEqual(float(terrain[-1]), 0.0)

    def test_rewrite_supercell_init_imposes_monotone_terrain_columns(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            init_nc = pathlib.Path(tmp) / "supercell_init.nc"
            write_synthetic_init(init_nc)

            stats = script.rewrite_supercell_init_terrain(
                init_nc, hill_height=1000.0, hill_half_width_m=20.0
            )

            with nc.Dataset(init_nc) as dataset:
                zgrid = np.asarray(dataset.variables["zgrid"][:])
                ter = np.asarray(dataset.variables["ter"][:])

        self.assertEqual(stats["nCells"], 5)
        self.assertAlmostEqual(stats["ztop"], 20000.0)
        self.assertAlmostEqual(float(np.max(ter)), 1000.0)
        self.assertTrue(np.allclose(zgrid[:, 0], ter))
        self.assertTrue(np.allclose(zgrid[:, -1], 20000.0))
        self.assertTrue(np.all(np.diff(zgrid, axis=1) > 0.0))

    def test_rewrite_supercell_init_supports_transposed_zgrid(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            init_nc = pathlib.Path(tmp) / "supercell_init.nc"
            write_synthetic_init_transposed(init_nc)

            stats = script.rewrite_supercell_init_terrain(
                init_nc, hill_height=1000.0, hill_half_width_m=20.0
            )

            with nc.Dataset(init_nc) as dataset:
                zgrid = np.asarray(dataset.variables["zgrid"][:])
                ter = np.asarray(dataset.variables["ter"][:])

        self.assertEqual(stats["nCells"], 5)
        self.assertAlmostEqual(float(np.max(ter)), 1000.0)
        self.assertTrue(np.allclose(zgrid[0, :], ter))
        self.assertTrue(np.allclose(zgrid[-1, :], 20000.0))
        self.assertTrue(np.all(np.diff(zgrid, axis=0) > 0.0))

    def test_rewrite_supercell_init_rejects_degenerate_terrain_without_writing(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            init_nc = pathlib.Path(tmp) / "supercell_init.nc"
            write_synthetic_init(init_nc)
            with nc.Dataset(init_nc) as dataset:
                original_zgrid = np.asarray(dataset.variables["zgrid"][:])
                original_ter = np.asarray(dataset.variables["ter"][:])

            with self.assertRaisesRegex(ValueError, "terrain"):
                script.rewrite_supercell_init_terrain(
                    init_nc, hill_height=20000.0, hill_half_width_m=20.0
                )

            with nc.Dataset(init_nc) as dataset:
                zgrid = np.asarray(dataset.variables["zgrid"][:])
                ter = np.asarray(dataset.variables["ter"][:])

        self.assertTrue(np.allclose(zgrid, original_zgrid))
        self.assertTrue(np.allclose(ter, original_ter))
        self.assertTrue(np.all(np.diff(zgrid, axis=1) > 0.0))

    def test_seed_terrain_run_dir_copies_init_file_before_editing(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            template = root / "template"
            run_dir = root / "run"
            template.mkdir()
            write_minimal_template(template)

            script.seed_terrain_run_dir(template, run_dir)

            init_path = run_dir / "supercell_init.nc"
            self.assertTrue(init_path.exists())
            self.assertFalse(init_path.is_symlink())
            self.assertTrue((run_dir / "supercell_grid.nc").exists())

    def test_seed_terrain_run_dir_rejects_same_dir_and_preserves_init(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            template = pathlib.Path(tmp) / "template"
            template.mkdir()
            write_minimal_template(template)
            init_path = template / "supercell_init.nc"
            original = init_path.read_bytes()

            with self.assertRaisesRegex(ValueError, "same directory"):
                script.seed_terrain_run_dir(template, template)

            self.assertTrue(init_path.exists())
            self.assertEqual(init_path.read_bytes(), original)

    def test_adds_terrain_namelist_keys(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = pathlib.Path(tmp)
            namelist = run_dir / "namelist.atmosphere"
            namelist.write_text("&nhyd_model\n/\n")

            script.add_terrain_namelist_keys(run_dir, hill_height=1000.0)

            text = namelist.read_text()

        self.assertIn("config_electrostatic_terrain_mode = 'zgrid'", text)
        self.assertIn("config_electrostatic_hill_height = 1000.000", text)


def write_synthetic_init(path):
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("nCells", 5)
        ds.createDimension("nVertLevelsP1", 4)
        ds.createVariable("xCell", "f8", ("nCells",))[:] = [0.0, 20.0, 40.0, 40.0, 80.0]
        ds.createVariable("yCell", "f8", ("nCells",))[:] = [0.0, 0.0, 40.0, 80.0, 80.0]
        ds.createVariable("zgrid", "f8", ("nCells", "nVertLevelsP1"))[:] = (
            np.linspace(0.0, 20000.0, 4)[None, :]
        )
        ds.createVariable("ter", "f8", ("nCells",))[:] = 0.0


def write_synthetic_init_transposed(path):
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("nCells", 5)
        ds.createDimension("nVertLevelsP1", 4)
        ds.createVariable("xCell", "f8", ("nCells",))[:] = [0.0, 20.0, 40.0, 40.0, 80.0]
        ds.createVariable("yCell", "f8", ("nCells",))[:] = [0.0, 0.0, 40.0, 80.0, 80.0]
        ds.createVariable("zgrid", "f8", ("nVertLevelsP1", "nCells"))[:] = (
            np.linspace(0.0, 20000.0, 4)[:, None]
        )
        ds.createVariable("ter", "f8", ("nCells",))[:] = 0.0


def write_minimal_template(template):
    write_synthetic_init(template / "supercell_init.nc")
    (template / "supercell_grid.nc").write_bytes(b"grid")
    (template / "namelist.atmosphere").write_text("&nhyd_model\n/\n")
    (template / "streams.atmosphere").write_text("<streams />\n")
    (template / "stream_list.atmosphere.output").write_text("xtime\n")
    (template / "supercell.graph.info.part.1").write_text("graph\n")


if __name__ == "__main__":
    unittest.main()
