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

    def test_terrain_namelist_enables_zgrid_without_inert_hill_key(self):
        # Finding #17: config_electrostatic_hill_height is read by the
        # driver only for source='mms_terrain'; the supercell runs
        # source='stub' with the hill baked into supercell_init.nc.
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = pathlib.Path(tmp)
            namelist = run_dir / "namelist.atmosphere"
            namelist.write_text("&nhyd_model\n/\n")

            script.add_terrain_namelist_keys(run_dir)

            text = namelist.read_text()

        self.assertIn("config_electrostatic_terrain_mode = 'zgrid'", text)
        self.assertNotIn("config_electrostatic_hill_height", text)

    def test_recompute_zgrid_metrics_matches_hand_values(self):
        # One edge between a flat column and a terrain column; every
        # expected number below is transcribed independently from
        # mpas_init_atm_cases.F (zz :1677, zxu :1685, zb/zb3 :1149-1181).
        script = load_script()
        zgrid = np.array(
            [
                [0.0, 1000.0, 2000.0, 3000.0],    # flat column
                [300.0, 1200.0, 2100.0, 3000.0],  # ter=300, ztop=3000
            ]
        )
        dzw = np.array([1000.0, 1000.0, 1000.0])
        mesh = {
            "cellsOnEdge": np.array([[1, 2]]),
            "cellsOnCell": np.array([[2], [1]]),
            "nEdgesOnCell": np.array([1, 1]),
            "derivTwo": np.zeros((1, 2, 15)),
            "dcEdge": np.array([500.0]),
            "dvEdge": np.array([300.0]),
            "areaCell": np.array([2.0e5, 2.0e5]),
        }

        metrics = script.recompute_zgrid_metrics(zgrid, dzw, mesh, theta_adv_order=3)

        self.assertTrue(np.allclose(metrics["zz"][0], 1.0))
        self.assertTrue(np.allclose(metrics["zz"][1], 1000.0 / 900.0))
        # zg2 - zg1 = [300, 200, 100, 0] -> zxu = [0.5, 0.3, 0.1]
        self.assertTrue(np.allclose(metrics["zxu"][0], [0.5, 0.3, 0.1]))
        # deriv_two = 0: z_edge = midpoint, zb = +/- 0.5*(zg2-zg1)*dv/area
        self.assertTrue(
            np.allclose(metrics["zb"][0, 0], [0.225, 0.15, 0.075, 0.0])
        )
        self.assertTrue(
            np.allclose(metrics["zb"][0, 1], [-0.225, -0.15, -0.075, 0.0])
        )
        self.assertTrue(np.allclose(metrics["zb3"], 0.0))

        # Nonzero deriv_two exercises the neighbor accumulation and zb3.
        mesh["derivTwo"][0, 0, :2] = [2.0e-6, -2.0e-6]   # cell1 side: self, nbr
        mesh["derivTwo"][0, 1, :2] = [1.0e-6, -1.0e-6]   # cell2 side: self, nbr
        metrics = script.recompute_zgrid_metrics(zgrid, dzw, mesh, theta_adv_order=3)
        for k in range(3):
            z1, z2 = zgrid[0, k], zgrid[1, k]
            d2_1 = 2.0e-6 * z1 - 2.0e-6 * z2
            d2_2 = 1.0e-6 * z2 - 1.0e-6 * z1
            z_edge = 0.5 * (z1 + z2) - 500.0**2 * (d2_1 + d2_2) / 12.0
            z_edge3 = -(500.0**2) * (d2_1 - d2_2) / 12.0
            self.assertAlmostEqual(
                metrics["zb"][0, 0, k], (z_edge - z1) * 300.0 / 2.0e5, places=12
            )
            self.assertAlmostEqual(
                metrics["zb"][0, 1, k], (z_edge - z2) * 300.0 / 2.0e5, places=12
            )
            self.assertAlmostEqual(
                metrics["zb3"][0, 0, k], z_edge3 * 300.0 / 2.0e5, places=12
            )

    def test_verify_init_metrics_flat_identity_and_corruption(self):
        # Golden identity: on an UNMODIFIED flat init file the recompute
        # must reproduce every stored metric field (finding #9 guard).
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            init_nc = pathlib.Path(tmp) / "supercell_init.nc"
            write_synthetic_init(init_nc)

            report = script.verify_init_metrics(init_nc)

            self.assertEqual(
                sorted(report),
                sorted(
                    ["zz", "zxu", "zb", "zb3", "dzu", "rdzu",
                     "fzm", "fzp", "cf1", "cf2", "cf3"]
                ),
            )
            self.assertLess(max(report.values()), 1.0e-9)

            with nc.Dataset(init_nc, "r+") as dataset:
                dataset.variables["zz"][0, 0] = 1.001
            with self.assertRaisesRegex(RuntimeError, "zz"):
                script.verify_init_metrics(init_nc)

    def test_rewrite_recomputes_zgrid_derived_metrics(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            init_nc = pathlib.Path(tmp) / "supercell_init.nc"
            write_synthetic_init(init_nc)

            stats = script.rewrite_supercell_init_terrain(
                init_nc, hill_height=1000.0, hill_half_width_m=60.0
            )

            # Self-consistency: stored metrics match a fresh recompute.
            report = script.verify_init_metrics(init_nc)
            self.assertLess(max(report.values()), 1.0e-9)

            with nc.Dataset(init_nc) as dataset:
                zz = np.asarray(dataset.variables["zz"][:])
                zxu = np.asarray(dataset.variables["zxu"][:])
                zb = np.asarray(dataset.variables["zb"][:])
                zb3 = np.asarray(dataset.variables["zb3"][:])

        # Finding #9: the hill must exert a coordinate slope (zxu was
        # identically 0) and compress the metric (zz > 1 under the hill).
        self.assertGreater(float(np.abs(zxu).max()), 0.0)
        self.assertGreater(float(zz.max()), 1.0)
        self.assertGreater(float(np.abs(zb).max()), 0.0)
        self.assertTrue(np.allclose(zb3, 0.0))  # deriv_two == 0 here
        self.assertGreater(stats["zxu_max_abs"], 0.0)
        self.assertGreater(stats["zz_max"], 1.0)


def write_synthetic_init(path):
    """Flat 5-cell chain mesh with self-consistent stored metrics.

    Stored values mirror what init_atmosphere writes for a flat plane:
    zz = 1, zxu = zb = zb3 = 0, uniform-zeta 1-D metrics with level-1
    entries 0 (Fortran loops start at k=2; Registry zero-initializes).
    """
    zeta = np.linspace(0.0, 20000.0, 4)
    dzw = np.diff(zeta)
    d = float(dzw[0])
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("nCells", 5)
        ds.createDimension("nEdges", 4)
        ds.createDimension("nVertLevels", 3)
        ds.createDimension("nVertLevelsP1", 4)
        ds.createDimension("TWO", 2)
        ds.createDimension("maxEdges", 2)
        ds.createDimension("FIFTEEN", 15)
        ds.config_theta_adv_order = 3
        ds.config_interface_projection = "linear_interpolation"
        ds.createVariable("xCell", "f8", ("nCells",))[:] = [0.0, 20.0, 40.0, 40.0, 80.0]
        ds.createVariable("yCell", "f8", ("nCells",))[:] = [0.0, 0.0, 40.0, 80.0, 80.0]
        ds.createVariable("zgrid", "f8", ("nCells", "nVertLevelsP1"))[:] = zeta[None, :]
        ds.createVariable("ter", "f8", ("nCells",))[:] = 0.0
        ds.createVariable("rdzw", "f8", ("nVertLevels",))[:] = 1.0 / dzw
        ds.createVariable("dzu", "f8", ("nVertLevels",))[:] = [0.0, d, d]
        ds.createVariable("rdzu", "f8", ("nVertLevels",))[:] = [0.0, 1.0 / d, 1.0 / d]
        ds.createVariable("fzm", "f8", ("nVertLevels",))[:] = [0.0, 0.5, 0.5]
        ds.createVariable("fzp", "f8", ("nVertLevels",))[:] = [0.0, 0.5, 0.5]
        ds.createVariable("cf1", "f8")[:] = 2.0
        ds.createVariable("cf2", "f8")[:] = -1.5
        ds.createVariable("cf3", "f8")[:] = 0.5
        ds.createVariable("zz", "f8", ("nCells", "nVertLevels"))[:] = 1.0
        ds.createVariable("zxu", "f8", ("nEdges", "nVertLevels"))[:] = 0.0
        ds.createVariable("zb", "f8", ("nEdges", "TWO", "nVertLevelsP1"))[:] = 0.0
        ds.createVariable("zb3", "f8", ("nEdges", "TWO", "nVertLevelsP1"))[:] = 0.0
        ds.createVariable("deriv_two", "f8", ("nEdges", "TWO", "FIFTEEN"))[:] = 0.0
        ds.createVariable("cellsOnEdge", "i4", ("nEdges", "TWO"))[:] = [
            [1, 2], [2, 3], [3, 4], [4, 5]]
        ds.createVariable("cellsOnCell", "i4", ("nCells", "maxEdges"))[:] = [
            [2, 0], [1, 3], [2, 4], [3, 5], [4, 0]]
        ds.createVariable("nEdgesOnCell", "i4", ("nCells",))[:] = [1, 2, 2, 2, 1]
        ds.createVariable("dcEdge", "f8", ("nEdges",))[:] = 20.0
        ds.createVariable("dvEdge", "f8", ("nEdges",))[:] = 20.0
        ds.createVariable("areaCell", "f8", ("nCells",))[:] = 400.0


def write_synthetic_init_transposed(path):
    """Same flat mesh as write_synthetic_init, level-bearing dims reversed."""
    zeta = np.linspace(0.0, 20000.0, 4)
    dzw = np.diff(zeta)
    d = float(dzw[0])
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("nCells", 5)
        ds.createDimension("nEdges", 4)
        ds.createDimension("nVertLevels", 3)
        ds.createDimension("nVertLevelsP1", 4)
        ds.createDimension("TWO", 2)
        ds.createDimension("maxEdges", 2)
        ds.createDimension("FIFTEEN", 15)
        ds.config_theta_adv_order = 3
        ds.config_interface_projection = "linear_interpolation"
        ds.createVariable("xCell", "f8", ("nCells",))[:] = [0.0, 20.0, 40.0, 40.0, 80.0]
        ds.createVariable("yCell", "f8", ("nCells",))[:] = [0.0, 0.0, 40.0, 80.0, 80.0]
        ds.createVariable("zgrid", "f8", ("nVertLevelsP1", "nCells"))[:] = zeta[:, None]
        ds.createVariable("ter", "f8", ("nCells",))[:] = 0.0
        ds.createVariable("rdzw", "f8", ("nVertLevels",))[:] = 1.0 / dzw
        ds.createVariable("dzu", "f8", ("nVertLevels",))[:] = [0.0, d, d]
        ds.createVariable("rdzu", "f8", ("nVertLevels",))[:] = [0.0, 1.0 / d, 1.0 / d]
        ds.createVariable("fzm", "f8", ("nVertLevels",))[:] = [0.0, 0.5, 0.5]
        ds.createVariable("fzp", "f8", ("nVertLevels",))[:] = [0.0, 0.5, 0.5]
        ds.createVariable("cf1", "f8")[:] = 2.0
        ds.createVariable("cf2", "f8")[:] = -1.5
        ds.createVariable("cf3", "f8")[:] = 0.5
        ds.createVariable("zz", "f8", ("nVertLevels", "nCells"))[:] = 1.0
        ds.createVariable("zxu", "f8", ("nVertLevels", "nEdges"))[:] = 0.0
        ds.createVariable("zb", "f8", ("nVertLevelsP1", "TWO", "nEdges"))[:] = 0.0
        ds.createVariable("zb3", "f8", ("nVertLevelsP1", "TWO", "nEdges"))[:] = 0.0
        ds.createVariable("deriv_two", "f8", ("FIFTEEN", "TWO", "nEdges"))[:] = 0.0
        ds.createVariable("cellsOnEdge", "i4", ("TWO", "nEdges"))[:] = np.transpose(
            [[1, 2], [2, 3], [3, 4], [4, 5]])
        ds.createVariable("cellsOnCell", "i4", ("maxEdges", "nCells"))[:] = np.transpose(
            [[2, 0], [1, 3], [2, 4], [3, 5], [4, 0]])
        ds.createVariable("nEdgesOnCell", "i4", ("nCells",))[:] = [1, 2, 2, 2, 1]
        ds.createVariable("dcEdge", "f8", ("nEdges",))[:] = 20.0
        ds.createVariable("dvEdge", "f8", ("nEdges",))[:] = 20.0
        ds.createVariable("areaCell", "f8", ("nCells",))[:] = 400.0


def write_minimal_template(template):
    write_synthetic_init(template / "supercell_init.nc")
    (template / "supercell_grid.nc").write_bytes(b"grid")
    (template / "namelist.atmosphere").write_text("&nhyd_model\n/\n")
    (template / "streams.atmosphere").write_text("<streams />\n")
    (template / "stream_list.atmosphere.output").write_text("xtime\n")
    (template / "supercell.graph.info.part.1").write_text("graph\n")


if __name__ == "__main__":
    unittest.main()
