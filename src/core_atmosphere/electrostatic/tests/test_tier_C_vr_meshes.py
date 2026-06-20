#!/usr/bin/env python3
"""Tests for Task 2/2.5: Tier C VR mesh helpers.

Run with:
    ~/miniconda3/envs/mpas/bin/python -m unittest \
        src.core_atmosphere.electrostatic.tests.test_tier_C_vr_meshes -v
"""

import importlib.util
import pathlib
import unittest
import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[4]
SETUP = ROOT / "src/core_atmosphere/electrostatic/scripts/setup_tier_C_vr_meshes.py"


def load():
    spec = importlib.util.spec_from_file_location("setup_tier_C_vr_meshes", SETUP)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class DensityTests(unittest.TestCase):
    def test_great_circle_known_pairs(self):
        m = load()
        self.assertAlmostEqual(float(m.great_circle_deg(0.0, 0.0, 0.0, 0.0)), 0.0, places=6)
        self.assertAlmostEqual(float(m.great_circle_deg(0.0, 90.0, 0.0, 0.0)), 90.0, places=4)
        self.assertAlmostEqual(float(m.great_circle_deg(90.0, 0.0, 0.0, 0.0)), 90.0, places=4)

    def test_smoothstep_endpoints_and_midpoint(self):
        m = load()
        self.assertEqual(float(m.smoothstep(np.array(5.0), 10.0, 30.0)), 0.0)
        self.assertEqual(float(m.smoothstep(np.array(35.0), 10.0, 30.0)), 1.0)
        self.assertAlmostEqual(float(m.smoothstep(np.array(20.0), 10.0, 30.0)), 0.5, places=6)

    def test_vr_cell_width_asymptotes_and_monotone(self):
        m = load()
        kw = dict(fine_km=120.0, coarse_km=480.0, lat0=0.0, lon0=0.0, r_in_deg=14.0, r_out_deg=34.0)
        lon2d, lat2d = np.meshgrid(np.linspace(-180, 180, 361), np.linspace(-90, 90, 181))
        cw = m.vr_cell_width(lat2d, lon2d, **kw)
        # center is fine, antipode is coarse
        self.assertAlmostEqual(float(m.vr_cell_width(np.array(0.0), np.array(0.0), **kw)), 120.0, places=6)
        self.assertAlmostEqual(float(m.vr_cell_width(np.array(0.0), np.array(180.0), **kw)), 480.0, places=6)
        self.assertGreaterEqual(cw.min(), 120.0 - 1e-9)
        self.assertLessEqual(cw.max(), 480.0 + 1e-9)
        # monotone non-decreasing along a meridian away from center
        d = np.linspace(0, 60, 200)
        prof = m.vr_cell_width(d, np.zeros_like(d), **kw)
        self.assertTrue(np.all(np.diff(prof) >= -1e-12))

    def test_build_cell_width_grid_shapes_and_range(self):
        m = load()
        cw, lon, lat = m.build_cell_width_grid(
            fine_km=120.0, coarse_km=480.0, lat0=0.0, lon0=0.0,
            r_in_deg=14.0, r_out_deg=34.0, dlat_deg=1.0)
        self.assertEqual(cw.shape, (lat.size, lon.size))
        self.assertGreater(cw.max() - cw.min(), 300.0)  # patch resolved on the grid

    def test_find_pentagons_on_synthetic(self):
        import netCDF4, tempfile, os
        m = load()
        path = os.path.join(tempfile.mkdtemp(), "grid.nc")
        with netCDF4.Dataset(path, "w") as ds:
            ds.createDimension("nCells", 5)
            ne = ds.createVariable("nEdgesOnCell", "i4", ("nCells",))
            la = ds.createVariable("latCell", "f8", ("nCells",))
            lo = ds.createVariable("lonCell", "f8", ("nCells",))
            ne[:] = [6, 5, 6, 5, 7]
            la[:] = np.radians([0.0, 10.0, 0.0, -10.0, 0.0])
            lo[:] = np.radians([0.0, 20.0, 0.0, 40.0, 0.0])
        idx, lat, lon = m.find_pentagons(path)
        self.assertEqual(list(idx), [1, 3])
        np.testing.assert_allclose(lat, [10.0, -10.0], atol=1e-6)
        np.testing.assert_allclose(lon, [20.0, 40.0], atol=1e-6)


if __name__ == "__main__":
    unittest.main()
