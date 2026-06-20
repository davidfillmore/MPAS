#!/usr/bin/env python3
"""Unit tests for the Tier C matched-CVT spatial/rate analysis (synthetic inputs).

Run:
    ~/miniconda3/envs/mpas/bin/python -m unittest \
        src.core_atmosphere.electrostatic.tests.test_tier_C_vr_hotspot -v
"""
import importlib.util
import os
import pathlib
import tempfile
import unittest

import numpy as np
import netCDF4

ROOT = pathlib.Path(__file__).resolve().parents[4]
ANA = ROOT / "src/core_atmosphere/electrostatic/scripts/analyze_tier_C_vr_hotspot.py"


def load():
    spec = importlib.util.spec_from_file_location("analyze_tier_C_vr_hotspot", ANA)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class ClassifyTests(unittest.TestCase):
    def test_classify_regions_by_distance(self):
        m = load()
        lat = np.array([0.0, 0.0, 0.0, 0.0])
        lon = np.array([0.0, 10.0, 25.0, 90.0])  # d ~ 0,10,25,90 deg from (0,0)
        region, d = m.classify_regions(lat, lon, lat0=0.0, lon0=0.0,
                                       r_in_deg=14.0, r_out_deg=34.0)
        self.assertEqual(list(region), ["core", "core", "transition", "far"])
        np.testing.assert_allclose(d[:2], [0.0, 10.0], atol=1e-4)


class AggregationTests(unittest.TestCase):
    def test_group_error_density_partition(self):
        m = load()
        err2 = np.array([10.0, 1.0, 1.0, 1.0])
        ref2 = np.ones(4)
        area = np.ones(4)
        labels = np.array(["a", "b", "b", "c"])
        rows = m.group_error_density(labels, err2, ref2, area)
        frac = {r["group"]: r["err2_frac"] for r in rows}
        self.assertAlmostEqual(frac["a"], 10.0 / 13.0)
        self.assertAlmostEqual(sum(r["err2_frac"] for r in rows), 1.0)
        dens = {r["group"]: r["err_density"] for r in rows}
        self.assertAlmostEqual(dens["b"], 1.0)  # (1+1)/(1+1)

    def test_radial_profile_bins_and_density(self):
        m = load()
        d = np.array([5.0, 6.0, 40.0, 41.0])
        err2 = np.array([1.0, 1.0, 4.0, 4.0])
        ref2 = np.ones(4)
        area = np.ones(4)
        out = m.radial_error_profile(d, err2, ref2, area, n_bins=2)
        self.assertEqual(int(out["count"].sum()), 4)
        self.assertGreater(out["err_density"][-1], out["err_density"][0])  # far bin denser

    def test_vr_rate_pure_power_law(self):
        m = load()
        h = np.array([120.0, 60.0, 30.0])
        l2 = (h / 120.0) ** 2 * 1e-6  # exact slope 2
        out = m.vr_rate(h, l2)
        self.assertAlmostEqual(out["finest_two"], 2.0, places=6)
        self.assertAlmostEqual(out["ols"], 2.0, places=6)

    def test_metric_regression_positive_correlation(self):
        m = load()
        metric = np.linspace(1.0, 5.0, 50)
        ref2 = np.ones(50)
        err2 = (0.1 * metric) ** 2  # relative error proportional to metric
        h = np.full(50, 100.0)
        rows = m.metric_regression(err2, ref2, metric, h, h_bands=[(50.0, 200.0)])
        self.assertEqual(len(rows), 1)
        self.assertGreater(rows[0]["spearman"], 0.9)
        self.assertGreater(rows[0]["ols_slope"], 0.0)


class MetricTests(unittest.TestCase):
    def test_distortion_metrics_hand_mesh(self):
        m = load()
        path = os.path.join(tempfile.mkdtemp(), "grid.nc")
        with netCDF4.Dataset(path, "w") as ds:
            ds.createDimension("nCells", 3)
            ds.createDimension("maxEdges", 2)
            ds.createDimension("nEdges", 3)
            ds.createVariable("nEdgesOnCell", "i4", ("nCells",))[:] = [2, 2, 2]
            ds.createVariable("cellsOnCell", "i4", ("nCells", "maxEdges"))[:] = [[2, 3], [1, 3], [1, 2]]
            ds.createVariable("edgesOnCell", "i4", ("nCells", "maxEdges"))[:] = [[1, 3], [1, 2], [2, 3]]
            ds.createVariable("areaCell", "f8", ("nCells",))[:] = [1.0, 4.0, 1.0]  # cell1 is 4x
            ds.createVariable("dvEdge", "f8", ("nEdges",))[:] = [2.0, 1.0, 1.0]
            ds.createVariable("dcEdge", "f8", ("nEdges",))[:] = [1.0, 1.0, 1.0]
        out = m.distortion_metrics(path)
        np.testing.assert_allclose(out["area_ratio"], [4.0, 4.0, 4.0])
        np.testing.assert_allclose(out["well_centred"][:2], [1.0, 1.0])
        self.assertTrue(bool(out["nonhex"].all()))  # all 2-sided -> non-hex


if __name__ == "__main__":
    unittest.main()
