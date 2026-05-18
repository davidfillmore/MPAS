#!/usr/bin/env python3
"""Tests for the free-space Gaussian charge benchmark runner."""

import contextlib
import importlib.util
import io
import pathlib
import tempfile
import unittest

import netCDF4 as nc
import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/run_free_space_charge_benchmarks.py"
)


def load_script():
    spec = importlib.util.spec_from_file_location("run_free_space_charge_benchmarks", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FreeSpaceChargeRunnerTests(unittest.TestCase):
    def test_configure_namelist_sets_gaussian_parameters(self):
        script = load_script()
        namelist = """&nhyd_model
    config_run_duration = '00:03:00'
/
"""
        updated = script.configure_namelist_text(
            namelist,
            source="gaussian_dipole_y",
            charge=12.5,
            sigma=2500.0,
            separation=9000.0,
            poisson_tol=1.0e-10,
            poisson_max_iter=6000,
        )

        self.assertIn("config_run_duration = '00_00:00:00'", updated)
        self.assertIn("config_electrostatic_source = 'gaussian_dipole_y'", updated)
        self.assertIn("config_electrostatic_source_charge = 12.5", updated)
        self.assertIn("config_electrostatic_source_sigma = 2500.0", updated)
        self.assertIn("config_electrostatic_dipole_separation = 9000.0", updated)
        self.assertIn("config_poisson_tol = 1e-10", updated)
        self.assertIn("config_poisson_max_iter = 6000", updated)

    def test_default_paths_follow_source_name(self):
        script = load_script()

        args = script.parse_args(["--source", "gaussian_monopole"])

        self.assertEqual(
            args.run_dir,
            pathlib.Path("~/Data/MPAS/poisson_free_space_charge/gaussian_monopole/run"),
        )
        self.assertEqual(args.summary.name, "gaussian_monopole_summary.json")
        self.assertEqual(args.plot.name, "gaussian_monopole_diagnostics.png")

    def test_prepare_run_dir_requires_seed_inputs(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            with self.assertRaises(FileNotFoundError):
                script.prepare_run_dir(
                    template_run_dir=tmp_path / "missing_template",
                    run_dir=tmp_path / "run",
                    model=tmp_path / "atmosphere_model",
                    ranks=8,
                    source="gaussian_dipole_y",
                    charge=20.0,
                    sigma=2000.0,
                    separation=8000.0,
                    poisson_tol=1.0e-10,
                    poisson_max_iter=5000,
                )

    def test_prepare_run_dir_seeds_and_configures_inputs(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            template = tmp_path / "template"
            template.mkdir()
            model = tmp_path / "atmosphere_model"
            model.write_text("#!/bin/sh\n")
            (template / "namelist.atmosphere").write_text(
                """&nhyd_model
    config_run_duration = '00:03:00'
    config_block_decomp_file_prefix = 'graph.info.part.'
/
"""
            )
            (template / "streams.atmosphere").write_text(
                """<streams>
  <stream name="output" type="output"/>
</streams>
"""
            )
            expected_partition = template / "graph.info.part.8"
            expected_partition.write_text("fake partition\n")

            run_dir, partition = script.prepare_run_dir(
                template_run_dir=template,
                run_dir=tmp_path / "run",
                model=model,
                ranks=8,
                source="gaussian_dipole_z",
                charge=15.0,
                sigma=3000.0,
                separation=10000.0,
                poisson_tol=1.23e-10,
                poisson_max_iter=7000,
            )

            self.assertEqual(run_dir, (tmp_path / "run").resolve())
            self.assertTrue(partition.exists())
            self.assertEqual(partition.resolve(), expected_partition.resolve())
            self.assertTrue((run_dir / "atmosphere_model").is_symlink())

            namelist = (run_dir / "namelist.atmosphere").read_text()
            self.assertIn("config_electrostatic_source = 'gaussian_dipole_z'", namelist)
            self.assertIn("config_electrostatic_source_charge = 15.0", namelist)
            self.assertIn("config_electrostatic_source_sigma = 3000.0", namelist)
            self.assertIn("config_electrostatic_dipole_separation = 10000.0", namelist)
            self.assertIn("config_poisson_tol = 1.23e-10", namelist)
            self.assertIn("config_poisson_max_iter = 7000", namelist)

            output_fields = (run_dir / "stream_list.atmosphere.output").read_text().splitlines()
            self.assertIn("rho_charge", output_fields)
            self.assertIn("phi", output_fields)
            self.assertIn("E_vector", output_fields)

    def test_task5_module_has_executable_main(self):
        script = load_script()

        self.assertTrue(hasattr(script, "main"))

    def test_analyze_output_writes_summary_for_synthetic_exact_field(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output = tmp_path / "output.nc"
            summary = tmp_path / "summary.json"
            plot = tmp_path / "plot.png"
            write_synthetic_gaussian_output(output, script)
            self.assertEqual(
                script.analytic_center(script.read_fields(output)),
                (5000.0, 5000.0, 5000.0),
            )

            result = script.analyze_output(
                output,
                summary,
                plot,
                source="gaussian_monopole",
                charge=2.0,
                sigma=1000.0,
                separation=4000.0,
                boundary_margin=0.0,
                core_radius=0.0,
                e_relative_floor=0.0,
            )

            self.assertEqual(result["source"], "gaussian_monopole")
            self.assertLess(result["phi_l2_relative"], 1.0e-12)
            self.assertLess(result["e_l2_relative"], 1.0e-12)
            self.assertTrue(summary.exists())
            self.assertTrue(plot.exists())

    def test_analyze_output_default_margins_keep_vertical_comparison_points(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output = tmp_path / "output.nc"
            summary = tmp_path / "summary.json"
            plot = tmp_path / "plot.png"
            args = script.parse_args(["--run-dir", str(tmp_path)])
            write_synthetic_gaussian_output(
                output,
                script,
                source=args.source,
                charge=args.charge,
                sigma=args.sigma,
                separation=args.separation,
                x=np.array([0.0, 25000.0, 50000.0, 25000.0, 25000.0, 25000.0]),
                y=np.array([0.0, 13000.0, 50000.0, 37000.0, 25000.0, 21000.0]),
                zgrid=np.array(
                    [
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        [20000.0, 20000.0, 20000.0, 20000.0, 20000.0, 20000.0],
                    ]
                ),
            )

            result = script.analyze_output(
                output,
                summary,
                plot,
                source=args.source,
                charge=args.charge,
                sigma=args.sigma,
                separation=args.separation,
                boundary_margin=args.boundary_margin,
                vertical_boundary_margin=args.vertical_boundary_margin,
                core_radius=args.core_radius,
                e_relative_floor=args.e_relative_floor,
            )

            self.assertGreater(result["mask_count"], 0)
            self.assertGreater(result["comparison_points"], 0)
            self.assertTrue(np.isfinite(result["phi_l2_relative"]))
            self.assertTrue(np.isfinite(result["e_l2_relative"]))
            self.assertEqual(result["plot_z"], 10000.0)

    def test_analysis_only_main_skips_prepare_run_dir(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            run_dir = tmp_path / "run"
            run_dir.mkdir()
            output = run_dir / "output.nc"
            summary = tmp_path / "summary.json"
            plot = tmp_path / "plot.png"
            write_synthetic_gaussian_output(output, script)

            with contextlib.redirect_stdout(io.StringIO()):
                status = script.main(
                    [
                        "--analysis-only",
                        "--run-dir",
                        str(run_dir),
                        "--template-run-dir",
                        str(tmp_path / "missing-template"),
                        "--model",
                        str(tmp_path / "missing-model"),
                        "--source",
                        "gaussian_monopole",
                        "--charge",
                        "2.0",
                        "--sigma",
                        "1000.0",
                        "--separation",
                        "4000.0",
                        "--boundary-margin",
                        "0.0",
                        "--core-radius",
                        "0.0",
                        "--summary",
                        str(summary),
                        "--plot",
                        str(plot),
                    ]
                )

            self.assertEqual(status, 0)
            self.assertTrue(summary.exists())
            self.assertTrue(plot.exists())


def write_synthetic_gaussian_output(
    path,
    script,
    *,
    source="gaussian_monopole",
    charge=2.0,
    sigma=1000.0,
    separation=4000.0,
    x=None,
    y=None,
    zgrid=None,
):
    analytic = script.analytic
    if x is None:
        x = np.array([3000.0, 5000.0, 7000.0, 5000.0])
    if y is None:
        y = np.array([3000.0, 5000.0, 7000.0, 5000.0])
    if zgrid is None:
        zgrid = np.array(
            [
                [4000.0, 4000.0, 4000.0, 4000.0],
                [6000.0, 6000.0, 6000.0, 6000.0],
            ]
        )
    center = (
        0.5 * (float(np.min(x)) + float(np.max(x))),
        0.5 * (float(np.min(y)) + float(np.max(y))),
        0.5 * (float(np.min(zgrid)) + float(np.max(zgrid))),
    )
    zmid = 0.5 * (zgrid[0, :] + zgrid[1, :])
    field = analytic.evaluate_gaussian_source(
        source,
        x,
        y,
        zmid,
        center=center,
        charge=charge,
        sigma=sigma,
        separation=separation,
    )
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("Time", 1)
        ds.createDimension("nCells", x.size)
        ds.createDimension("nVertLevels", 1)
        ds.createDimension("nVertLevelsP1", 2)
        ds.createDimension("R3", 3)
        ds.createVariable("xCell", "f8", ("nCells",))[:] = x
        ds.createVariable("yCell", "f8", ("nCells",))[:] = y
        ds.createVariable("areaCell", "f8", ("nCells",))[:] = np.ones(x.size)
        ds.createVariable("zgrid", "f8", ("nCells", "nVertLevelsP1"))[:, :] = zgrid.T
        ds.createVariable("rho_charge", "f8", ("Time", "nCells", "nVertLevels"))[0, :, 0] = field.rho
        ds.createVariable("phi", "f8", ("Time", "nCells", "nVertLevels"))[0, :, 0] = field.phi
        e_vector = np.zeros((1, x.size, 1, 3))
        e_vector[0, :, 0, 0] = field.ex
        e_vector[0, :, 0, 1] = field.ey
        e_vector[0, :, 0, 2] = field.ez
        ds.createVariable("E_vector", "f8", ("Time", "nCells", "nVertLevels", "R3"))[:, :, :, :] = e_vector
        ds.createVariable("cg_residual_final", "f8", ("Time",))[:] = [1.0e-12]


if __name__ == "__main__":
    unittest.main()
