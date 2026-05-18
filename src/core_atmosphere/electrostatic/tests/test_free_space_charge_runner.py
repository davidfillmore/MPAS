#!/usr/bin/env python3
"""Tests for the free-space Gaussian charge benchmark runner."""

import importlib.util
import pathlib
import tempfile
import unittest


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

    def test_task4_module_has_no_executable_main(self):
        script = load_script()

        self.assertFalse(hasattr(script, "main"))


if __name__ == "__main__":
    unittest.main()
