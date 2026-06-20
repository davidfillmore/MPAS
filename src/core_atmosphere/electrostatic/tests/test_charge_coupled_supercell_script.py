#!/usr/bin/env python3
"""Tests for the diagnostic charge-coupled supercell runner."""

import importlib.util
import pathlib
import tempfile
import unittest
import xml.etree.ElementTree as ET


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/run_charge_coupled_supercell.py"
)
REGISTRY_PATH = REPO_ROOT / "src/core_atmosphere/Registry.xml"


def load_script():
    spec = importlib.util.spec_from_file_location("run_charge_coupled_supercell", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ChargeCoupledSupercellScriptTests(unittest.TestCase):
    def test_configure_namelist_selects_dynamic_stub_source(self):
        script = load_script()
        namelist = """&nhyd_model
    config_run_duration = '00_00:00:00'
/
"""

        updated = script.configure_namelist_text(
            namelist,
            run_duration="00_00:10:00",
            interval=60.0,
            stub_alpha=2.5e-8,
            stub_beta=1.5e-8,
            poisson_tol=1.0e-10,
            poisson_max_iter=6000,
            solve_at_init=False,
        )

        self.assertIn("config_run_duration = '00_00:10:00'", updated)
        self.assertIn("config_electrostatic_enable = .true.", updated)
        self.assertIn("config_electrostatic_solve_at_init = .false.", updated)
        self.assertIn("config_electrostatic_source = 'stub'", updated)
        self.assertIn("config_electrostatic_interval = 60.0", updated)
        self.assertIn("config_stub_alpha = 2.5e-08", updated)
        self.assertIn("config_stub_beta = 1.5e-08", updated)
        self.assertIn("config_poisson_tol = 1e-10", updated)
        self.assertIn("config_poisson_max_iter = 6000", updated)

    def test_output_stream_list_is_compact_and_includes_coupling_fields(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "stream_list.atmosphere.output"
            path.write_text("theta\n")

            script.ensure_coupled_output_stream_list(path)
            script.ensure_coupled_output_stream_list(path)
            fields = [line.strip() for line in path.read_text().splitlines() if line.strip()]

        for field in (
            "xtime",
            "rho",
            "w",
            "scalars",
            "rho_charge",
            "phi",
            "E_vector",
            "cg_residual_final",
        ):
            self.assertIn(field, fields)
            self.assertEqual(fields.count(field), 1)
        self.assertNotIn("theta", fields)
        self.assertNotIn("qg", fields)
        self.assertNotIn("qi", fields)

    def test_configure_streams_sets_dynamic_output_interval(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "streams.atmosphere"
            path.write_text(
                """<streams>
  <stream name="output" type="output" filename_template="old.nc" output_interval="initial_only">
    <file name="stream_list.atmosphere.output" />
  </stream>
  <stream name="diagnostics" type="output" output_interval="none" />
</streams>
"""
            )

            script.configure_streams(path, "00:10:00")

            text = path.read_text()

        self.assertIn('io_type="netcdf"', text)
        self.assertIn('filename_template="output.nc"', text)
        self.assertIn('filename_interval="none"', text)
        self.assertIn('output_interval="00:10:00"', text)

    def test_default_paths_target_charge_coupled_run_root(self):
        script = load_script()

        args = script.parse_args([])

        self.assertEqual(
            args.run_dir,
            pathlib.Path("~/Data/MPAS/poisson_charge_coupled_supercell/run"),
        )
        self.assertEqual(args.run_duration, "00_00:10:00")
        self.assertEqual(args.electrostatic_interval, 60.0)
        self.assertEqual(args.output_interval, "00:10:00")

    def test_registry_electrostatic_fields_are_time_dependent_for_dynamic_output(self):
        tree = ET.parse(REGISTRY_PATH)
        electrostatic = tree.find(".//var_struct[@name='electrostatic']")
        self.assertIsNotNone(electrostatic)

        for name in ("rho_charge", "phi", "E_normal", "E_vector"):
            var = electrostatic.find(f"./var[@name='{name}']")
            self.assertIsNotNone(var)
            self.assertIn("Time", var.attrib["dimensions"].split())

    def test_registry_has_poisson_altitude_coordinate_fields(self):
        tree = ET.parse(REGISTRY_PATH)
        electrostatic = tree.find(".//var_struct[@name='electrostatic']")
        self.assertIsNotNone(electrostatic)

        expected = {
            "poisson_zmid": "nVertLevels nCells Time",
            "poisson_air_thickness": "nVertLevels nCells Time",
        }
        for name, dims in expected.items():
            var = electrostatic.find(f"./var[@name='{name}']")
            self.assertIsNotNone(var, f"missing {name}")
            self.assertEqual(var.attrib["dimensions"], dims)


if __name__ == "__main__":
    unittest.main()
