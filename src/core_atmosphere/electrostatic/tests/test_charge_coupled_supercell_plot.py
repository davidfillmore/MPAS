#!/usr/bin/env python3
"""Unit tests for charge-coupled supercell figure helpers."""

import importlib.util
import pathlib
import tempfile
import unittest
from unittest import mock

import netCDF4 as nc
import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/plot_charge_coupled_supercell.py"
)


def load_script():
    spec = importlib.util.spec_from_file_location("plot_charge_coupled_supercell", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ChargeCoupledSupercellPlotTests(unittest.TestCase):
    def test_read_figure_fields_uses_final_time_and_liquid_water(self):
        script = load_script()

        with tempfile.TemporaryDirectory() as tmp:
            output = pathlib.Path(tmp) / "output.nc"
            write_synthetic_output(output)
            fields = script.read_figure_fields(output)

        self.assertEqual(fields["liquid_water_content"].shape, (3, 4))
        self.assertEqual(fields["liquid_water_mixing_ratio"].shape, (3, 4))
        self.assertEqual(fields["rho_charge"].shape, (3, 4))
        self.assertEqual(fields["phi"].shape, (3, 4))
        self.assertTrue(np.allclose(fields["air_density"], 1.2))
        self.assertTrue(np.allclose(fields["liquid_water_mixing_ratio"], 5.0))
        self.assertTrue(np.allclose(fields["liquid_water_content"], 6.0))
        self.assertTrue(np.allclose(fields["rho_charge"], 4.0))
        self.assertTrue(np.allclose(fields["phi"], -7.0))
        self.assertEqual(fields["time_index"], -1)
        self.assertEqual(fields["simulation_time_label"], "Simulation Time: 10 min")

    def test_select_cross_section_passes_through_peak_liquid_water_column(self):
        script = load_script()
        fields = {
            "xCell": np.array([0.0, 1.0, 0.0, 1.0]),
            "yCell": np.array([0.0, 0.0, 1.0, 1.0]),
            "liquid_water_content": np.array(
                [
                    [0.0, 4.0, 0.0, 1.0],
                    [0.0, 5.0, 0.0, 1.0],
                ]
            ),
            "rho_charge": np.zeros((2, 4)),
        }

        indices = script.select_cross_section(fields)

        self.assertEqual(indices.tolist(), [1, 3])

    def test_section_coordinates_are_level_cell_samples(self):
        script = load_script()

        fields = {
            "yCell": np.array([0.0, 1.0]),
            "zMid": np.array([[0.5, 0.6], [1.5, 1.6], [2.5, 2.6]]),
        }

        y, z = script.section_coordinates(fields, np.array([0, 1]))

        self.assertEqual(y.shape, (3, 2))
        self.assertEqual(z.shape, (3, 2))
        self.assertTrue(np.allclose(y[:, 0], 0.0))
        self.assertTrue(np.allclose(y[:, 1], 1.0))
        self.assertTrue(np.allclose(z, fields["zMid"]))

    def test_signed_line_levels_scale_each_charge_sign_independently(self):
        script = load_script()

        negative, positive = script.signed_line_levels(
            np.array([-1.647, -0.7, -0.1, 0.0, 0.0042])
        )

        self.assertTrue(np.allclose(negative, [-1.0, -0.5, -0.2, -0.1, -0.05, -0.02, -0.01]))
        self.assertTrue(np.allclose(positive, [0.001, 0.002]))

    def test_charge_contour_style_uses_red_negative_and_gray_positive_ramps(self):
        script = load_script()

        self.assertEqual(script.NEGATIVE_CHARGE_CONTOUR_COLORS[0], "#fcae91")
        self.assertEqual(script.NEGATIVE_CHARGE_CONTOUR_COLORS[-1], "#67000d")
        self.assertEqual(script.POSITIVE_CHARGE_CONTOUR_COLORS[0], "#bdbdbd")
        self.assertEqual(script.POSITIVE_CHARGE_CONTOUR_COLORS[-1], "#000000")
        self.assertEqual(script.CHARGE_CONTOUR_LINESTYLE, "solid")

    def test_lwc_colormap_starts_at_white(self):
        script = load_script()

        cmap = script.lwc_colormap()

        self.assertTrue(np.allclose(cmap(0.0)[:3], (1.0, 1.0, 1.0)))

    def test_positive_filled_levels_use_round_upper_bound_and_ticks(self):
        script = load_script()

        levels, ticks = script.positive_filled_levels(np.array([0.0, 2.557]))

        self.assertAlmostEqual(float(levels[0]), 0.0)
        self.assertAlmostEqual(float(levels[-1]), 3.0)
        self.assertTrue(np.allclose(ticks, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]))

    def test_signed_filled_levels_are_centered_on_zero(self):
        script = load_script()

        levels, ticks, norm = script.signed_filled_levels(np.array([-351.3, -0.02]))

        self.assertAlmostEqual(float(levels[0]), -400.0)
        self.assertAlmostEqual(float(levels[-1]), 0.0)
        self.assertTrue(np.allclose(ticks, [-400.0, -300.0, -200.0, -100.0, 0.0]))
        self.assertEqual(norm.vcenter, 0.0)

    def test_plot_labels_use_formatted_superscripts_and_subscripts(self):
        script = load_script()

        labels = (
            script.LWC_COLORBAR_LABEL,
            script.RHO_CONTOUR_LEGEND_TITLE,
            script.PHI_COLORBAR_LABEL,
        )

        for label in labels:
            self.assertNotIn("^", label)
            self.assertNotIn("^-", label)

    def test_plot_titles_use_title_case_except_short_conjunctions(self):
        script = load_script()

        self.assertEqual(script.LEFT_PANEL_TITLE, "Liquid Water Content and Charge Source")
        self.assertEqual(script.RIGHT_PANEL_TITLE, "Potential and Electric Field")

    def test_plot_charge_coupled_output_omits_panel_letter_annotations(self):
        script = load_script()

        from matplotlib.axes import Axes

        labels = []
        original_text = Axes.text

        def record_text(ax, x, y, label, *args, **kwargs):
            labels.append(str(label))
            return original_text(ax, x, y, label, *args, **kwargs)

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output = tmp_path / "output.nc"
            plot = tmp_path / "charge_coupled.png"
            write_synthetic_output(output)

            with mock.patch.object(Axes, "text", record_text):
                script.plot_charge_coupled_output(output, plot)

        self.assertNotIn("(a)", labels)
        self.assertNotIn("(b)", labels)

    def test_xtime_to_simulation_time_label_formats_hours_and_minutes(self):
        script = load_script()

        label = script.xtime_to_simulation_time_label("0000-01-01_01:30:00")

        self.assertEqual(label, "Simulation Time: 1 h 30 min")

    def test_plot_charge_coupled_output_writes_png(self):
        script = load_script()

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output = tmp_path / "output.nc"
            plot = tmp_path / "charge_coupled.png"
            write_synthetic_output(output)

            result = script.plot_charge_coupled_output(output, plot)

            self.assertEqual(result, plot)
            self.assertTrue(plot.exists())
            self.assertGreater(plot.stat().st_size, 0)


def write_synthetic_output(path):
    x = np.array([0.0, 1.0, 0.0, 1.0])
    y = np.array([0.0, 0.0, 1.0, 1.0])
    zgrid = np.array(
        [
            [0.0, 1.0, 2.0, 3.0],
            [0.0, 1.2, 2.2, 3.2],
            [0.0, 1.0, 2.0, 3.0],
            [0.0, 1.2, 2.2, 3.2],
        ]
    )
    shape = (2, 4, 3)
    qc = np.zeros(shape)
    qr = np.zeros(shape)
    rho = np.zeros(shape)
    phi = np.zeros(shape)
    qc[1, :, :] = 2.0
    qr[1, :, :] = 3.0
    density = np.ones(shape)
    density[1, :, :] = 1.2
    rho[1, :, :] = 4.0
    phi[1, :, :] = -7.0
    e_vector = np.zeros((2, 4, 3, 3))
    e_vector[1, :, :, 0] = 1.0
    e_vector[1, :, :, 1] = 2.0
    e_vector[1, :, :, 2] = -1.0

    with nc.Dataset(path, "w") as ds:
        ds.createDimension("Time", 2)
        ds.createDimension("StrLen", 64)
        ds.createDimension("nCells", 4)
        ds.createDimension("nVertLevels", 3)
        ds.createDimension("nVertLevelsP1", 4)
        ds.createDimension("R3", 3)
        xtime = np.full((2, 64), b"\x00", dtype="S1")
        for row, text in enumerate(("0000-01-01_00:00:00", "0000-01-01_00:10:00")):
            xtime[row, : len(text)] = np.frombuffer(text.encode("ascii"), dtype="S1")
        ds.createVariable("xtime", "S1", ("Time", "StrLen"))[:] = xtime
        ds.createVariable("xCell", "f8", ("nCells",))[:] = x
        ds.createVariable("yCell", "f8", ("nCells",))[:] = y
        ds.createVariable("zgrid", "f8", ("nCells", "nVertLevelsP1"))[:] = zgrid
        ds.createVariable("qc", "f8", ("Time", "nCells", "nVertLevels"))[:] = qc
        ds.createVariable("qr", "f8", ("Time", "nCells", "nVertLevels"))[:] = qr
        ds.createVariable("rho", "f8", ("Time", "nCells", "nVertLevels"))[:] = density
        ds.createVariable("rho_charge", "f8", ("Time", "nCells", "nVertLevels"))[:] = rho
        ds.createVariable("phi", "f8", ("Time", "nCells", "nVertLevels"))[:] = phi
        ds.createVariable("E_vector", "f8", ("Time", "nCells", "nVertLevels", "R3"))[
            :
        ] = e_vector


if __name__ == "__main__":
    unittest.main()
