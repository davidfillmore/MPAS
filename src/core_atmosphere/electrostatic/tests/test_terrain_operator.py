#!/usr/bin/env python3
"""Unit tests for the 2-D terrain-following coordinate prototype (Task A1).

Run with:
    ~/miniconda3/envs/mpas/bin/python -m unittest \
        src.core_atmosphere.electrostatic.tests.test_terrain_operator -v
"""

import importlib.util
import pathlib
import unittest

import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/prototype_terrain_operator.py"
)


def load_script():
    spec = importlib.util.spec_from_file_location(
        "prototype_terrain_operator", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TerrainGridTests(unittest.TestCase):
    """Tests for terrain_grid(nx, nz, hill_height, L, H)."""

    def test_dzdx_matches_fd_of_interface_heights_in_x(self):
        """Analytic dzdx at interfaces agrees with centred periodic FD of z_int in x
        to better than 1e-6.

        Uses nx=256 cells with hill_height=0.001 so the 2nd-order FD truncation
        error O(dx^2 * hill_height * (2pi/L)^3 / 6) ≈ 6.3e-7 is safely below 1e-6.
        """
        script = load_script()

        nx, nz = 256, 16
        hill_height, L, H = 0.001, 1.0, 1.0
        g = script.terrain_grid(nx, nz, hill_height, L, H)

        # Centred periodic FD of interface heights in x
        dzdx_fd = (
            np.roll(g.z_int, -1, axis=0) - np.roll(g.z_int, 1, axis=0)
        ) / (2.0 * g.dx)

        max_err = np.max(np.abs(dzdx_fd - g.dzdx_int))
        self.assertLess(max_err, 1e-6)

    def test_dzdx_at_top_interface_is_exactly_zero(self):
        """dzdx vanishes at the top interface (zeta = H) because (1 - H/H) = 0."""
        script = load_script()

        nx, nz = 64, 32
        hill_height, L, H = 0.1, 1.0, 1.0
        g = script.terrain_grid(nx, nz, hill_height, L, H)

        # At the top interface k = nz: zeta_int = H, factor (1 - zeta/H) = 0 exactly
        np.testing.assert_array_equal(g.dzdx_int[:, -1], 0.0)

    def test_flat_grid_has_zero_dzdx_and_unit_zz(self):
        """A hill_height=0 grid has dzdx=0 and zz=1 everywhere (cell centres)."""
        script = load_script()

        nx, nz = 32, 16
        g = script.terrain_grid(nx, nz, hill_height=0.0, L=1.0, H=1.0)

        np.testing.assert_array_equal(g.dzdx, np.zeros((nx, nz)))
        np.testing.assert_array_equal(g.zz, np.ones((nx, nz)))


if __name__ == "__main__":
    unittest.main()
