#!/usr/bin/env python3
"""Unit tests for the 2-D terrain-following coordinate prototype (Task A1).

Run with:
    ~/miniconda3/envs/mpas/bin/python -m unittest \
        src.core_atmosphere.electrostatic.tests.test_terrain_operator -v
"""

import importlib.util
import math
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


class TerrainOperatorTests(unittest.TestCase):
    """Tests for Task A2: A = Gᵀ W G on the terrain grid.

    Checks (a) symmetric + SPD on a non-trivial hill and (b) that the operator
    reduces EXACTLY (< 1e-12) to the standard finite-volume 5-point Laplacian
    when the slope is zero (hill_height = 0).
    """

    GROUND_CELL = 0

    @staticmethod
    def _flat_5point_reference(nx, nz, dx, dzeta, eps, ground_cell):
        """Independent finite-volume 5-point Laplacian on the flat grid.

        Periodic in x (every cell has two x-faces), no-flux (Neumann) at the
        top and bottom in zeta (boundary cells have only one z-face), then
        grounded by decoupling one cell (zero its row/column, keep diagonal).

        Cell ordering matches the operator: c(i, k) = i*nz + k.
        Face weights are eps * physical control volume = eps * dx * dzeta, so
        the coefficients are cx = eps*dzeta/dx (x) and cz = eps*dx/dzeta (z).
        """
        n = nx * nz
        A = np.zeros((n, n))
        cx = eps * dzeta / dx
        cz = eps * dx / dzeta

        def c(i, k):
            return i * nz + k

        for i in range(nx):
            im = (i - 1) % nx
            ip = (i + 1) % nx
            for k in range(nz):
                cc = c(i, k)
                # x: periodic, always two faces
                A[cc, cc] += 2.0 * cx
                A[cc, c(ip, k)] -= cx
                A[cc, c(im, k)] -= cx
                # z: Neumann (no flux) at top/bottom
                if k > 0:
                    A[cc, cc] += cz
                    A[cc, c(i, k - 1)] -= cz
                if k < nz - 1:
                    A[cc, cc] += cz
                    A[cc, c(i, k + 1)] -= cz

        # Ground one cell: decouple it, keep its (positive) diagonal.
        diag_g = A[ground_cell, ground_cell]
        A[ground_cell, :] = 0.0
        A[:, ground_cell] = 0.0
        A[ground_cell, ground_cell] = diag_g
        return A

    def test_operator_symmetric_and_spd_on_hill(self):
        """A is symmetric and strictly positive-definite on a 0.3*H hill."""
        script = load_script()

        nx, nz = 24, 16
        hill_height, L, H = 0.3, 1.0, 1.0
        eps = 1.0
        g = script.terrain_grid(nx, nz, hill_height, L, H)

        A = script.assemble_terrain_operator(g, eps, ground_cell=self.GROUND_CELL)

        # Symmetry: no entry of |A - A.T| exceeds 1e-10.
        self.assertEqual((abs(A - A.T) > 1e-10).nnz, 0)

        # SPD: smallest eigenvalue strictly positive. Dense eigvalsh is exact
        # and deterministic for this small grid (n = nx*nz = 384).
        min_eig = np.linalg.eigvalsh(A.toarray()).min()
        self.assertGreater(min_eig, 0.0)

    def test_flat_grid_reduces_to_5point_laplacian(self):
        """hill_height=0 => A equals the standard FV 5-point Laplacian (< 1e-12).

        Uses dx != dzeta (nx=16, nz=12 on a unit square) so the x/z scaling is
        exercised independently.
        """
        script = load_script()

        nx, nz = 16, 12
        L, H = 1.0, 1.0
        eps = 1.0
        g = script.terrain_grid(nx, nz, hill_height=0.0, L=L, H=H)

        A = script.assemble_terrain_operator(g, eps, ground_cell=self.GROUND_CELL)

        ref = self._flat_5point_reference(
            nx, nz, g.dx, g.dzeta, eps, self.GROUND_CELL
        )

        max_err = np.max(np.abs(A.toarray() - ref))
        self.assertLess(max_err, 1e-12)


class TerrainMMSOrderTests(unittest.TestCase):
    """Task A3 — THE GATE: operator-only MMS convergence order.

    The flat case (hill_fraction=0) must reproduce ~2nd order; it is the
    honesty guard on the residual measure (V_cell scaling + interior mask).
    The hill case (0.3*H) is the gate proper: GO if the slope-corrected
    operator stays ~2nd order on the sloped mesh.
    """

    def test_flat_operator_is_second_order(self):
        """Honesty guard: flat grid reproduces ~2nd order (validates b/V_cell/mask)."""
        script = load_script()
        slope = script.terrain_mms_order(hill_fraction=0.0)
        self.assertGreater(slope, 1.9)
        self.assertLess(slope, 2.1)

    @unittest.expectedFailure  # NO-GO 2026-06-18; see notes/2026-06-18-poisson-terrain-prototype.md
    def test_terrain_operator_is_second_order_on_hill(self):
        """THE GATE (NO-GO): the naive column-averaged slope-corrected GtWG is SPD
        but NOT consistent on a hill (achieved slope ~-0.5 vs the >=1.9 target).
        Expected to PASS once approach A (a consistent mimetic cross-term grounded
        in MPAS zb/zb3) lands -- remove this decorator then."""
        script = load_script()
        slope = script.terrain_mms_order(hill_fraction=0.3)
        self.assertGreaterEqual(slope, 1.9)


class MFDTerrainTests(unittest.TestCase):
    """Bake-off: the MFD operator on the x-z terrain mesh."""

    @unittest.expectedFailure  # BLOCKED 2026-06-18; see notes/2026-06-18-poisson-mfd-bakeoff.md
    def test_mfd_flat_honesty_guard_is_second_order(self):
        """FLAT HONESTY GUARD (BLOCKED): the A0 mimetic block (calibrated to the
        HALF-cell offset C = fc-centroid, |C|=dx/2) composed with the plain-
        difference d0 (phi_hi-phi_lo over the FULL cell distance dx) does NOT
        reduce to the flat 5-point Laplacian: the shared-face coupling is over-
        counted ~4x and T0's off-diagonals add a spurious distance-2 stencil.
        Flat slope ~ -0.01 (errors stuck ~0.83). Threshold 1.9-2.1 is correct and
        left intact; passing it requires a redesign of the A0 d0/Hodge scaling
        contract (out of A1 scope -- do not guess). See the bake-off note."""
        script = load_script()
        slope = script.terrain_mms_order_mfd(hill_fraction=0.0, variant="mfd")
        self.assertGreater(slope, 1.9)
        self.assertLess(slope, 2.1)

    def test_mfd_operator_is_spd_on_hill(self):
        """Structural: each cell's mimetic block is SPD, so the grounded
        A = d0ᵀ H d0 is symmetric + SPD. This holds even though the operator is
        mis-scaled vs the Laplacian (the consistency defect is in the magnitude/
        stencil, not in symmetry or definiteness)."""
        script = load_script()
        nx, nz = 24, 16
        g = script.terrain_grid(nx, nz, hill_height=0.3, L=1.0, H=1.0)
        A = script.mfd_terrain_operator(g, eps=1.0, variant="mfd")
        self.assertEqual((abs(A - A.T) > 1e-9).nnz, 0)
        from numpy.linalg import eigvalsh
        self.assertGreater(eigvalsh(A.toarray()).min(), 0.0)

    @unittest.expectedFailure  # NO-GO/BLOCKED 2026-06-18; see notes/2026-06-18-poisson-mfd-bakeoff.md
    def test_mfd_terrain_is_second_order_on_hill(self):
        """THE TERRAIN GATE (NO-GO/BLOCKED): with the A0 block + plain-difference
        d0 as briefed, the MFD operator is not even consistent on the FLAT grid
        (see the flat guard), so this hill gate cannot validate the slope physics.
        Hill slope ~ -1.56 (errors DIVERGE 2.9->7.6->22.2 on refinement). Root
        cause is the half-cell-vs-full-difference scaling mismatch at the A0/A1
        interface, not the terrain geometry. Threshold >=1.9 left intact. Resolve
        the A0 scaling contract first, then re-run. See the bake-off note."""
        script = load_script()
        slope = script.terrain_mms_order_mfd(hill_fraction=0.3, variant="mfd")
        self.assertGreaterEqual(slope, 1.9)


class WhitneyTerrainTests(unittest.TestCase):
    """Bake-off pivot: the lowest-order Whitney/P1 scalar Hodge operator on the
    x-z terrain mesh, built on a triangulation of the physical cell centres.

    A = d0ᵀ ⋆₁ d0 with PLAIN ±1 incidence d0 and ⋆₁ the DIAGONAL scalar edge
    Hodge (l_e/d_e, the cotangent weight) assembled per triangle from the proven
    lowest-order Whitney/DEC primitive.  Unlike variant "mfd", the slope enters
    through the physical triangle geometry and the operator is consistent.
    """

    def test_whitney_flat_honesty_guard_is_second_order(self):
        """Flat (hill=0) Whitney operator must reproduce slope ~2.0 (validates the
        residual measure + triangulation + diagonal-Hodge assembly)."""
        script = load_script()
        slope = script.terrain_mms_order_mfd(hill_fraction=0.0, variant="whitney")
        self.assertGreater(slope, 1.9)
        self.assertLess(slope, 2.1)

    def test_whitney_operator_is_spd_on_hill(self):
        """The grounded A = d0ᵀ ⋆₁ d0 is symmetric and SPD on the 0.3*H hill."""
        script = load_script()
        nx, nz = 24, 16
        g = script.terrain_grid(nx, nz, hill_height=0.3, L=1.0, H=1.0)
        A = script.mfd_terrain_operator(g, eps=1.0, variant="whitney")
        self.assertEqual((abs(A - A.T) > 1e-9).nnz, 0)
        from numpy.linalg import eigvalsh
        self.assertGreater(eigvalsh(A.toarray()).min(), 0.0)

    def test_whitney_terrain_is_second_order_on_hill(self):
        """THE TERRAIN GATE (GO). The Whitney/P1 scalar Hodge operator stays
        ~2nd order on the 0.3*H hill (slope ~1.96), vs the naive GtWG's -0.52 and
        the mis-formulated variant-M's divergence."""
        script = load_script()
        slope = script.terrain_mms_order_mfd(hill_fraction=0.3, variant="whitney")
        self.assertGreaterEqual(slope, 1.9)


class TerrainSolutionErrorTests(unittest.TestCase):
    """Solution-error MMS gate for the terrain operator (Phase B / solution accuracy)."""

    def test_terrain_solution_error_flat_is_second_order(self):
        """Flat (hill=0) Whitney terrain solution-error slope must be ~2nd order."""
        script = load_script()
        slope = script.terrain_solution_error_order(hill_fraction=0.0, variant="whitney")
        self.assertGreater(slope, 1.9)
        self.assertLess(slope, 2.1)

    def test_terrain_solution_error_hill_slope_is_positive(self):
        """Hill (0.3*H) Whitney terrain solution-error slope is finite and positive (diagnostic)."""
        script = load_script()
        slope = script.terrain_solution_error_order(hill_fraction=0.3, variant="whitney")
        self.assertTrue(math.isfinite(slope), f"hill solution-error slope is NaN/inf: {slope}")
        self.assertGreater(slope, 0.0, f"hill solution-error slope non-positive: {slope}")


if __name__ == "__main__":
    unittest.main()
