"""Task 7: 3-D shell mesh tests (mesh part).

Tests cover both icosahedral and SCVT-dual shell mesh families.
Each surface triangulation is extruded radially into conforming tetrahedra
via the Freudenthal–Kuhn sorted-vertex prism split.

Checks:
  1. topology.dim == 3, geometry.dim == 3
  2. Every tet has positive signed volume
  3. Total tet volume is within the polyhedral band of the analytic shell volume
  4. Both inner (r≈R) and outer (r≈R+H) boundary facets exist
"""
import math
import pathlib

import numpy as np
import pytest

from fecore.mesh import icosa, scvt_dual

ROOT = pathlib.Path("~/Data/MPAS/poisson_tier_A2_scvt/meshes").expanduser()
GRID_480 = ROOT / "480km" / "grid.nc"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _shell_volume(R: float, H: float) -> float:
    """Analytic shell volume (4/3)π((R+H)³ − R³)."""
    return (4.0 / 3.0) * math.pi * ((R + H) ** 3 - R ** 3)


def _tet_signed_volumes(coords: np.ndarray, cells: np.ndarray) -> np.ndarray:
    """Return signed volume of each tet.  Positive means correct orientation.

    V = det([v1-v0, v2-v0, v3-v0]) / 6
    """
    v0 = coords[cells[:, 0]]
    v1 = coords[cells[:, 1]]
    v2 = coords[cells[:, 2]]
    v3 = coords[cells[:, 3]]
    a = v1 - v0
    b = v2 - v0
    c = v3 - v0
    return (a[:, 0] * (b[:, 1] * c[:, 2] - b[:, 2] * c[:, 1])
            - a[:, 1] * (b[:, 0] * c[:, 2] - b[:, 2] * c[:, 0])
            + a[:, 2] * (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0])) / 6.0


def _get_coords_and_cells(msh):
    """Extract coordinates and cell→vertex connectivity from a tet mesh."""
    msh.topology.create_connectivity(3, 0)
    c2v = msh.topology.connectivity(3, 0)
    dofmap = msh.geometry.dofmap      # geometry node indices per cell
    coords = msh.geometry.x           # (n_geo_nodes, 3)
    n_cells = msh.topology.index_map(3).size_local
    cells = np.array([dofmap[i] for i in range(n_cells)], dtype=np.int64)
    return coords, cells


# ---------------------------------------------------------------------------
# icosa shell tests
# ---------------------------------------------------------------------------

class TestIcosaShell:
    """Tests for icosa_shell_mesh (unit sphere, R=1, H=0.25, n_layers=2)."""

    @pytest.fixture(scope="class")
    @classmethod
    def msh(cls):
        return icosa.icosa_shell_mesh(refine=2, n_layers=2, H=0.25)

    def test_topology_dim_3(self, msh):
        assert msh.topology.dim == 3

    def test_gdim_3(self, msh):
        assert msh.geometry.dim == 3

    def test_all_tets_positive_volume(self, msh):
        coords, cells = _get_coords_and_cells(msh)
        vols = _tet_signed_volumes(coords, cells)
        n_neg = int(np.sum(vols <= 0))
        assert n_neg == 0, (
            f"{n_neg} tets have non-positive signed volume (max neg = "
            f"{vols[vols <= 0].min() if n_neg else 'n/a'})"
        )

    def test_total_volume_in_shell_band(self, msh):
        R, H = 1.0, 0.25
        V_analytic = _shell_volume(R, H)
        coords, cells = _get_coords_and_cells(msh)
        vols = _tet_signed_volumes(coords, cells)
        V_tet = float(np.sum(np.abs(vols)))
        ratio = V_tet / V_analytic
        assert 0.80 < ratio < 1.01, (
            f"Volume ratio {ratio:.4f} outside (0.80, 1.01). "
            f"V_tet={V_tet:.6f}, V_analytic={V_analytic:.6f}"
        )

    def test_inner_and_outer_boundary_facets_exist(self, msh):
        R, H = 1.0, 0.25
        r_tol = 0.02 * R  # 2 % tolerance

        coords = msh.geometry.x
        radii = np.linalg.norm(coords, axis=1)

        has_inner = bool(np.any(np.abs(radii - R) < r_tol))
        has_outer = bool(np.any(np.abs(radii - (R + H)) < r_tol))

        assert has_inner, (
            f"No nodes found near inner boundary r={R}. "
            f"Min radius = {radii.min():.4f}"
        )
        assert has_outer, (
            f"No nodes found near outer boundary r={R+H}. "
            f"Max radius = {radii.max():.4f}"
        )


# ---------------------------------------------------------------------------
# SCVT shell tests
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not GRID_480.exists(), reason="SCVT 480km mesh not present")
class TestScvtShell:
    """Tests for scvt_shell_mesh (Earth-radius SCVT, n_layers=2)."""

    @pytest.fixture(scope="class")
    @classmethod
    def msh_and_params(cls):
        H = 2.0e4  # 20 km shell
        msh = scvt_dual.scvt_shell_mesh(GRID_480, n_layers=2, H=H)
        R = scvt_dual.sphere_radius(GRID_480)
        return msh, R, H

    def test_topology_dim_3(self, msh_and_params):
        msh, R, H = msh_and_params
        assert msh.topology.dim == 3

    def test_gdim_3(self, msh_and_params):
        msh, R, H = msh_and_params
        assert msh.geometry.dim == 3

    def test_all_tets_positive_volume(self, msh_and_params):
        msh, R, H = msh_and_params
        coords, cells = _get_coords_and_cells(msh)
        vols = _tet_signed_volumes(coords, cells)
        n_neg = int(np.sum(vols <= 0))
        assert n_neg == 0, (
            f"{n_neg} tets have non-positive signed volume (max neg = "
            f"{vols[vols <= 0].min() if n_neg else 'n/a'})"
        )

    def test_total_volume_in_shell_band(self, msh_and_params):
        msh, R, H = msh_and_params
        V_analytic = _shell_volume(R, H)
        coords, cells = _get_coords_and_cells(msh)
        vols = _tet_signed_volumes(coords, cells)
        V_tet = float(np.sum(np.abs(vols)))
        ratio = V_tet / V_analytic
        assert 0.80 < ratio < 1.01, (
            f"Volume ratio {ratio:.4f} outside (0.80, 1.01). "
            f"V_tet={V_tet:.6e}, V_analytic={V_analytic:.6e}"
        )

    def test_inner_and_outer_boundary_facets_exist(self, msh_and_params):
        msh, R, H = msh_and_params
        r_tol = 0.005 * R  # 0.5 % tolerance on Earth radius

        coords = msh.geometry.x
        radii = np.linalg.norm(coords, axis=1)

        has_inner = bool(np.any(np.abs(radii - R) < r_tol))
        has_outer = bool(np.any(np.abs(radii - (R + H)) < r_tol))

        assert has_inner, (
            f"No nodes found near inner boundary r={R:.2e}. "
            f"Min radius = {radii.min():.6e}"
        )
        assert has_outer, (
            f"No nodes found near outer boundary r={R+H:.2e}. "
            f"Max radius = {radii.max():.6e}"
        )
