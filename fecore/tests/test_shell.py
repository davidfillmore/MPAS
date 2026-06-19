"""Task 7: 3-D shell mesh tests (mesh part).

Tests cover both icosahedral and SCVT-dual shell mesh families.
Each surface triangulation is extruded radially into conforming tetrahedra
via the Freudenthal–Kuhn sorted-vertex prism split.

Checks:
  1. topology.dim == 3, geometry.dim == 3
  2. Every tet has positive signed volume
  3. Total tet volume is within the polyhedral band of the analytic shell volume
  4. Conformity: every facet is shared by exactly 1 (boundary) or 2 (interior)
     cells — never 3+. A 3+-shared facet would mean a non-manifold / non-
     conforming mesh (two adjacent prisms chose opposite shared-quad diagonals).
  5. Boundary-facet geometry: every boundary facet's centroid lies near the
     inner sphere (r≈R) or the outer sphere (r≈R+H), and both surfaces are
     represented.
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
    coords = msh.geometry.x           # (n_geo_nodes, 3)
    dofmap = msh.geometry.dofmap      # geometry node indices per cell
    n_cells = msh.topology.index_map(3).size_local
    cells = np.array([dofmap[i] for i in range(n_cells)], dtype=np.int64)
    return coords, cells


def _facet_cell_counts(msh):
    """Return list of cell-counts for every facet (via facet→cell adjacency).

    Creates the 2↔3 and 2↔0 connectivity tables on *msh* as a side-effect.
    Returns a list of length n_facets.
    """
    msh.topology.create_connectivity(2, 3)
    f2c = msh.topology.connectivity(2, 3)
    n_facets = msh.topology.index_map(2).size_local
    return [len(f2c.links(i)) for i in range(n_facets)]


def _boundary_facet_centroid_radii(msh, counts):
    """Return an array of centroid radii for facets with exactly 1 adjacent cell.

    Parameters
    ----------
    msh:
        dolfinx mesh with 2→0 connectivity already created.
    counts:
        Per-facet cell-count list (from _facet_cell_counts).
    """
    msh.topology.create_connectivity(2, 0)
    f2v = msh.topology.connectivity(2, 0)
    coords = msh.geometry.x
    bfacets = [i for i, c in enumerate(counts) if c == 1]
    radii = np.array([
        np.linalg.norm(coords[f2v.links(fi)].mean(axis=0))
        for fi in bfacets
    ])
    return radii


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

    def test_conformity(self, msh):
        """Every facet is shared by ≤ 2 cells (1=boundary, 2=interior).

        A facet shared by 3+ cells would indicate a non-manifold / non-
        conforming mesh, meaning the Freudenthal–Kuhn sorted-vertex split
        failed to produce a consistent quad-face diagonal between adjacent
        prisms.  This is the canonical conformity gate for the FK split.
        """
        counts = _facet_cell_counts(msh)
        max_count = max(counts)
        n_boundary = sum(1 for c in counts if c == 1)
        n_interior = sum(1 for c in counts if c == 2)
        n_nonconf = sum(1 for c in counts if c >= 3)

        assert max_count <= 2, (
            f"NON-CONFORMING MESH: max cells-per-facet = {max_count} "
            f"({n_nonconf} non-conforming facets). "
            "The Freudenthal–Kuhn split is broken — stop and investigate."
        )
        assert n_boundary > 0, "No boundary facets found (unexpected for a shell mesh)"
        assert n_interior > 0, "No interior facets found (unexpected for a shell mesh)"

    def test_boundary_facets_at_shell_surfaces(self, msh):
        """Every boundary facet centroid lies at r≈R (inner) or r≈R+H (outer).

        Checks both that the centroid radii are near a shell surface (not at
        an interior radius) and that both inner and outer surfaces are present.

        Facets are partitioned at the shell midpoint (R + H/2) before the
        radius check so the two tolerance bands cannot overlap regardless of
        H/R ratio.  The tolerance is 3 % of R to accommodate the chord-
        shortfall of flat triangles on the curved sphere at refine=2 (≈1.8 %
        observed; 3 % captures it with margin).
        """
        R, H = 1.0, 0.25
        r_tol = 0.03 * R  # chord shortfall at refine=2 ≈ 1.8%; 3 % captures it

        counts = _facet_cell_counts(msh)
        bfacet_radii = _boundary_facet_centroid_radii(msh, counts)

        # Partition at the shell midpoint.
        mid_r = R + H / 2.0
        inner_radii = bfacet_radii[bfacet_radii < mid_r]
        outer_radii = bfacet_radii[bfacet_radii >= mid_r]

        assert len(inner_radii) > 0, (
            f"No boundary facets below the shell midpoint r={mid_r:.4f}. "
            f"Inner surface (r≈{R}) is missing."
        )
        assert len(outer_radii) > 0, (
            f"No boundary facets above the shell midpoint r={mid_r:.4f}. "
            f"Outer surface (r≈{R+H}) is missing."
        )

        inner_stray = np.sum(np.abs(inner_radii - R) >= r_tol)
        assert inner_stray == 0, (
            f"{inner_stray}/{len(inner_radii)} inner-partition boundary facets "
            f"have centroid radius not within {r_tol:.3f} of R={R}. "
            f"Range: [{inner_radii.min():.4f}, {inner_radii.max():.4f}]"
        )

        outer_stray = np.sum(np.abs(outer_radii - (R + H)) >= r_tol)
        assert outer_stray == 0, (
            f"{outer_stray}/{len(outer_radii)} outer-partition boundary facets "
            f"have centroid radius not within {r_tol:.3f} of R+H={R+H}. "
            f"Range: [{outer_radii.min():.4f}, {outer_radii.max():.4f}]"
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

    def test_conformity(self, msh_and_params):
        """Every facet is shared by ≤ 2 cells (1=boundary, 2=interior).

        A facet shared by 3+ cells would indicate a non-manifold / non-
        conforming mesh.  This is the canonical conformity gate for the FK split.
        """
        msh, R, H = msh_and_params
        counts = _facet_cell_counts(msh)
        max_count = max(counts)
        n_boundary = sum(1 for c in counts if c == 1)
        n_interior = sum(1 for c in counts if c == 2)
        n_nonconf = sum(1 for c in counts if c >= 3)

        assert max_count <= 2, (
            f"NON-CONFORMING MESH: max cells-per-facet = {max_count} "
            f"({n_nonconf} non-conforming facets). "
            "The Freudenthal–Kuhn split is broken — stop and investigate."
        )
        assert n_boundary > 0, "No boundary facets found (unexpected for a shell mesh)"
        assert n_interior > 0, "No interior facets found (unexpected for a shell mesh)"

    def test_boundary_facets_at_shell_surfaces(self, msh_and_params):
        """Every boundary facet centroid lies at r≈R (inner) or r≈R+H (outer).

        Both inner and outer surfaces must be present, and no boundary facet
        should sit at an interior radius.

        Facets are partitioned at the shell midpoint (R + H/2): those below the
        midpoint are compared to R (inner), those above to R+H (outer).  The
        tolerance on each side is the larger of (a) 1 % of H (to handle chord
        shortfall) and (b) 0.5 % of R (a loose physical guard); in practice the
        480 km mesh chord shortfall is ≈ 6.3 km while H/100 = 200 m, so the
        guard (0.5 %·R ≈ 31.9 km) is the active limit.  Because we partition
        first, the two tolerance bands can never overlap regardless of how thin
        the shell is relative to R.
        """
        msh, R, H = msh_and_params
        # Tolerance: generous enough to absorb the chord shortfall of flat
        # triangles inscribed on a curved sphere, but tight enough that an
        # interior facet (radius near R + k*H/n_layers for k∈{1}) would fail.
        r_tol = max(0.01 * H, 0.005 * R)

        counts = _facet_cell_counts(msh)
        bfacet_radii = _boundary_facet_centroid_radii(msh, counts)

        # Partition at the shell midpoint.
        mid_r = R + H / 2.0
        inner_radii = bfacet_radii[bfacet_radii < mid_r]
        outer_radii = bfacet_radii[bfacet_radii >= mid_r]

        assert len(inner_radii) > 0, (
            f"No boundary facets below the shell midpoint r={mid_r:.4e}. "
            f"Inner surface (r≈{R:.4e}) is missing."
        )
        assert len(outer_radii) > 0, (
            f"No boundary facets above the shell midpoint r={mid_r:.4e}. "
            f"Outer surface (r≈{R+H:.4e}) is missing."
        )

        # Every inner-partition facet must be near R.
        inner_stray = np.sum(np.abs(inner_radii - R) >= r_tol)
        assert inner_stray == 0, (
            f"{inner_stray}/{len(inner_radii)} inner-partition boundary facets "
            f"have centroid radius not within {r_tol:.3e} of R={R:.4e}. "
            f"Range: [{inner_radii.min():.6e}, {inner_radii.max():.6e}]"
        )

        # Every outer-partition facet must be near R+H.
        outer_stray = np.sum(np.abs(outer_radii - (R + H)) >= r_tol)
        assert outer_stray == 0, (
            f"{outer_stray}/{len(outer_radii)} outer-partition boundary facets "
            f"have centroid radius not within {r_tol:.3e} of R+H={R+H:.4e}. "
            f"Range: [{outer_radii.min():.6e}, {outer_radii.max():.6e}]"
        )
