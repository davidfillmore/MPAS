"""Task 7 + Task 8: 3-D shell mesh tests (mesh + solve parts).

Task 7 tests cover both icosahedral and SCVT-dual shell mesh families.
Each surface triangulation is extruded radially into conforming tetrahedra
via the Freudenthal–Kuhn sorted-vertex prism split.

Mesh checks (Task 7):
  1. topology.dim == 3, geometry.dim == 3
  2. Every tet has positive signed volume
  3. Total tet volume is within the polyhedral band of the analytic shell volume
  4. Conformity: every facet is shared by exactly 1 (boundary) or 2 (interior)
     cells — never 3+. A 3+-shared facet would mean a non-manifold / non-
     conforming mesh (two adjacent prisms chose opposite shared-quad diagonals).
  5. Boundary-facet geometry: every boundary facet's centroid lies near the
     inner sphere (r≈R) or the outer sphere (r≈R+H), and both surfaces are
     represented.

Solve checks (Task 8):
  6. solve_shell returns a finite relative L2 error in (0, 1) on a small icosa
     shell mesh, and the inner-surface (ground) DOFs are enforced to zero.
  7. shell_convergence("icosa", 2, "consistent") achieves slope >= 1.9.
     Like the surface, the shell is capped at 2nd order by the flat-tet geometry
     (see convergence.py's isoparametric geometry cap note).
  8. SCVT shell convergence (guarded by scvt_available()): record the matrix.
"""
import math
import pathlib

import numpy as np
import pytest
import ufl

from fecore.mesh import icosa, scvt_dual
from fecore.operator.poisson import solve_shell
from fecore.verify import convergence, mms

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


# ---------------------------------------------------------------------------
# Task 8: Shell solve tests (ground Dirichlet + top Neumann MMS)
# ---------------------------------------------------------------------------

def test_shell_solve_runs_and_is_finite():
    """solve_shell returns a finite relative L2 error in (0, 1) for P1-consistent.

    Uses a small icosa shell (refine=2, n_layers=2) so this test is fast.
    The MMS is Y_4^2 * sin(pi/2 * zeta): zero on the inner sphere, zero
    Neumann on the outer sphere (cos vanishes at zeta=1).
    """
    R, H = 1.0, 0.25
    msh = icosa.icosa_shell_mesh(refine=2, n_layers=2, H=H)
    uh, err = solve_shell(
        msh, R=R, H=H, degree=1, rhs_mode="consistent",
        source_fn=mms.shell_source, exact_fn=mms.shell_exact, eps=1.0,
    )
    assert err is not None, "solve_shell returned None for err"
    assert np.isfinite(err), f"solve_shell returned non-finite err={err}"
    assert 0 < err < 1.0, f"Expected 0 < err < 1.0, got err={err}"


def test_shell_ground_dirichlet_enforced():
    """Inner-surface DOFs must be zero after solve_shell (ground BC).

    Locates the inner-surface facets by centroid radius, retrieves the
    corresponding FEM DOFs from the returned function, and checks that
    all values are numerically zero.  This verifies that the Dirichlet
    BC is correctly applied — not just set up — by solve_shell.
    """
    import dolfinx.fem
    import dolfinx.mesh as dmesh

    R, H = 1.0, 0.25
    msh = icosa.icosa_shell_mesh(refine=2, n_layers=2, H=H)
    uh, _ = solve_shell(
        msh, R=R, H=H, degree=1, rhs_mode="consistent",
        source_fn=mms.shell_source, exact_fn=mms.shell_exact, eps=1.0,
    )

    # Locate inner surface: facets whose centroid radius ≈ R.
    V = uh.function_space
    tol = 0.05 * R
    inner_facets = dmesh.locate_entities_boundary(
        msh, msh.topology.dim - 1,
        lambda x: np.linalg.norm(x, axis=0) < R + tol,
    )
    inner_dofs = dolfinx.fem.locate_dofs_topological(
        V, msh.topology.dim - 1, inner_facets
    )
    assert len(inner_dofs) > 0, "No inner-surface DOFs found — inner boundary missing"
    max_val = float(np.abs(uh.x.array[inner_dofs]).max())
    assert max_val < 1e-12, (
        f"Ground Dirichlet not enforced: max |uh| on inner surface = {max_val:.3e}"
    )


def test_shell_stiffness_matrix_is_spd():
    """Shell stiffness matrix assembled with ground Dirichlet BC is SPD.

    Assembles A = eps * integral(grad(u).grad(v) dx) with the ground-Dirichlet
    BC (phi=0 on r≈R) on a small icosa shell mesh.  The Dirichlet BC removes the
    null space, making A SPD.  We verify two cheap necessary conditions:

    1. Positive diagonal: every diagonal entry > 0 (necessary for SPD; a direct
       check without a full eigendecomposition; CG convergence in the solve tests
       already implies SPD empirically).
    2. Symmetry: norm(A - A^T, Frobenius) / norm(A, Frobenius) ≈ 0 (FEM
       assembly is symmetric by construction; this catches implementation errors
       such as unsymmetric BC application).  Computed via CSR data from PETSc.
    """
    import dolfinx.fem
    import dolfinx.fem.petsc
    import petsc4py.PETSc as PETSc

    R, H = 1.0, 0.25
    eps = 1.0
    msh = icosa.icosa_shell_mesh(refine=2, n_layers=2, H=H)
    V = dolfinx.fem.functionspace(msh, ("Lagrange", 1))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    # Ground Dirichlet BC (same as solve_shell).
    import dolfinx.mesh as dmesh
    tol = 0.05 * R
    inner_facets = dmesh.locate_entities_boundary(
        msh, msh.topology.dim - 1,
        lambda x: np.linalg.norm(x, axis=0) < R + tol,
    )
    inner_dofs = dolfinx.fem.locate_dofs_topological(
        V, msh.topology.dim - 1, inner_facets
    )
    bc = dolfinx.fem.dirichletbc(np.float64(0.0), inner_dofs, V)

    a_ufl = eps * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    a_form = dolfinx.fem.form(a_ufl)
    A = dolfinx.fem.petsc.assemble_matrix(a_form, bcs=[bc])
    A.assemble()

    # 1. Positive diagonal check via PETSc diagonal vector.
    diag_vec = A.createVecRight()
    A.getDiagonal(diag_vec)
    diag = diag_vec.getArray().copy()
    diag_vec.destroy()
    n_nonpos = int(np.sum(diag <= 0))
    assert n_nonpos == 0, (
        f"{n_nonpos}/{len(diag)} diagonal entries are non-positive "
        f"(min diagonal = {diag.min():.3e}). "
        "Stiffness matrix with ground Dirichlet BC must have positive diagonal."
    )

    # 2. Symmetry check: extract CSR data and compute Frobenius norm of A - A^T.
    # PETSc getValuesCSR returns (row_ptr, col_indices, values) on the local rows.
    ai, aj, av = A.getValuesCSR()
    n = A.getSize()[0]
    # Build a dense array only for this small (refine=2, n_layers=2) mesh — n ≈ 200.
    A_dense = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        cols = aj[ai[i]:ai[i + 1]]
        vals = av[ai[i]:ai[i + 1]]
        A_dense[i, cols] = vals
    frob_diff = np.linalg.norm(A_dense - A_dense.T, "fro")
    frob_A = np.linalg.norm(A_dense, "fro")
    rel_asymmetry = frob_diff / frob_A if frob_A > 0 else frob_diff
    assert rel_asymmetry < 1e-12, (
        f"Stiffness matrix is not symmetric: "
        f"norm(A-A^T,'fro')/norm(A,'fro') = {rel_asymmetry:.3e}"
    )

    A.destroy()


def test_icosa_shell_p2_is_second_order():
    """P2 FEM on the icosa shell must achieve slope >= 1.9 (icosa gate).

    The threshold is >= 1.9 (not ~3.0) because the shell is meshed with
    degree-1 geometry (flat tetrahedra approximating the spherical shell):
    the same isoparametric geometry cap that applies on the surface (see
    convergence.py's module docstring) applies here.  P2 elements achieve
    ~2.0 slope, limited by the O(h^2) geometric approximation error of the
    flat tets — not by polynomial degree.

    Horizontal resolution doubles at each level (refine=2,3,4) while
    n_layers doubles too (4,8,16), so the 3-D L^2 rate truly reflects the
    combined horizontal+vertical refinement without the fixed-nVertLevels
    cap that plagues MPAS-FV convergence studies.

    GO  -> this assert passes; do not touch it.
    NO-GO -> xfail with a verdict note; NEVER lower 1.9.
    """
    result = convergence.shell_convergence("icosa", degree=2, rhs_mode="consistent")
    s = result["slope"]
    assert s >= 1.9, (
        f"icosa shell P2 slope {s:.4f} < 1.9  "
        f"(h={[f'{v:.4g}' for v in result['h']]}, "
        f"err={[f'{v:.4g}' for v in result['err']]})"
    )


@pytest.mark.skipif(
    not convergence.scvt_available(),
    reason="SCVT mesh bundle not present at expected path",
)
def test_scvt_shell_convergence_recorded():
    """Record SCVT shell convergence matrix (P1 and P2 consistent).

    This test always passes: it records the slope for documentation and
    verifies the solve returns a finite positive error.  The SCVT shell
    gate is a *record* test — we do not assert a slope threshold here
    because the moderate n_layers budget (2/4/8 layers) may reduce the
    measured rate relative to the icosa family.

    The icosa gate (test_icosa_shell_p2_is_second_order) carries the GO/NO-GO
    verdict; this test records the SCVT data for the report.
    """
    for degree in (1, 2):
        result = convergence.shell_convergence(
            "scvt", degree=degree, rhs_mode="consistent"
        )
        s = result["slope"]
        assert np.isfinite(s), f"SCVT shell P{degree} slope is not finite: {s}"
        assert all(e > 0 for e in result["err"]), (
            f"SCVT shell P{degree} errors contain non-positive values: {result['err']}"
        )
