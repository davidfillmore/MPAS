#!/usr/bin/env python3
import importlib.util, pathlib, unittest
import numpy as np
REPO = pathlib.Path(__file__).resolve().parents[4]
SCRIPT = REPO / "src/core_atmosphere/electrostatic/scripts/prototype_terrain_3d.py"
def load():
    s = importlib.util.spec_from_file_location("prototype_terrain_3d", SCRIPT)
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m

class TerrainMesh3DTests(unittest.TestCase):
    def test_flat_mesh_has_uniform_layer_thickness(self):
        m = load()
        g = m.terrain_mesh_3d(n_side=8, nz=8, hill_fraction=0.0)
        dz = np.diff(g.zcol, axis=1)            # (nCellsH, nz-1)
        self.assertLess(float(dz.std()), 1e-12)  # flat: identical thickness everywhere

    def test_hill_mesh_varies_thickness_but_positive(self):
        m = load()
        g = m.terrain_mesh_3d(n_side=8, nz=8, hill_fraction=0.3)
        dz = np.diff(g.zcol, axis=1)
        self.assertGreater(float(dz.std()), 0.0)   # terrain: thickness varies
        self.assertGreater(float(dz.min()), 0.0)   # Jacobian stays positive


class TerrainOperator3DTests(unittest.TestCase):
    """B2 operator tests on the STRUCTURED conforming prism->tet mesh.

    INTENT (unchanged from the Delaunay-era tests): the assembled operator is SPD,
    annihilates constants exactly, and its horizontal coupling reach is reported.
    The tet source is now `_structured_prism_tets` (full coverage, no filters).
    """

    def test_operator_spd_on_hill(self):
        m = load(); g = m.terrain_mesh_3d(8, 8, hill_fraction=0.3)
        A = m.terrain_operator(g, ground_cell=0)
        self.assertEqual((abs(A - A.T) > 1e-9).nnz, 0)
        self.assertGreater(m.spd_min_eig_of(A), 0.0)

    def test_ungrounded_annihilates_constants(self):
        m = load(); g = m.terrain_mesh_3d(8, 8, hill_fraction=0.3)
        A = m.terrain_operator_ungrounded(g)
        r = A @ np.ones(A.shape[0])
        self.assertLess(float(np.max(np.abs(r))), 1e-9)

    def test_structured_mesh_has_full_coverage(self):
        """The conforming prism->tet mesh tiles the column extrusion with no gaps
        or overlaps: tet-volume sum == analytic extrusion volume (coverage ~1.0).
        This is the structural fix that replaced the old ~71%-of-box filtered
        Delaunay support."""
        m = load()
        for hf in (0.0, 0.3):
            g = m.terrain_mesh_3d(8, 8, hill_fraction=hf)
            cov = m.coverage_fraction(g)
            self.assertAlmostEqual(cov, 1.0, places=6,
                                   msg=f"coverage {cov} != 1.0 (non-conforming?)")

    def test_halo_width_is_reported(self):
        m = load(); g = m.terrain_mesh_3d(8, 8, hill_fraction=0.3)
        rings = m.halo_rings(g)
        self.assertGreaterEqual(rings, 1)      # records the actual coupling reach
        # Structured prism->tet stencil couples only direct (1-ring) column
        # neighbours in the deep interior; expect exactly 1 ring.
        self.assertEqual(rings, 1)
        print(f"\n[terrain] measured halo rings = {rings}")


class TerrainMMSGateTests(unittest.TestCase):
    """AXIS-1 MMS consistency gate, HONEST pointwise truncation metric.

    Gate metric (no ||phi|| normalization): the finest-two log-log slope of
    `max_i |r_i / b_i|` over interior cells with |b_i| > 1e-3*max|b|, where
    r = A_ung @ phi - b and b = eps*(kx^2+ky^2+kz^2)*phi*V_cell.

    VERDICT: NO-GO.  See notes/2026-06-18-poisson-terrain-vr-verdict.md (rev 2).

    With the genuinely conforming, full-coverage structured prism->tet mesh, the
    pointwise relative truncation max_i|r_i/b_i| still DIVERGES at ~h^-2 (slope
    trending to -2 under refinement) on BOTH the flat and hill domains.  The
    prior GO was spurious: it used a ||phi|| denominator (~h^-3/2) that inflated
    the slope by ~+3.  The structured mesh fixed coverage (71%->100%), SPD,
    exact constant-annihilation, and halo (2->1 ring), and cut the absolute error
    by ~10x, but did NOT make the strong-form pointwise truncation 2nd order:
    P1 finite-element stiffness with lumped mass is 2nd-order in the SOLUTION
    (Galerkin/L2) sense -- the grounded Dirichlet solve converges at clean order
    ~2.05 (documented in the verdict note) -- while its pointwise FV consistency
    error on a sheared terrain mesh is only ~1st order, and the relative pointwise
    measure blows up where b ~ 0.  The honest pointwise gate prescribed for this
    task therefore fails: NO-GO.  These two tests are marked expectedFailure to
    record the NO-GO without forcing a pass by switching normalizations.
    """

    @unittest.expectedFailure
    def test_flat_honesty_guard_is_second_order(self):
        """FLAT honesty guard (pointwise metric).  Expected NO-GO: the pointwise
        slope is ~-1.7 (diverging), not in [1.9, 2.1].  Debugged with a verified
        conforming full-coverage mesh; the divergence is a genuine property of the
        strong-form pointwise truncation of P1+lumped-mass, not a mesh defect.
        Verdict note: notes/2026-06-18-poisson-terrain-vr-verdict.md."""
        m = load()
        slope = m.terrain_mms_order(hill_fraction=0.0)
        self.assertGreater(slope, 1.9)
        self.assertLess(slope, 2.1)

    @unittest.expectedFailure
    def test_terrain_gate_is_second_order_on_hill(self):
        """AXIS-1 GATE (pointwise metric).  GO iff the pointwise truncation slope
        >= 1.9.  VERDICT: NO-GO -- slope ~-1.6 (diverging) on the hill.  Never
        weakened below 1.9.  Verdict note:
        notes/2026-06-18-poisson-terrain-vr-verdict.md."""
        m = load()
        slope = m.terrain_mms_order(hill_fraction=0.3)
        self.assertGreaterEqual(slope, 1.9)


if __name__ == "__main__":
    unittest.main()
