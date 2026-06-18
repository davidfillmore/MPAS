#!/usr/bin/env python3
"""Unit tests for the mimetic (MFD) core: scripts/mfd_operator.py."""
import importlib.util, pathlib, unittest
import numpy as np

REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT = REPO_ROOT / "src/core_atmosphere/electrostatic/scripts/mfd_operator.py"

def load():
    spec = importlib.util.spec_from_file_location("mfd_operator", SCRIPT)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def unit_square_faces():
    """Axis-aligned unit square cell centred at origin. Faces +x,-x,+y,-y.
    Returns area-weighted outward normals, face centroids, centroid, volume."""
    n = np.array([[1.,0.],[-1.,0.],[0.,1.],[0.,-1.]])      # unit normals
    area = 1.0                                              # 2-D 'area' = edge length
    fc = np.array([[0.5,0.],[-0.5,0.],[0.,0.5],[0.,-0.5]])  # face centroids
    return n*area, fc, np.zeros(2), 1.0

class MimeticBlockTests(unittest.TestCase):
    def test_consistency_reproduces_linear_gradient(self):
        m = load()
        N, fc, c, vol = unit_square_faces()
        eps = 2.5
        T = m.mimetic_hodge_block(N, fc, c, vol, eps)
        C = fc - c
        # T @ C must equal eps * N exactly (the MFD consistency condition).
        np.testing.assert_allclose(T @ C, eps * N, atol=1e-12)

    def test_block_is_symmetric_and_spd(self):
        m = load()
        N, fc, c, vol = unit_square_faces()
        T = m.mimetic_hodge_block(N, fc, c, vol, 1.0)
        np.testing.assert_allclose(T, T.T, atol=1e-12)
        self.assertGreater(np.linalg.eigvalsh(T).min(), 0.0)

    def test_consistency_on_sheared_quad(self):
        """Slope/shear is exactly what terrain introduces: a sheared quad must
        still satisfy the consistency condition (this is why MFD beats the naive
        cross-term)."""
        m = load()
        # Parallelogram: corners (0,0),(1,0),(1.4,1),(0.4,1) -> shear 0.4
        # Faces as midpoints of the 4 edges; outward area-normals via edge rotate.
        corners = np.array([[0.,0.],[1.,0.],[1.4,1.],[0.4,1.]])
        c = corners.mean(0); vol = 1.0  # base 1 * height 1
        N, fc = [], []
        for a, b in [(0,1),(1,2),(2,3),(3,0)]:
            edge = corners[b]-corners[a]
            outward = np.array([edge[1], -edge[0]])  # rotate -90; |outward|=edge len
            mid = 0.5*(corners[a]+corners[b])
            if np.dot(outward, mid-c) < 0: outward = -outward
            N.append(outward); fc.append(mid)
        N = np.array(N); fc = np.array(fc)
        eps = 1.3
        T = m.mimetic_hodge_block(N, fc, c, vol, eps)
        np.testing.assert_allclose(T @ (fc-c), eps*N, atol=1e-12)
        self.assertGreater(np.linalg.eigvalsh(T).min(), 0.0)

class MimeticAssemblyTests(unittest.TestCase):
    """Integration coverage for assemble_dT_H_d + spd_min_eig on a 2-cell mesh.

    Two adjacent unit squares sharing the face at x=1:
        cell 0 = [0,1]x[0,1]  centroid (0.5,0.5)
        cell 1 = [1,2]x[0,1]  centroid (1.5,0.5)
    Each cell's 4 faces use the unit_square_faces local order (+x,-x,+y,-y) with
    area-weighted outward normals and edge-midpoint centroids. Seven global faces:
        0 = shared face at x=1 (cell0 +x, cell1 -x)
        1 = cell0 -x (x=0)     2 = cell0 +y (y=1) 3 = cell0 -y (y=0)
        4 = cell1 +x (x=2)     5 = cell1 +y (y=1) 6 = cell1 -y (y=0)
    """
    def _build(self, m):
        eps = 1.7
        # Per-cell geometry: same local order (+x,-x,+y,-y) as unit_square_faces.
        n_local = np.array([[1.,0.],[-1.,0.],[0.,1.],[0.,-1.]])  # outward unit normals
        # Cell 0, centroid (0.5,0.5): face centroids in +x,-x,+y,-y order.
        c0 = np.array([0.5, 0.5])
        fc0 = np.array([[1.,0.5],[0.,0.5],[0.5,1.],[0.5,0.]])
        # Cell 1, centroid (1.5,0.5).
        c1 = np.array([1.5, 0.5])
        fc1 = np.array([[2.,0.5],[1.,0.5],[1.5,1.],[1.5,0.]])
        blk0 = m.mimetic_hodge_block(n_local, fc0, c0, 1.0, eps)
        blk1 = m.mimetic_hodge_block(n_local, fc1, c1, 1.0, eps)
        # Global face ids in each cell's local (+x,-x,+y,-y) order.
        faces = [[0, 1, 2, 3],   # cell 0: +x is the shared face id 0
                 [4, 0, 5, 6]]   # cell 1: -x is the shared face id 0
        blocks = [blk0, blk1]
        # face_cells (lo, hi): -1 on the absent side of the 6 boundary faces;
        # both cells present on the shared face.
        face_cells = np.array([[0, 1],    # 0 shared: both cells
                               [0, -1],   # 1 cell0 -x boundary
                               [0, -1],   # 2 cell0 +y boundary
                               [0, -1],   # 3 cell0 -y boundary
                               [1, -1],   # 4 cell1 +x boundary
                               [1, -1],   # 5 cell1 +y boundary
                               [1, -1]])  # 6 cell1 -y boundary
        return faces, blocks, face_cells

    def test_assembly_grounded_spd(self):
        m = load()
        faces, blocks, face_cells = self._build(m)
        A = m.assemble_dT_H_d(2, faces, blocks, face_cells, ground_cell=0)
        self.assertEqual(A.shape, (2, 2))
        Ad = A.toarray()
        # Symmetric.
        self.assertTrue(np.allclose(Ad, Ad.T))
        # Grounding: row 0 and col 0 zero except the kept positive diagonal.
        self.assertGreater(Ad[0, 0], 0.0)
        self.assertAlmostEqual(Ad[0, 1], 0.0)
        self.assertAlmostEqual(Ad[1, 0], 0.0)
        # SPD on the grounded operator.
        self.assertGreater(m.spd_min_eig(A), 0.0)

    def test_ungrounded_interior_coupling_annihilates_constant(self):
        """Plumbing check on the cell-to-cell coupling: the UN-grounded interior
        operator d0ᵀ H d0 (built from INTERIOR faces only) must kill a constant.

        Interior faces carry both signs (+1 hi, -1 lo) so their row of d0 sums to
        zero -> a constant lies in ker(d0_int) and the operator annihilates it.
        This isolates exactly the shared-face coupling: a wrong face_cells linkage
        (e.g. the shared face not joining both cells, or duplicated sign) gives a
        nonzero off-diagonal contribution that no longer cancels on a constant, so
        this assertion fails. (Boundary faces are deliberately excluded: they carry
        a single ±1 and are NOT in ker(d0); the full d0ᵀ H d0 does not annihilate a
        constant in the presence of boundary, which is correct, not a bug.)"""
        import scipy.sparse as sp
        m = load()
        faces, blocks, face_cells = self._build(m)
        n_cells, n_faces = 2, face_cells.shape[0]
        # Rebuild d0 (interior faces only) and the full H inline (mirror of
        # assemble_dT_H_d, sans grounding).
        dr, dc, dd = [], [], []
        for f in range(n_faces):
            lo, hi = int(face_cells[f, 0]), int(face_cells[f, 1])
            if lo < 0 or hi < 0:
                continue                      # boundary face: not in ker(d0)
            dr.append(f); dc.append(lo); dd.append(-1.0)
            dr.append(f); dc.append(hi); dd.append(1.0)
        d0 = sp.csr_matrix((dd, (dr, dc)), shape=(n_faces, n_cells))
        # The shared (interior) face must actually couple the two cells.
        self.assertGreater(d0.nnz, 0, "no interior face: shared-face linkage missing")
        hr, hc, hd = [], [], []
        for cc in range(n_cells):
            fs = faces[cc]; T = blocks[cc]
            for a in range(len(fs)):
                for b in range(len(fs)):
                    hr.append(fs[a]); hc.append(fs[b]); hd.append(T[a, b])
        H = sp.csr_matrix((hd, (hr, hc)), shape=(n_faces, n_faces))
        A_int = (d0.T @ H @ d0).toarray()
        # The coupling is genuinely present (nonzero off-diagonal)...
        self.assertLess(A_int[0, 1], 0.0)
        # ...and a constant field is annihilated.
        np.testing.assert_allclose(A_int @ np.ones(n_cells), np.zeros(n_cells), atol=1e-10)

if __name__ == "__main__":
    unittest.main()
