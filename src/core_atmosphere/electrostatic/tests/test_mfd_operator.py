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

if __name__ == "__main__":
    unittest.main()
