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

    def test_halo_width_is_reported(self):
        m = load(); g = m.terrain_mesh_3d(8, 8, hill_fraction=0.3)
        rings = m.halo_rings(g)
        self.assertGreaterEqual(rings, 1)      # records the actual coupling reach
        print(f"\n[terrain] measured halo rings = {rings}")
