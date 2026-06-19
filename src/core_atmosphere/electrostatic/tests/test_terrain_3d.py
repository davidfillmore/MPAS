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
