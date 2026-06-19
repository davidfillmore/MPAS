"""Task 4: Fresh icosahedral surface mesh.

The mesh is built by midpoint subdivision of the regular icosahedron,
projecting all vertices to the unit sphere.  It is a topological 2-sphere
embedded in R^3; the Euler characteristic test guards that property.
"""
import numpy as np
from fecore.mesh import icosa


def test_icosa_surface_is_closed_sphere_and_refines():
    for refine in (1, 2):
        msh = icosa.icosa_surface_mesh(refine)
        nv = msh.topology.index_map(0).size_global
        msh.topology.create_entities(1)
        ne = msh.topology.index_map(1).size_global
        nc = msh.topology.index_map(2).size_global
        assert nv - ne + nc == 2
        r = np.linalg.norm(msh.geometry.x, axis=1)
        assert np.ptp(r) / np.mean(r) < 1e-9          # on the unit sphere
    assert icosa.icosa_ncells(2) > icosa.icosa_ncells(1)
