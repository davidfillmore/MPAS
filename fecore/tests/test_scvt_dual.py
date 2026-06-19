"""Task 3: SCVT-dual surface mesh from MPAS grid.nc.

The mesh has vertices = SCVT cell centres (xCell/yCell/zCell) and
triangles = cellsOnVertex (1-based -> 0-based).  It is a topological
2-sphere embedded in R^3; the Euler characteristic test guards that
property.
"""
import numpy as np
import pathlib
import pytest

from fecore.mesh import scvt_dual

ROOT = pathlib.Path("~/Data/MPAS/poisson_tier_A2_scvt/meshes").expanduser()


@pytest.mark.skipif(
    not (ROOT / "480km" / "grid.nc").exists(),
    reason="SCVT mesh not present",
)
def test_scvt_surface_mesh_is_a_closed_sphere():
    msh = scvt_dual.scvt_surface_mesh(ROOT / "480km" / "grid.nc")
    assert msh.topology.dim == 2 and msh.geometry.dim == 3
    nv = msh.topology.index_map(0).size_global
    msh.topology.create_entities(1)
    ne = msh.topology.index_map(1).size_global
    nc = msh.topology.index_map(2).size_global
    assert nv - ne + nc == 2  # Euler characteristic of a sphere
    # vertices lie on a sphere (radius ~ constant)
    r = np.linalg.norm(msh.geometry.x, axis=1)
    assert np.ptp(r) / np.mean(r) < 1e-6
