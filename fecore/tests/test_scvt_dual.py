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


@pytest.mark.skipif(
    not (ROOT / "480km" / "grid.nc").exists(),
    reason="SCVT mesh not present",
)
def test_scvt_surface_mesh_orientation_outward():
    """Assert consistent outward orientation via signed volume.

    For a closed, consistently outward-oriented triangulated surface the
    signed volume V = (1/6) Σ_tri  v0 · (v1 × v2) is positive and close
    to the true sphere volume (4/3)πR³.  A mesh with inconsistent
    orientation (roughly half inward, half outward) yields |V| ≈ 0 because
    the positive and negative tetrahedral contributions cancel.

    The 2562-vertex inscribed triangulation slightly underestimates the
    sphere volume, so we check 0.80*(4/3)πR³ < V < 1.01*(4/3)πR³.
    """
    path = ROOT / "480km" / "grid.nc"
    msh = scvt_dual.scvt_surface_mesh(path)
    R = scvt_dual.sphere_radius(path)

    # Build cell→geometry-node connectivity so we can read triangle coords
    # from the geometry dofmap (for P1 elements this matches topology order).
    msh.topology.create_connectivity(2, 0)
    c2v = msh.topology.connectivity(2, 0)
    geo_dofmap = msh.geometry.dofmap  # shape (n_cells, 3) for P1
    coords = msh.geometry.x           # shape (n_geo_nodes, 3)

    n_cells = msh.topology.index_map(2).size_local
    signed_vol = 0.0
    for i in range(n_cells):
        nodes = geo_dofmap[i]         # 3 geometry node indices for cell i
        v0 = coords[nodes[0]]
        v1 = coords[nodes[1]]
        v2 = coords[nodes[2]]
        signed_vol += np.dot(v0, np.cross(v1, v2))
    signed_vol /= 6.0

    V_sphere = (4.0 / 3.0) * np.pi * R**3

    assert signed_vol > 0, (
        f"Signed volume is NEGATIVE ({signed_vol:.4e}): mesh is consistently "
        f"INWARD-oriented — _orient_outward sign convention is inverted."
    )
    ratio = signed_vol / V_sphere
    assert 0.80 < ratio < 1.01, (
        f"Signed volume ratio {ratio:.4f} out of expected band (0.80, 1.01). "
        f"V_signed={signed_vol:.4e}, V_sphere={V_sphere:.4e}. "
        f"A ratio near 0 means inconsistent orientation; > 1.01 is unexpected."
    )
