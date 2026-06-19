"""Task 4: Fresh icosahedral surface mesh.

The mesh is built by midpoint subdivision of the regular icosahedron,
projecting all vertices to the unit sphere.  It is a topological 2-sphere
embedded in R^3; the Euler characteristic test guards that property.
"""
import numpy as np
from fecore.mesh import icosa


def test_icosa_surface_is_closed_sphere_and_refines():
    for refine in (0, 1, 2):
        msh = icosa.icosa_surface_mesh(refine)
        nv = msh.topology.index_map(0).size_global
        msh.topology.create_entities(1)
        ne = msh.topology.index_map(1).size_global
        nc = msh.topology.index_map(2).size_global
        assert nv - ne + nc == 2
        r = np.linalg.norm(msh.geometry.x, axis=1)
        assert np.ptp(r) / np.mean(r) < 1e-9          # on the unit sphere
    assert icosa.icosa_ncells(0) == 20
    assert icosa.icosa_ncells(1) == 80
    assert icosa.icosa_ncells(2) > icosa.icosa_ncells(1)


def test_icosa_surface_mesh_orientation_outward():
    """Assert consistent outward orientation via signed volume.

    For a closed, consistently outward-oriented triangulated surface the
    signed volume V = (1/6) Σ_tri  v0 · (v1 × v2) is positive and close
    to the true unit-sphere volume (4/3)π.  A mesh with inconsistent
    orientation (roughly half inward, half outward) yields |V| ≈ 0 because
    the positive and negative tetrahedral contributions cancel.

    A coarse icosphere (refine=2, 320 faces) slightly underestimates the
    sphere volume, so we check 0.85*(4/3)π < V < 1.01*(4/3)π.
    """
    msh = icosa.icosa_surface_mesh(2)

    msh.topology.create_connectivity(2, 0)
    geo_dofmap = msh.geometry.dofmap  # shape (n_cells, 3) for P1
    coords = msh.geometry.x           # shape (n_geo_nodes, 3)

    n_cells = msh.topology.index_map(2).size_local
    signed_vol = 0.0
    for i in range(n_cells):
        nodes = geo_dofmap[i]
        v0 = coords[nodes[0]]
        v1 = coords[nodes[1]]
        v2 = coords[nodes[2]]
        signed_vol += np.dot(v0, np.cross(v1, v2))
    signed_vol /= 6.0

    V_sphere = (4.0 / 3.0) * np.pi  # unit sphere

    assert signed_vol > 0, (
        f"Signed volume is NEGATIVE ({signed_vol:.4e}): mesh is consistently "
        f"INWARD-oriented — _orient_outward sign convention is inverted."
    )
    ratio = signed_vol / V_sphere
    assert 0.85 < ratio < 1.01, (
        f"Signed volume ratio {ratio:.4f} out of expected band (0.85, 1.01). "
        f"V_signed={signed_vol:.4e}, V_sphere={V_sphere:.4e}. "
        f"A ratio near 0 means inconsistent orientation; > 1.01 is unexpected."
    )
