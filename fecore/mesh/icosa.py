"""Build a dolfinx 2-sphere surface mesh via icosahedral midpoint subdivision.

Starting from the regular icosahedron (12 vertices, 20 faces), each level of
refinement splits every triangle into four by inserting the midpoint of each
edge and projecting it to the unit sphere.  A shared-edge cache (keyed by
sorted vertex-index pair) ensures each midpoint vertex is created exactly once.

dolfinx 0.10.0 API notes
-------------------------
  create_mesh(comm, cells, element, coords)  -- element BEFORE coords
  ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(3,)))
      describes a 2-D triangle embedded in R^3 (gdim=3).
"""

import numpy as np
from mpi4py import MPI
import dolfinx
import ufl
import basix

# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------


def icosa_surface_mesh(refine: int) -> dolfinx.mesh.Mesh:
    """Return a dolfinx surface mesh of the unit sphere (tdim=2, gdim=3).

    Parameters
    ----------
    refine:
        Number of midpoint-subdivision passes.  Level 0 is the raw 20-face
        icosahedron; each level quadruples the face count.

    Returns
    -------
    dolfinx.mesh.Mesh
        Triangulated unit sphere with 20*4**refine faces.
    """
    pts, tris = _subdivided_icosa(refine)
    el = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(3,)))
    return dolfinx.mesh.create_mesh(MPI.COMM_WORLD, tris, el, pts)


def icosa_ncells(refine: int) -> int:
    """Return the number of triangular faces at subdivision level *refine*."""
    return 20 * (4 ** refine)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _base_icosahedron():
    """Return (vertices, faces) for the regular icosahedron on the unit sphere.

    Vertices are the 12 cyclic permutations of (0, ±1, ±φ) with φ=(1+√5)/2,
    normalized to the unit sphere.  All 20 faces are wound consistently
    outward (centroid · normal > 0).
    """
    phi = (1.0 + np.sqrt(5.0)) / 2.0

    raw = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            raw.append([0.0,     s1,       s2 * phi])
            raw.append([s1,      s2 * phi, 0.0     ])
            raw.append([s2 * phi, 0.0,     s1      ])
    verts = np.array(raw, dtype=np.float64)
    verts /= np.linalg.norm(verts, axis=1)[:, None]   # project to unit sphere

    # Derive faces: all vertex triples whose pairwise chord-lengths equal the
    # icosahedron edge length (the shortest nonzero inter-vertex distance).
    dists = np.linalg.norm(verts[None, :, :] - verts[:, None, :], axis=2)
    np.fill_diagonal(dists, np.inf)
    edge_len = dists.min()
    adj = np.abs(dists - edge_len) < 1e-10

    n = len(verts)
    faces_list = []
    for i in range(n):
        for j in range(i + 1, n):
            if not adj[i, j]:
                continue
            for k in range(j + 1, n):
                if adj[i, k] and adj[j, k]:
                    faces_list.append([i, j, k])
    faces = np.array(faces_list, dtype=np.int64)

    # Orient all faces outward.
    _orient_outward(faces, verts)
    return verts, faces


def _midpoint_on_sphere(i: int, j: int,
                        pts: list, cache: dict) -> int:
    """Return the index of the unit-sphere midpoint of edge (i, j).

    If the midpoint has already been created (shared edge), return its cached
    index; otherwise create it, append to *pts*, cache, and return the new index.
    """
    key = (min(i, j), max(i, j))
    if key in cache:
        return cache[key]
    mid = (pts[i] + pts[j]) * 0.5
    mid /= np.linalg.norm(mid)
    idx = len(pts)
    pts.append(mid)
    cache[key] = idx
    return idx


def _subdivided_icosa(levels: int):
    """Return (vertices, faces) for *levels* rounds of midpoint subdivision."""
    base_verts, base_faces = _base_icosahedron()

    # Work with a list for dynamic vertex appending.
    pts = list(base_verts)
    tris = list(base_faces)

    for _ in range(levels):
        new_tris = []
        cache: dict = {}
        for tri in tris:
            a, b, c = int(tri[0]), int(tri[1]), int(tri[2])
            ab = _midpoint_on_sphere(a, b, pts, cache)
            bc = _midpoint_on_sphere(b, c, pts, cache)
            ac = _midpoint_on_sphere(a, c, pts, cache)
            # Each triangle → four child triangles (orientation preserved).
            new_tris.extend([
                [a,  ab, ac],
                [b,  bc, ab],
                [c,  ac, bc],
                [ab, bc, ac],
            ])
        tris = new_tris

    return (np.array(pts, dtype=np.float64),
            np.array(tris, dtype=np.int64))


def _orient_outward(tris: np.ndarray, pts: np.ndarray) -> None:
    """Flip triangle vertex order in-place so every face normal is outward.

    A triangle's outward normal is (v1-v0) × (v2-v0).  If the dot product
    with the face centroid is negative the triangle is inward-facing; we
    correct it by swapping vertices 1 and 2.
    """
    v0 = pts[tris[:, 0]]
    v1 = pts[tris[:, 1]]
    v2 = pts[tris[:, 2]]
    normals = np.cross(v1 - v0, v2 - v0)
    centroids = (v0 + v1 + v2) / 3.0
    dots = np.einsum("ij,ij->i", normals, centroids)
    inward = dots < 0
    tris[inward, 1], tris[inward, 2] = (
        tris[inward, 2].copy(),
        tris[inward, 1].copy(),
    )
