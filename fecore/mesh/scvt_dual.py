"""Build a dolfinx 2-sphere surface mesh from an MPAS SCVT grid.nc (dual Delaunay).

Vertices  = SCVT cell centres (xCell/yCell/zCell), snapped to sphere_radius.
Triangles = cellsOnVertex rows (1-based -> 0-based), dropped if any index < 0.

Triangle orientation is corrected so that every face normal points outward
(centroid · normal > 0).  Without this fix the raw MPAS cellsOnVertex has
alternating inward/outward normals, which causes the Euler check to fail
because dolfinx rejects a non-orientable connectivity.

dolfinx 0.10.0 API notes
-------------------------
  create_mesh(comm, cells, element, coords)  -- element BEFORE coords
  ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(3,)))
      describes a 2-D triangle embedded in R^3 (gdim=3).
"""

import numpy as np
import netCDF4 as nc
from mpi4py import MPI
import dolfinx
import ufl
import basix


def sphere_radius(grid_nc_path: str | object) -> float:
    """Return the sphere radius stored in *grid_nc_path* (metres).

    MPAS writes ``sphere_radius`` as a global attribute.  Unit-sphere
    meshes (radius ~ 1.0) use a dimensionless convention; in that case
    we return the standard Earth radius 6 371 229 m so that coordinates
    are scaled to physical metres before mesh creation.
    """
    with nc.Dataset(grid_nc_path) as ds:
        v = getattr(ds, "sphere_radius", None)
        r = float(v) if v is not None else 0.0
    return r if r > 1000.0 else 6_371_229.0


def scvt_surface_mesh(grid_nc_path: str | object) -> dolfinx.mesh.Mesh:
    """Return a dolfinx surface mesh (tdim=2, gdim=3) from *grid_nc_path*.

    Parameters
    ----------
    grid_nc_path:
        Path to an MPAS ``grid.nc`` file that carries the SCVT dual mesh.

    Returns
    -------
    dolfinx.mesh.Mesh
        Triangulated sphere with vertices at the SCVT cell centres.
    """
    R = sphere_radius(grid_nc_path)

    with nc.Dataset(grid_nc_path) as ds:
        xc = np.asarray(ds["xCell"][:])
        yc = np.asarray(ds["yCell"][:])
        zc = np.asarray(ds["zCell"][:])
        # cellsOnVertex is 1-based in MPAS; convert to 0-based
        cov = np.asarray(ds["cellsOnVertex"][:]).astype(np.int64) - 1

    # Stack coordinates and snap all vertices to the target radius.
    pts = np.column_stack((xc, yc, zc)).astype(np.float64)
    pts *= R / np.linalg.norm(pts, axis=1)[:, None]

    # Drop any triangle whose index conversion produced a value < 0
    # (i.e. the original MPAS index was 0, a fill / missing-cell sentinel).
    tris = cov[np.all(cov >= 0, axis=1)].copy()

    # Orient all triangles outward.
    # The raw MPAS cellsOnVertex has alternating CCW/CW triangles (tested
    # empirically: exactly half are inward-facing on the 480 km mesh).
    # dolfinx requires consistent orientation for a manifold mesh.
    _orient_outward(tris, pts)

    # dolfinx 0.10.0: create_mesh(comm, cells, element, coords) — element
    # must come BEFORE coords (the brief snippet has them swapped; the
    # confirmed working signature is the one used here).
    el = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(3,)))
    return dolfinx.mesh.create_mesh(MPI.COMM_WORLD, tris, el, pts)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _orient_outward(tris: np.ndarray, pts: np.ndarray) -> None:
    """Flip triangle vertex order in-place so every face normal is outward.

    A triangle's outward normal is the cross product (v1-v0) × (v2-v0).
    If the dot product of that normal with the face centroid is negative the
    triangle is inward-facing; we correct it by swapping vertices 1 and 2.
    """
    v0 = pts[tris[:, 0]]
    v1 = pts[tris[:, 1]]
    v2 = pts[tris[:, 2]]
    normals = np.cross(v1 - v0, v2 - v0)               # (nTri, 3)
    centroids = (v0 + v1 + v2) / 3.0
    dots = np.einsum("ij,ij->i", normals, centroids)    # (nTri,)
    inward = dots < 0
    tris[inward, 1], tris[inward, 2] = (
        tris[inward, 2].copy(),
        tris[inward, 1].copy(),
    )
