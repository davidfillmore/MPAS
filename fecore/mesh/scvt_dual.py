"""Build dolfinx sphere-surface and spherical-shell meshes from MPAS SCVT grid.nc.

Surface mesh (scvt_surface_mesh):
  Vertices  = SCVT cell centres (xCell/yCell/zCell), snapped to sphere_radius.
  Triangles = cellsOnVertex rows (1-based -> 0-based), dropped if any index < 0.
  Triangle orientation is corrected so that every face normal points outward
  (centroid · normal > 0).

Shell mesh (scvt_shell_mesh, 3-D):
  Extrudes the surface triangulation radially through *n_layers* uniform layers
  from R to R+H.  Each triangular prism is split into 3 tetrahedra via the
  Freudenthal–Kuhn sorted-vertex rule for a conforming mesh.

dolfinx 0.10.0 API notes
-------------------------
  create_mesh(comm, cells, element, coords)  -- element BEFORE coords
  ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(3,)))
      describes a 2-D triangle embedded in R^3 (gdim=3).
  ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
      describes a 3-D tet embedded in R^3 (gdim=3).
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

    Assumption: a stored ``sphere_radius`` ≤ 1000 is treated as a
    unit-sphere convention and replaced with Earth radius (6 371 229 m).
    A genuine sub-kilometre physical mesh would be misread, but this is
    not a concern for MPAS atmospheric grids.
    """
    with nc.Dataset(grid_nc_path) as ds:
        v = getattr(ds, "sphere_radius", None)
        r = float(v) if v is not None else 0.0
    return r if r > 1000.0 else 6_371_229.0


def scvt_shell_mesh(
    grid_nc_path: str | object,
    n_layers: int = 4,
    H: float = 2.0e4,
) -> dolfinx.mesh.Mesh:
    """Return a dolfinx 3-D tetrahedral shell mesh from an MPAS SCVT grid.

    The SCVT-dual surface triangulation is extruded radially outward from R
    (the sphere radius stored in *grid_nc_path*) to R+H through *n_layers*
    uniform layers.  Each triangular prism is split into 3 conforming
    tetrahedra via the Freudenthal–Kuhn sorted-vertex rule.

    Parameters
    ----------
    grid_nc_path:
        Path to an MPAS ``grid.nc`` file.
    n_layers:
        Number of radial layers (prism stacks) in the shell.
    H:
        Shell thickness in metres (default 20 km = 2e4 m).

    Returns
    -------
    dolfinx.mesh.Mesh
        Tetrahedral mesh (tdim=3, gdim=3).
    """
    R = sphere_radius(grid_nc_path)

    with nc.Dataset(grid_nc_path) as ds:
        xc = np.asarray(ds["xCell"][:])
        yc = np.asarray(ds["yCell"][:])
        zc = np.asarray(ds["zCell"][:])
        cov = np.asarray(ds["cellsOnVertex"][:]).astype(np.int64) - 1

    # Surface vertices on the sphere at radius R.
    surf_pts = np.column_stack((xc, yc, zc)).astype(np.float64)
    surf_pts *= R / np.linalg.norm(surf_pts, axis=1)[:, None]

    # Drop triangles with fill/missing indices.
    surf_tris = cov[np.all(cov >= 0, axis=1)].copy()

    # Orient surface consistently outward (required for consistent tet sign).
    _orient_outward(surf_tris, surf_pts)

    n_surf = len(surf_pts)

    # Build 3-D node coordinates: (n_layers+1) shells of surface vertices.
    # Node global index: k * n_surf + v_surf,  k = 0..n_layers, v = 0..n_surf-1.
    # Unit direction vectors for each surface vertex.
    unit_pts = surf_pts / np.linalg.norm(surf_pts, axis=1)[:, None]
    radii = R + np.arange(n_layers + 1) * H / n_layers   # (n_layers+1,)
    coords = (unit_pts[None, :, :] * radii[:, None, None]).reshape(-1, 3)
    coords = np.ascontiguousarray(coords, dtype=np.float64)

    tets = _extrude_to_tets(surf_tris, n_surf, n_layers)
    _fix_tet_orientations(tets, coords)

    el = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    return dolfinx.mesh.create_mesh(MPI.COMM_WORLD, tets, el, coords)


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

def _extrude_to_tets(
    surf_tris: np.ndarray,
    n_surf: int,
    n_layers: int,
) -> np.ndarray:
    """Extrude surface triangles radially into conforming tetrahedra.

    Each surface triangle at radial level k defines a triangular prism with
    bottom nodes at level k and top nodes at level k+1.  The prism is split
    into 3 tetrahedra using the Freudenthal–Kuhn sorted-vertex rule:

      Sort the three prism columns by the surface vertex's global index.
      Let (b0,t0), (b1,t1), (b2,t2) be the sorted bottom/top node pairs.
      Emit exactly:
        (b0, b1, b2, t2)
        (b0, b1, t1, t2)
        (b0, t0, t1, t2)

    Because the sort key is the global surface vertex index, any two adjacent
    prisms that share a vertical quad face agree on the face diagonal, giving a
    conforming tetrahedral mesh.

    Parameters
    ----------
    surf_tris:
        Surface triangle connectivity, shape (n_surf_tris, 3), int.
    n_surf:
        Number of surface vertices.
    n_layers:
        Number of radial layers.

    Returns
    -------
    np.ndarray, shape (3 * n_surf_tris * n_layers, 4), int64
    """
    n_surf_tris = len(surf_tris)
    total_tets = 3 * n_surf_tris * n_layers
    tets = np.empty((total_tets, 4), dtype=np.int64)

    idx = 0
    for k in range(n_layers):
        base_k = k * n_surf
        base_kp1 = (k + 1) * n_surf
        for tri in surf_tris:
            v0, v1, v2 = int(tri[0]), int(tri[1]), int(tri[2])

            # Sort columns by surface vertex global index for conformity.
            col = sorted(
                [(v0, base_k + v0, base_kp1 + v0),
                 (v1, base_k + v1, base_kp1 + v1),
                 (v2, base_k + v2, base_kp1 + v2)],
                key=lambda c: c[0],
            )
            b0, b1, b2 = col[0][1], col[1][1], col[2][1]
            t0, t1, t2 = col[0][2], col[1][2], col[2][2]

            tets[idx + 0] = [b0, b1, b2, t2]
            tets[idx + 1] = [b0, b1, t1, t2]
            tets[idx + 2] = [b0, t0, t1, t2]
            idx += 3

    return tets


def _fix_tet_orientations(tets: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Swap last two vertices of any tet with non-positive signed volume.

    det([v1-v0, v2-v0, v3-v0]) / 6 must be > 0.  Swapping vertices 2 and 3
    negates the determinant without breaking the Freudenthal–Kuhn conformity
    (the shared face between adjacent tets uses vertices 0,1,2 or 0,1,3 which
    are not disturbed by swapping positions 2↔3).

    Parameters
    ----------
    tets:
        Shape (n_tets, 4), int64 — modified in-place.
    coords:
        Shape (n_nodes, 3), float64.

    Returns
    -------
    tets (the same array, modified in-place for convenience).
    """
    v0 = coords[tets[:, 0]]
    v1 = coords[tets[:, 1]]
    v2 = coords[tets[:, 2]]
    v3 = coords[tets[:, 3]]
    a, b, c = v1 - v0, v2 - v0, v3 - v0
    vols = (a[:, 0] * (b[:, 1] * c[:, 2] - b[:, 2] * c[:, 1])
            - a[:, 1] * (b[:, 0] * c[:, 2] - b[:, 2] * c[:, 0])
            + a[:, 2] * (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0])) / 6.0
    neg = vols <= 0
    tets[neg, 2], tets[neg, 3] = tets[neg, 3].copy(), tets[neg, 2].copy()
    return tets


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
