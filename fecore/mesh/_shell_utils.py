"""Shared helpers for spherical-shell mesh extrusion (Freudenthal–Kuhn split).

Used by both ``fecore.mesh.icosa`` and ``fecore.mesh.scvt_dual``.
"""

import numpy as np


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


def _fix_tet_orientations(tets: np.ndarray, coords: np.ndarray) -> tuple[np.ndarray, int]:
    """Swap last two vertices of any tet with non-positive signed volume.

    det([v1-v0, v2-v0, v3-v0]) / 6 must be > 0.  Swapping vertices 2 and 3
    negates the determinant without breaking the Freudenthal–Kuhn conformity:
    swapping positions 2↔3 changes only orientation (sign), not the tet's
    vertex SET, so every facet's node-set — hence inter-cell facet
    sharing/conformity — is preserved; only the signed volume flips to positive.

    In practice roughly half of all tets (those whose columns are wound such
    that the outward-sphere radial direction gives a negative determinant for
    the FK template) require flipping.

    Parameters
    ----------
    tets:
        Shape (n_tets, 4), int64 — modified in-place.
    coords:
        Shape (n_nodes, 3), float64.

    Returns
    -------
    tets:
        The same array, modified in-place.
    n_flipped:
        Number of tetrahedra whose orientation was corrected.
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
    n_flipped = int(np.sum(neg))
    tets[neg, 2], tets[neg, 3] = tets[neg, 3].copy(), tets[neg, 2].copy()
    return tets, n_flipped
