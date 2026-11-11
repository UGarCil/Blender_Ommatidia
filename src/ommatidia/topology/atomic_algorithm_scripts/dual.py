"""Atomic: the dual of the centre triangulation = the outlines of the facets.

Each triangle of centres has one cell vertex (its circumcentre): the corner
where its three facets meet. Two triangles that share an edge are joined by a
cell edge: the wall between the two facets on either side of that edge.
"""

from __future__ import annotations

import numpy as np

# DD. CELL_EDGES
# cell_edges = np.ndarray, shape (E, 2), int64
# interp. facet boundaries: an edge joins the cell vertices of two triangles
#         that share a triangulation edge (the dual graph)
# Invariant: i < j; indices < T; no duplicate rows
cell_edges0 = np.zeros((0, 2), dtype=np.int64)  # stub
cell_edges1 = np.array([[0, 1]])


# FD. circumcenters
# Signature: POINTS, TRIANGLES -> POINTS
def circumcenters(points: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Circumcentre of each triangle, in the triangle's plane.

    Algorithm:
        with a = B − A, b = C − A:  A + (|a|² b − |b|² a) × (a × b) / (2 |a × b|²)
    """
    if len(triangles) == 0:
        return np.zeros((0, 3))
    A = points[triangles[:, 0]]
    a = points[triangles[:, 1]] - A
    b = points[triangles[:, 2]] - A
    axb = np.cross(a, b)
    num = np.cross(np.sum(a**2, 1)[:, None] * b - np.sum(b**2, 1)[:, None] * a, axb)
    return A + num / (2.0 * np.sum(axb**2, 1)[:, None])


# FD. circumradii
# Signature: POINTS, TRIANGLES, POINTS (circumcentres) -> DISTANCES
def circumradii(points: np.ndarray, triangles: np.ndarray, centers: np.ndarray) -> np.ndarray:
    """Distance from each circumcentre to its triangle's first corner."""
    if len(triangles) == 0:
        return np.zeros(0)
    return np.linalg.norm(centers - points[triangles[:, 0]], axis=1)


def _shared_edge_table(triangles: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(primal edges shared by two triangles, the two triangle ids), aligned and sorted."""
    if len(triangles) == 0:
        z = np.zeros((0, 2), dtype=np.int64)
        return z, z
    edges = np.sort(triangles[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1)
    owner = np.repeat(np.arange(len(triangles)), 3)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    edges, owner = edges[order], owner[order]
    same = np.all(edges[1:] == edges[:-1], axis=1)
    t1, t2 = owner[:-1][same], owner[1:][same]
    duals = np.column_stack([np.minimum(t1, t2), np.maximum(t1, t2)])
    primal = edges[:-1][same]
    order = np.lexsort((duals[:, 1], duals[:, 0]))
    return primal[order].astype(np.int64), duals[order].astype(np.int64)


# FD. dual_edges
# Signature: TRIANGLES -> CELL_EDGES
def dual_edges(triangles: np.ndarray) -> np.ndarray:
    """Join every two triangles that share an edge.

    Algorithm:
        list (sorted edge, triangle id) for the three edges of each triangle
        sort by edge; consecutive equal edges belong to two adjacent triangles
        emit (smaller id, larger id), sorted
    """
    return _shared_edge_table(triangles)[1]


# FD. crossed_edges
# Signature: TRIANGLES -> (E, 2) int64
def crossed_edges(triangles: np.ndarray) -> np.ndarray:
    """For each row of dual_edges(triangles), the pair of centres whose shared
    triangulation edge that cell edge crosses (i.e. the two facets it separates)."""
    return _shared_edge_table(triangles)[0]
