"""Atomic: lay out the vertices and faces of the facet composite (before projection).

Every facet is a fan of rings from its centre to its walls. Walls and corners
are shared between neighbouring facets, so the composite is one connected,
watertight sheet. Heights are not set here (see dome_profile.py): this layout
is built once per Generate, and moving the sliders only moves vertices along
their normals.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# DD. FACET_RESOLUTION
# resolution = (sides: int, rings: int)
# interp. each facet wall is cut into `sides` segments and each facet has
#         `rings` rings from centre to wall (more = smoother domes)
# Invariant: sides >= 1, rings >= 1
# SUBDIVISION_LEVEL L (CONFIG) maps to sides = L + 1, rings = 2·(L + 1); -1 = auto (L = 2)
AUTO_SUBDIVISION = 2


def facet_resolution(subdivision_level: int) -> tuple[int, int]:
    """(sides, rings) for a SUBDIVISION_LEVEL."""
    level = AUTO_SUBDIVISION if subdivision_level < 0 else subdivision_level
    return level + 1, 2 * (level + 1)


# CDD. FACET_LAYOUT
# layout = FacetLayout(points (V,3), owner (V,) int, rho_hex (V,), rho_round (V,),
#                      loop_start (F,) int, loop_vertices (L,) int)
# interp. the unprojected vertices of the composite and its faces (polygon p
#         uses loop_vertices[loop_start[p]:loop_start[p+1]], counter-clockwise
#         about the surface normal). owner = the centre whose facet a vertex is
#         inside (-1 on walls and corners, which are shared). rho_hex = fraction
#         of the way from the centre to the wall along the ray (0 centre, 1 wall);
#         rho_round = distance from the centre / the facet's inner radius.
@dataclass
class FacetLayout:
    points: np.ndarray
    owner: np.ndarray
    rho_hex: np.ndarray
    rho_round: np.ndarray
    loop_start: np.ndarray
    loop_vertices: np.ndarray


# FD. facet_layout
# Signature: POINTS (centres), POINTS (cell points), CELL_LOOPS, CELL_EDGES, FACET_RESOLUTION
#            -> FACET_LAYOUT, (N,) float (inner radius per centre; 0 without a facet)
def facet_layout(centers: np.ndarray, cell_points: np.ndarray, loops: list,
                 cell_edges: np.ndarray, sides: int, rings: int) -> tuple[FacetLayout, np.ndarray]:
    """Vertices: used corners, then wall points, then each facet's centre and rings.

    Algorithm:
        facets = loops of >= 3 corners whose consecutive corners all share a wall
            (overlapping triangles at hard creases can break that: no facet there)
        used corners and used walls (consecutive corners of a loop share a wall)
        wall points: sides−1 points evenly along each used wall (lower corner id first)
        for each group of facets with the same number of corners m (vectorised):
            boundary = per side: its first corner, then its wall points in loop order
            inner radius = mean distance from the centre to the wall midpoints
            ring r (1..rings−1) = centre + r/rings · (boundary − centre)
            faces: centre fan to ring 1, quads between rings, last ring to boundary
    """
    edge_row = {(int(a), int(b)): e for e, (a, b) in enumerate(np.asarray(cell_edges).tolist())}
    built = [k for k, loop in enumerate(loops) if len(loop) >= 3 and all(
        _key(loop[j], loop[(j + 1) % len(loop)]) in edge_row for j in range(len(loop)))]
    radius = np.zeros(len(centers))
    if not built:
        empty_i = np.zeros(0, dtype=np.int64)
        return FacetLayout(np.zeros((0, 3)), empty_i, np.zeros(0), np.zeros(0), empty_i, empty_i), radius

    used_corners = np.unique(np.concatenate([loops[k] for k in built]))
    used_walls = sorted({edge_row[_key(loop[j], loop[(j + 1) % len(loop)])]
                         for k in built for loop in [loops[k]] for j in range(len(loop))})

    corner_vid = np.full(len(cell_points), -1, dtype=np.int64)
    corner_vid[used_corners] = np.arange(len(used_corners))
    wall_vid = np.full((len(cell_edges), max(sides - 1, 0)), -1, dtype=np.int64)
    n_walls = len(used_walls)
    if sides > 1 and n_walls:
        wall_vid[used_walls] = len(used_corners) + np.arange(n_walls * (sides - 1)).reshape(n_walls, sides - 1)
    walls = np.asarray(cell_edges)[used_walls]
    f = (np.arange(1, sides) / sides)[None, :, None]
    wall_points = cell_points[walls[:, 0], None] * (1 - f) + cell_points[walls[:, 1], None] * f

    points = [cell_points[used_corners], wall_points.reshape(-1, 3)]
    owner = [np.full(len(used_corners) + n_walls * (sides - 1), -1, dtype=np.int64)]
    rho_hex = [np.ones(len(owner[0]))]
    rho_round = [np.ones(len(owner[0]))]
    faces = []
    next_vid = len(owner[0])

    shared_points = np.concatenate(points)        # corners and wall points
    by_size = {}
    for k in built:
        by_size.setdefault(len(loops[k]), []).append(k)
    for m, cells in sorted(by_size.items()):
        cells = np.asarray(cells)
        loop = np.stack([loops[k] for k in cells])                 # (G, m)
        g = len(cells)
        a, b = loop, np.roll(loop, -1, axis=1)
        rows = np.vectorize(lambda x, y: edge_row[_key(x, y)])(a, b) if g else np.zeros((0, m), int)
        forward = a < b
        # boundary ids and positions, (G, m*sides)
        ids = np.empty((g, m, sides), dtype=np.int64)
        ids[:, :, 0] = corner_vid[a]
        if sides > 1:
            w = wall_vid[rows]                                      # (G, m, sides-1)
            ids[:, :, 1:] = np.where(forward[..., None], w, w[..., ::-1])
        boundary_ids = ids.reshape(g, m * sides)
        boundary = shared_points[boundary_ids]                       # (G, M, 3)
        center = centers[cells]                                      # (G, 3)
        mids = 0.5 * (cell_points[a] + cell_points[b])
        r_in = np.linalg.norm(mids - center[:, None], axis=2).mean(axis=1)
        radius[cells] = r_in

        M = m * sides
        center_vid = next_vid + np.arange(g)
        next_vid += g
        ring_vid = next_vid + np.arange(g * (rings - 1) * M).reshape(g, rings - 1, M)
        next_vid += g * (rings - 1) * M
        fr = (np.arange(1, rings) / rings)[None, :, None, None]
        ring_points = center[:, None, None] + fr * (boundary[:, None] - center[:, None, None])
        dist = fr[..., 0] * np.linalg.norm(boundary - center[:, None], axis=2)[:, None]   # (G, R-1, M)

        points += [center, ring_points.reshape(-1, 3)]
        owner += [cells, np.repeat(cells, (rings - 1) * M)]
        rho_hex += [np.zeros(g), np.broadcast_to(fr[..., 0], (g, rings - 1, M)).ravel()]
        rho_round += [np.zeros(g), (dist / r_in[:, None, None]).ravel()]

        # rings as id arrays, innermost first, boundary last: (G, rings, M)
        levels = np.concatenate([ring_vid, boundary_ids[:, None]], axis=1)
        nxt = np.roll(np.arange(M), -1)
        fan = np.stack([np.repeat(center_vid, M), levels[:, 0].ravel(), levels[:, 0][:, nxt].ravel()], axis=1)
        faces.append(fan)
        if rings > 1:
            inner, outer = levels[:, :-1], levels[:, 1:]
            quads = np.stack([inner, outer, outer[..., nxt], inner[..., nxt]], axis=-1).reshape(-1, 4)
            faces.append(quads)

    tris = [x for x in faces if x.shape[1] == 3]
    quads = [x for x in faces if x.shape[1] == 4]
    tris = np.concatenate(tris) if tris else np.zeros((0, 3), dtype=np.int64)
    quads = np.concatenate(quads) if quads else np.zeros((0, 4), dtype=np.int64)
    loop_vertices = np.concatenate([tris.ravel(), quads.ravel()]).astype(np.int64)
    loop_start = np.concatenate([np.arange(len(tris)) * 3,
                                 3 * len(tris) + np.arange(len(quads)) * 4]).astype(np.int64)
    layout = FacetLayout(np.concatenate(points), np.concatenate(owner), np.concatenate(rho_hex),
                         np.concatenate(rho_round), loop_start, loop_vertices)
    return layout, radius


def _key(a, b) -> tuple[int, int]:
    a, b = int(a), int(b)
    return (a, b) if a < b else (b, a)
