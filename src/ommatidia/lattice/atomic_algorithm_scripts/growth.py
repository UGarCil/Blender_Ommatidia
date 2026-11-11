"""Atomic: grow a perfect hexagonal patch outward from a centre, ring by ring.

A patch is a hexagonal lattice indexed by axial coordinates (q, r). Each new
node is placed at the apex of an equilateral triangle built on two placed
neighbours, so every facet wall starts exactly one SPACING long. Curvature
still compresses the patch; spring steps on the patch edges spread that strain
instead of letting it turn into 5/7 defects. Where the patch would overlap
itself (fronts meeting on the far side of a curved surface) or leave the
surface, the node dies and the gap becomes a seam for the compound to fill.
"""

from __future__ import annotations

import numpy as np

from .neighbors import _candidate_pairs, neighbor_pairs

# Axial neighbour offsets, counter-clockwise; consecutive entries are
# neighbours of each other.
HEX_DIRECTIONS = np.array([(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)])
_DIRECTION_SET = {tuple(d) for d in HEX_DIRECTIONS}
_SQRT3_2 = np.sqrt(3.0) / 2.0

COLLISION_DISTANCE = 0.7   # x spacing; a node this close to a non-neighbour overlaps the patch
# Lattice neighbours may come closer: curvature compresses the newest ring
# (~30% at ring 8 on a 750-sphere) until the next spring steps even it out.
# Measured on the 750-sphere: 0.7 -> perfect to 49°, 0.55 -> 62° (logbook).
NEIGHBOR_COLLISION_DISTANCE = 0.55  # x spacing
OFF_SURFACE_DISTANCE = 0.5  # x spacing; a proposal this far from the surface fell off it
RIM_DISTANCE = 0.5          # x spacing; same rule as trim_boundary
# A new node whose normal turns more than ~95° from its parents' jumped to another
# surface (the other skin of a thin solid faces ~180° away). 90° creases (a raw
# cube) still pass: a stricter 75° limit dropped the cube from 85% to 72% hexagons.
MIN_TURN_COSINE = -0.1
SPRING_STEP = 0.25          # fraction of the spring force applied per step
SPRING_MOMENTUM = 0.8

# DD. HEX_INDEX
# hex_index = dict[(q: int, r: int) -> int]
# interp. the nodes of one patch in axial coordinates; the value is the node's
#         row in the point arrays, or DEAD if it could not be placed
# Invariant: (0, 0) is the patch centre; ring(q, r) = hex distance from (0, 0)
DEAD = -1
hex_index0 = {}                       # stub
hex_index1 = {(0, 0): 0, (1, 0): 1, (0, 1): DEAD}

# DD. PATCH_EDGES
# patch_edges = np.ndarray, shape (E, 2), int64
# interp. pairs of point rows that are lattice neighbours in a patch (facet walls
#         of the perfect lattice); each pair appears once
patch_edges0 = np.zeros((0, 2), dtype=np.int64)  # stub


# FD. hex_ring
# Signature: HEX_INDEX key -> int
def hex_ring(node: tuple[int, int]) -> int:
    """Hex distance of an axial node from the patch centre."""
    q, r = node
    return (abs(q) + abs(r) + abs(q + r)) // 2


# Facing direction used to pick the patch centre on a closed mesh: Blender's
# front view looks along +Y, so the visible side faces -Y.
FRONT_DIRECTION = np.array([0.0, -1.0, 0.0])


# FD. auto_origin
# Signature: POINTS (surface samples), DISTANCES (boundary) -> POINT
def auto_origin(samples: np.ndarray, boundary_distance: np.ndarray) -> np.ndarray:
    """Default patch centre: the sample farthest from an open edge; on a closed
    surface (all distances infinite), the sample farthest toward FRONT_DIRECTION
    from the centroid."""
    boundary_distance = np.asarray(boundary_distance)
    if np.isfinite(boundary_distance).any():
        return samples[int(np.argmax(np.where(np.isfinite(boundary_distance), boundary_distance, -1)))]
    return samples[int(np.argmax((samples - samples.mean(axis=0)) @ FRONT_DIRECTION))]


# FD. tangent_frame
# Signature: NORMAL, float (degrees) -> (3,), (3,)
def tangent_frame(normal: np.ndarray, rotation_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Two orthonormal tangent directions at the patch centre.

    Algorithm:
        reference = world +Z projected onto the tangent plane (world +X if the
            normal is nearly vertical)
        rotate it about the normal by rotation_deg -> e1; e2 = normal × e1
    """
    n = normal / np.linalg.norm(normal)
    ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = ref - (ref @ n) * n
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    a = np.radians(rotation_deg)
    return np.cos(a) * e1 + np.sin(a) * e2, -np.sin(a) * e1 + np.cos(a) * e2


# FD. first_ring
# Signature: POINT, (3,), (3,), SPACING -> POINTS (6, 3)
def first_ring(center: np.ndarray, e1: np.ndarray, e2: np.ndarray, spacing: float) -> np.ndarray:
    """Unprojected positions of the six neighbours of the centre, in HEX_DIRECTIONS order."""
    d = HEX_DIRECTIONS
    flat = np.column_stack([d[:, 0] + 0.5 * d[:, 1], _SQRT3_2 * d[:, 1]])
    return center + spacing * (flat[:, :1] * e1 + flat[:, 1:] * e2)


# FD. ring_candidates
# Signature: [HEX_INDEX key], HEX_INDEX, int -> [HEX_INDEX key]
def ring_candidates(front: list, index: dict, ring: int) -> list:
    """Unvisited nodes of the given ring next to a live node of the front.

    Algorithm:
        for each live front node, for each of its 6 neighbours:
            keep it if it is in `ring` and not yet in the index (in first-seen order)
    """
    seen = {}
    for q, r in front:
        if index.get((q, r), DEAD) == DEAD:
            continue
        for dq, dr in HEX_DIRECTIONS:
            node = (q + dq, r + dr)
            if node not in index and node not in seen and hex_ring(node) == ring:
                seen[node] = None
    return list(seen)


# FD. apex_triples
# Signature: [HEX_INDEX key], HEX_INDEX -> (M,), (M,), (M,), (M,) int64
def apex_triples(nodes: list, index: dict) -> tuple[np.ndarray, ...]:
    """For each node, every pair (a, b) of placed neighbours that are neighbours of
    each other, with c = their other common neighbour (also placed).

    Algorithm:
        for each node j and each consecutive direction pair (d1, d2):
            a = node + d1, b = node + d2, c = node + d1 + d2
            if a, b and c are all placed: emit (j, a, b, c)
    """
    owner, a_rows, b_rows, c_rows = [], [], [], []
    for j, (q, r) in enumerate(nodes):
        for t in range(6):
            d1, d2 = HEX_DIRECTIONS[t], HEX_DIRECTIONS[(t + 1) % 6]
            a = index.get((q + d1[0], r + d1[1]), DEAD)
            b = index.get((q + d2[0], r + d2[1]), DEAD)
            c = index.get((q + d1[0] + d2[0], r + d1[1] + d2[1]), DEAD)
            if a != DEAD and b != DEAD and c != DEAD:
                owner.append(j)
                a_rows.append(a)
                b_rows.append(b)
                c_rows.append(c)
    as_int = lambda x: np.asarray(x, dtype=np.int64)  # noqa: E731
    return as_int(owner), as_int(a_rows), as_int(b_rows), as_int(c_rows)


# FD. apex_points
# Signature: POINTS, NORMALS, rows a, rows b, rows c, SPACING -> POINTS (M, 3)
def apex_points(points: np.ndarray, normals: np.ndarray, a: np.ndarray, b: np.ndarray,
                c: np.ndarray, spacing: float) -> np.ndarray:
    """Apex of the equilateral triangle on edge (a, b), on the side away from c.

    Algorithm:
        m = midpoint of a, b; tangent plane normal = mean normal of a, b
        u = unit vector in that plane, perpendicular to (b − a), pointing away from c
        h = sqrt(spacing² − |b − a|²/4)   (0 if the edge is longer than 2·spacing)
        return m + h·u
    """
    m = 0.5 * (points[a] + points[b])
    nm = normals[a] + normals[b]
    ab = points[b] - points[a]
    u = np.cross(nm, ab)
    u /= np.maximum(np.linalg.norm(u, axis=1, keepdims=True), 1e-300)
    u *= np.sign(np.sum((m - points[c]) * u, axis=1))[:, None]
    half = 0.5 * np.linalg.norm(ab, axis=1)
    h = np.sqrt(np.maximum(spacing**2 - half**2, 0.0))
    return m + h[:, None] * u


# FD. mean_by_owner
# Signature: POINTS (M, 3), (M,) int, int -> POINTS (K, 3), (K,) int
def mean_by_owner(proposals: np.ndarray, owner: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Average the proposals of each owner; return the averages and the owners that had any."""
    count = np.bincount(owner, minlength=n)
    has = np.flatnonzero(count)
    total = np.column_stack([np.bincount(owner, proposals[:, k], n) for k in range(3)])
    return total[has] / count[has, None], has


# FD. placement_mask
# DD. LATTICE_LINKS
# lattice_links = np.ndarray, shape (L, 2), int64
# interp. (new node, other point) pairs that are lattice neighbours, so they may
#         be closer than COLLISION_DISTANCE (down to NEIGHBOR_COLLISION_DISTANCE)
lattice_links0 = np.zeros((0, 2), dtype=np.int64)  # stub


# FD. placement_mask
# Signature: POINTS (new, projected), POINTS (proposals), NORMALS, DISTANCES (boundary),
#            POINTS (existing), NORMALS (existing), SPACING,
#            LATTICE_LINKS (new -> existing), LATTICE_LINKS (new -> new) -> (K,) bool
def placement_mask(new: np.ndarray, proposals: np.ndarray, new_normals: np.ndarray,
                   boundary_distance: np.ndarray, existing: np.ndarray,
                   existing_normals: np.ndarray, spacing: float,
                   links: np.ndarray = lattice_links0,
                   batch_links: np.ndarray = lattice_links0,
                   parent_normals: np.ndarray | None = None) -> np.ndarray:
    """Which new nodes may be placed (the rest die and leave a seam).

    Algorithm:
        die if closer than RIM_DISTANCE·spacing to an open edge
        die if the projection moved the proposal more than OFF_SURFACE_DISTANCE·spacing
        die if its normal turned away from its parents' (cosine < MIN_TURN_COSINE)
        die if an existing point facing the same side is closer than
            COLLISION_DISTANCE·spacing (NEIGHBOR_COLLISION_DISTANCE for lattice neighbours)
        among the survivors, greedily, by the same rule against earlier survivors
    """
    alive = np.asarray(boundary_distance) >= RIM_DISTANCE * spacing
    alive &= np.linalg.norm(new - proposals, axis=1) <= OFF_SURFACE_DISTANCE * spacing
    if parent_normals is not None:
        unit = parent_normals / np.maximum(np.linalg.norm(parent_normals, axis=1, keepdims=True), 1e-300)
        alive &= np.sum(new_normals * unit, axis=1) >= MIN_TURN_COSINE
    radius = COLLISION_DISTANCE * spacing
    i, j = _candidate_pairs(new, existing, radius)
    hit = _collides(new[i], existing[j], new_normals[i], existing_normals[j],
                    _is_link(i, j, links, len(existing)), spacing)
    alive[i[hit]] = False
    pairs = neighbor_pairs(new, radius)
    if len(pairs):
        p, q = pairs[:, 0], pairs[:, 1]
        linked = _is_link(p, q, batch_links, len(new)) | _is_link(q, p, batch_links, len(new))
        hit = _collides(new[p], new[q], new_normals[p], new_normals[q], linked, spacing)
        for a, b in pairs[hit]:
            if alive[a] and alive[b]:
                alive[b] = False
    return alive


def _is_link(i, j, links, n_other):
    """Which (i, j) pairs appear in links."""
    if len(links) == 0 or len(i) == 0:
        return np.zeros(len(i), dtype=bool)
    return np.isin(i * (n_other + 1) + j, links[:, 0] * (n_other + 1) + links[:, 1])


def _collides(a, b, na, nb, linked, spacing):
    d2 = np.sum((a - b) ** 2, axis=1)
    limit = np.where(linked, NEIGHBOR_COLLISION_DISTANCE, COLLISION_DISTANCE) * spacing
    return (d2 < limit**2) & (np.sum(na * nb, axis=1) > 0)


# FD. node_links
# Signature: [HEX_INDEX key], HEX_INDEX, int -> LATTICE_LINKS, LATTICE_LINKS
def node_links(nodes: list, index: dict, row_offset: int) -> tuple[np.ndarray, np.ndarray]:
    """Lattice-neighbour links of new nodes: to placed rows (shifted by row_offset)
    and to other new nodes of the same batch.

    Algorithm:
        for each new node j and each of its 6 neighbours:
            placed -> link (j, row + row_offset); in the batch -> link (j, its batch number)
    """
    batch = {node: k for k, node in enumerate(nodes)}
    placed, among = [], []
    for j, (q, r) in enumerate(nodes):
        for dq, dr in HEX_DIRECTIONS:
            other = (q + dq, r + dr)
            row = index.get(other, DEAD)
            if row != DEAD:
                placed.append((j, row + row_offset))
            elif other in batch:
                among.append((j, batch[other]))
    as_links = lambda x: np.asarray(x, dtype=np.int64).reshape(-1, 2)  # noqa: E731
    return as_links(placed), as_links(among)


# FD. patch_edges
# Signature: HEX_INDEX -> PATCH_EDGES
def patch_edges(index: dict) -> np.ndarray:
    """Every pair of placed lattice neighbours, once each.

    Algorithm:
        for each placed node, look at the first three directions only (the other
        three are the same edges seen from the other end); keep placed pairs
    """
    edges = []
    for (q, r), row in index.items():
        if row == DEAD:
            continue
        for dq, dr in HEX_DIRECTIONS[:3]:
            other = index.get((q + dq, r + dr), DEAD)
            if other != DEAD:
                edges.append((row, other))
    return np.asarray(edges, dtype=np.int64).reshape(-1, 2)


# FD. spring_step
# Signature: POINTS, VELOCITIES, NORMALS, PATCH_EDGES, (N,) bool, SPACING -> POINTS, VELOCITIES
def spring_step(points: np.ndarray, velocities: np.ndarray, normals: np.ndarray,
                edges: np.ndarray, movable: np.ndarray, spacing: float) -> tuple[np.ndarray, np.ndarray]:
    """Pull every patch edge toward length SPACING, along the surface.

    The edges fix the topology, so this evens out the strain of a curved patch
    without creating defects (plain repulsion would let 5/7 pairs in).

    Algorithm:
        for each edge (i, j): f = (|p_i − p_j| − spacing) · unit(p_i − p_j)
            force on i −= f; force on j += f
        move = MOMENTUM · (tangential previous move) + STEP · (tangential force)
        fixed points do not move
        return points + move, move
    """
    force = np.zeros_like(points)
    if len(edges):
        i, j = edges[:, 0], edges[:, 1]
        diff = points[i] - points[j]
        length = np.linalg.norm(diff, axis=1)
        f = diff * ((length - spacing) / np.maximum(length, 1e-300))[:, None]
        for k in range(3):
            force[:, k] = np.bincount(j, f[:, k], len(points)) - np.bincount(i, f[:, k], len(points))
    tangential = lambda v: v - np.sum(v * normals, axis=1, keepdims=True) * normals  # noqa: E731
    move = SPRING_MOMENTUM * tangential(velocities) + SPRING_STEP * tangential(force)
    move[~movable] = 0.0
    return points + move, move
