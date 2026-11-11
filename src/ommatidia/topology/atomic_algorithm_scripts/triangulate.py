"""Atomic: triangulate ommatidium centres on a surface without SciPy.

Each centre builds its own Delaunay fan in its tangent plane from nearby
centres; a triangle is kept only if all three of its corners proposed it.
On a relaxed lattice the fans agree almost everywhere, and the few that do
not (near-co-circular quadruples) show up as small holes, which we count.
"""

from __future__ import annotations

import numpy as np

# DD. TRIANGLES
# triangles = np.ndarray, shape (T, 3), int64
# interp. triangulation of the LATTICE centres; rows index LATTICE.centers, and
#         vertex order is counter-clockwise seen from outside (along the normal)
# Invariant: every edge is shared by at most 2 triangles; no triangle has a
#            circumradius above SPACING
triangles0 = np.zeros((0, 3), dtype=np.int64)  # stub
triangles1 = np.array([[0, 1, 2], [0, 2, 3]])

# DD. VALENCE
# valence = np.ndarray, shape (N,), int64
# interp. number of triangulation neighbours of each centre = sides of its cell;
#         6 = hexagon, 5 and 7 are lattice defects
valence0 = np.zeros(0, dtype=np.int64)  # stub
valence1 = np.array([6, 6, 5, 7, 6])

_CHUNK = 256  # centres per vectorised block; bounds memory at ~C·m³ floats


def _adjacency(pairs: np.ndarray, n: int) -> np.ndarray:
    """Padded neighbour table (n, m_max), -1 where empty."""
    both = np.concatenate([pairs, pairs[:, ::-1]]) if len(pairs) else np.zeros((0, 2), np.int64)
    both = both[np.argsort(both[:, 0], kind="stable")]
    degree = np.bincount(both[:, 0], minlength=n)
    table = np.full((n, max(int(degree.max(initial=0)), 2)), -1, dtype=np.int64)
    slot = np.arange(len(both)) - np.repeat(np.cumsum(degree) - degree, degree)
    table[both[:, 0], slot] = both[:, 1]
    return table


def _tangent_basis(normals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    helper = np.where(np.abs(normals[:, :1]) < 0.9, [[1.0, 0, 0]], [[0, 1.0, 0]])
    u = np.cross(normals, helper)
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    return u, np.cross(normals, u)


# FD. local_delaunay_fans
# Signature: POINTS, NORMALS, NEIGHBOR_PAIRS, SPACING -> (F, 3) int64 (proposed triangles)
def local_delaunay_fans(points: np.ndarray, normals: np.ndarray, pairs: np.ndarray,
                        max_circumradius: float) -> np.ndarray:
    """For each centre i, propose every triangle (i, j, k) of its local Delaunay fan.

    The empty-circle test of a triangle is done in one plane that depends only
    on the triangle: the plane perpendicular to the mean surface normal of its
    three corners. All three corners therefore reach the same verdict (as long
    as each sees every point near the circle, which holds when the neighbour
    radius is >= 2 · max_circumradius), so fans agree. On smooth surfaces this
    plane is ~the triangle's own plane; across a sharp crease it is the
    bisecting plane, which "unfolds" the crease symmetrically.

    Algorithm:
        for each centre i (in vectorised chunks), with i at the origin:
            for each neighbour pair (j, k):
                n = normalise(n_i + n_j + n_k); project everything onto the plane ⟂ n
                circumcircle of the projected (i, j, k): centre c, radius R
                keep if R <= max_circumradius (drops open-rim slivers)
                     and no other neighbour l lies inside or on the circle
                     (neighbours whose normal faces away from n, e.g. the other
                     side of a thin sheet, are ignored)
            orient (i, j, k) counter-clockwise about i's normal
    """
    n = len(points)
    table = _adjacency(pairs, n)
    m = table.shape[1]
    jj, kk = np.triu_indices(m, 1)
    corner = (np.arange(m)[None, None, :] == jj[None, :, None]) | \
             (np.arange(m)[None, None, :] == kk[None, :, None])      # (1, P, m)
    out = []
    for start in range(0, n, _CHUNK):
        rows = np.arange(start, min(start + _CHUNK, n))
        nb = table[rows]                                   # (C, m)
        valid = nb >= 0
        safe_nb = np.where(valid, nb, 0)
        d = points[safe_nb] - points[rows, None, :]        # (C, m, 3)
        nrm_nb = normals[safe_nb]                          # (C, m, 3)
        n_t = normals[rows, None, :] + nrm_nb[:, jj] + nrm_nb[:, kk]   # (C, P, 3)
        n_t /= np.maximum(np.linalg.norm(n_t, axis=-1, keepdims=True), 1e-300)
        # Project every neighbour onto each triangle's plane (C, P, m, 3).
        q = d[:, None, :, :] - np.einsum("cmk,cpk->cpm", d, n_t)[..., None] * n_t[:, :, None, :]
        P = len(jj)
        a = q[:, np.arange(P), jj]                          # (C, P, 3)
        b = q[:, np.arange(P), kk]
        axb = np.cross(a, b)
        axb2 = np.sum(axb**2, axis=-1)
        safe = np.where(axb2 > 1e-24 * max_circumradius**4, axb2, np.inf)
        # Circumcentre of (0, a, b): (|a|² b − |b|² a) × (a × b) / (2 |a × b|²)
        c = np.cross(np.sum(a**2, -1)[..., None] * b - np.sum(b**2, -1)[..., None] * a, axb)
        c /= 2.0 * safe[..., None]
        r2 = np.sum(c**2, axis=-1)                         # (C, P)
        ok = valid[:, jj] & valid[:, kk] & np.isfinite(safe) & (r2 <= max_circumradius**2)
        dist2 = np.sum((q - c[:, :, None, :]) ** 2, axis=-1)          # (C, P, m)
        same_side = np.einsum("cmk,cpk->cpm", nrm_nb, n_t) > 0
        # A point ON the circle counts as inside, so an exact tie (co-circular
        # quadruple) rejects both diagonals: a counted hole, never an overlap.
        inside = ((dist2 <= r2[..., None] * (1 + 1e-9)) & same_side
                  & valid[:, None, :] & ~corner)
        ok &= ~inside.any(axis=-1)
        c_idx, p_idx = np.nonzero(ok)
        i = rows[c_idx]
        j = nb[c_idx, jj[p_idx]]
        k = nb[c_idx, kk[p_idx]]
        ccw = np.einsum("tk,tk->t", axb[c_idx, p_idx], normals[i]) > 0
        out.append(np.column_stack([i, np.where(ccw, j, k), np.where(ccw, k, j)]))
    return np.concatenate(out) if out else np.zeros((0, 3), dtype=np.int64)


# FD. agree_triangles
# Signature: (F, 3) proposed triangles -> TRIANGLES
def agree_triangles(proposals: np.ndarray) -> np.ndarray:
    """Keep the triangles that all three of their corners proposed.

    Algorithm:
        canonical key = the sorted corner indices
        count proposals per key; keep keys proposed exactly 3 times
        keep the orientation of the first proposal
    """
    if len(proposals) == 0:
        return np.zeros((0, 3), dtype=np.int64)
    keys = np.sort(proposals, axis=1)
    _, first, counts = np.unique(keys, axis=0, return_index=True, return_counts=True)
    return proposals[np.sort(first[counts == 3])].astype(np.int64)


# FD. fill_small_holes
# Signature: TRIANGLES, POINTS, float -> TRIANGLES
def fill_small_holes(triangles: np.ndarray, points: np.ndarray, max_circumradius: float) -> np.ndarray:
    """Close holes bounded by 3 or 4 open edges (near-co-circular quadruples, whose
    candidate triangles were each tested in a slightly different plane).

    Algorithm:
        directed boundary edges = edges of triangles with no twin, reversed
            (so a hole's boundary runs counter-clockwise around the hole)
        follow them from vertex to vertex (skip vertices with two outgoing edges)
        cycle of 3: add it as a triangle; cycle of 4: add two triangles split
            along the shorter diagonal
        skip a cycle if a new triangle's circumradius exceeds max_circumradius,
            or if 3 of its corners already form a triangle (the rim of a sheet
            of one or two triangles is not a hole)
    """
    if len(triangles) == 0:
        return triangles
    directed = triangles[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2)
    present = set(map(tuple, directed.tolist()))
    existing = {frozenset(tri) for tri in triangles.tolist()}
    boundary = [(b, a) for a, b in directed.tolist() if (b, a) not in present]
    outgoing = {}
    for a, b in boundary:
        outgoing.setdefault(a, []).append(b)
    added, used = [], set()
    for start in outgoing:
        if start in used or len(outgoing[start]) != 1:
            continue
        cycle = [start]
        while len(cycle) <= 4:
            nxt = outgoing.get(cycle[-1], [])
            if len(nxt) != 1:
                break
            if nxt[0] == start:
                break
            cycle.append(nxt[0])
        if outgoing.get(cycle[-1], [None]) != [start] or len(cycle) not in (3, 4):
            continue
        if len(cycle) == 3:
            new = [cycle]
        else:
            a, b, c, d = cycle
            short_ac = np.linalg.norm(points[a] - points[c]) <= np.linalg.norm(points[b] - points[d])
            new = [[a, b, c], [a, c, d]] if short_ac else [[a, b, d], [b, c, d]]
        corners = [frozenset(cycle[:k] + cycle[k + 1:]) for k in range(len(cycle))]             if len(cycle) == 4 else [frozenset(cycle)]
        if any(c in existing for c in corners):
            continue
        if all(_circumradius(points[tri]) <= max_circumradius for tri in new):
            added += new
            used.update(cycle)
    if not added:
        return triangles
    return np.concatenate([triangles, np.asarray(added, dtype=np.int64)])


def _circumradius(corners: np.ndarray) -> float:
    a, b, c = (np.linalg.norm(corners[i] - corners[j]) for i, j in ((1, 2), (2, 0), (0, 1)))
    area2 = np.linalg.norm(np.cross(corners[1] - corners[0], corners[2] - corners[0]))
    return a * b * c / max(2.0 * area2, 1e-300)


# FD. open_edges
# Signature: TRIANGLES -> (K, 2) int64
def open_edges(triangles: np.ndarray) -> np.ndarray:
    """Undirected edges used by exactly one triangle: the open rim of a sheet,
    or a hole on a closed surface.

    Algorithm:
        list the three edges of every triangle (sorted endpoints)
        count each undirected edge; keep those counted once
    """
    if len(triangles) == 0:
        return np.zeros((0, 2), dtype=np.int64)
    edges = np.sort(triangles[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    return unique[counts == 1].astype(np.int64)


# FD. overused_edges
# Signature: TRIANGLES -> (K, 2) int64
def overused_edges(triangles: np.ndarray) -> np.ndarray:
    """Undirected edges used by more than two triangles (overlapping triangles).
    Zero on smooth surfaces; appears at sharp creases of the input mesh."""
    if len(triangles) == 0:
        return np.zeros((0, 2), dtype=np.int64)
    edges = np.sort(triangles[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    return unique[counts > 2].astype(np.int64)


# FD. valence
# Signature: TRIANGLES, int -> VALENCE
def valence(triangles: np.ndarray, n: int) -> np.ndarray:
    """Number of distinct triangulation neighbours of each centre.

    Algorithm:
        list the three edges of every triangle; deduplicate undirected edges
        count edge endpoints per centre
    """
    if len(triangles) == 0:
        return np.zeros(n, dtype=np.int64)
    edges = np.sort(triangles[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1)
    edges = np.unique(edges, axis=0)
    return np.bincount(edges.ravel(), minlength=n).astype(np.int64)
