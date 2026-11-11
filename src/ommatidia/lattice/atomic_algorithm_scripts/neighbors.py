"""Atomic: fixed-radius neighbour queries with a numpy spatial-hash grid.

Replaces a KD-tree: SciPy is not available in Blender and mathutils is not
allowed in the pure modules. Points are binned into cubes of side `radius`;
any pair closer than `radius` lies in the same or in adjacent cubes (27 cells).
"""

from __future__ import annotations

import numpy as np

# DD. NEIGHBOR_PAIRS
# neighbor_pairs = np.ndarray, shape (M, 2), int64
# interp. every unordered pair (i, j) of POINTS closer than a cut-off radius
# Invariant: i < j in every row; no duplicate rows; rows sorted lexicographically
neighbor_pairs0 = np.zeros((0, 2), dtype=np.int64)  # stub
neighbor_pairs1 = np.array([[0, 1], [0, 2], [1, 2]])

_OFFSETS = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])


def _candidate_pairs(a: np.ndarray, b: np.ndarray, cell: float) -> tuple[np.ndarray, np.ndarray]:
    """All (ia, ib) with b[ib] in a grid cell adjacent to a[ia]'s cell (a superset
    of the pairs closer than `cell`)."""
    if len(a) == 0 or len(b) == 0:
        empty = np.zeros(0, dtype=np.int64)
        return empty, empty
    origin = np.minimum(a.min(axis=0), b.min(axis=0))
    ia_cell = np.floor((a - origin) / cell).astype(np.int64) + 1
    ib_cell = np.floor((b - origin) / cell).astype(np.int64) + 1
    dims = np.maximum(ia_cell.max(axis=0), ib_cell.max(axis=0)) + 2

    def key(c):
        return (c[..., 0] * dims[1] + c[..., 1]) * dims[2] + c[..., 2]

    b_keys = key(ib_cell)
    order = np.argsort(b_keys, kind="stable")
    sorted_keys = b_keys[order]
    out_a, out_b = [], []
    for off in _OFFSETS:
        k = key(ia_cell + off)
        start = np.searchsorted(sorted_keys, k, side="left")
        stop = np.searchsorted(sorted_keys, k, side="right")
        counts = stop - start
        total = int(counts.sum())
        if total == 0:
            continue
        rep_a = np.repeat(np.arange(len(a)), counts)
        within = np.arange(total) - np.repeat(np.cumsum(counts) - counts, counts)
        out_a.append(rep_a)
        out_b.append(order[np.repeat(start, counts) + within])
    if not out_a:
        empty = np.zeros(0, dtype=np.int64)
        return empty, empty
    return np.concatenate(out_a), np.concatenate(out_b)


# FD. neighbor_pairs
# Signature: POINTS, float -> NEIGHBOR_PAIRS
def neighbor_pairs(points: np.ndarray, radius: float) -> np.ndarray:
    """Every pair of points closer than radius.

    Algorithm:
        bin points into cubes of side radius
        for each of the 27 neighbouring cube offsets: list candidate pairs
        keep pairs with i < j and distance < radius
        sort rows lexicographically
    """
    points = np.asarray(points, dtype=float)
    i, j = _candidate_pairs(points, points, radius)
    keep = i < j
    i, j = i[keep], j[keep]
    d2 = np.sum((points[i] - points[j]) ** 2, axis=1)
    keep = d2 < radius**2
    i, j = i[keep], j[keep]
    order = np.lexsort((j, i))
    return np.column_stack([i[order], j[order]]).astype(np.int64)


# FD. assign_nearest
# Signature: POINTS (samples), POINTS (centres), float -> (K,) int64
def assign_nearest(samples: np.ndarray, centers: np.ndarray, radius: float) -> np.ndarray:
    """Index of the nearest centre for each sample (used for coverage checks).

    Algorithm:
        candidate pairs within the grid neighbourhood of side radius
        keep only candidates closer than radius (a farther candidate may lose to
            a centre outside the searched cells)
        for each sample keep the candidate with the smallest distance
        samples with no candidate (a hole wider than radius): brute force, in chunks
    """
    samples = np.asarray(samples, dtype=float)
    centers = np.asarray(centers, dtype=float)
    if len(centers) == 0:
        raise ValueError("assign_nearest needs at least one centre")
    s, c = _candidate_pairs(samples, centers, radius)
    d2 = np.sum((samples[s] - centers[c]) ** 2, axis=1)
    close = d2 < radius**2
    s, c, d2 = s[close], c[close], d2[close]
    # Sort by (sample, distance) and take the first row of each sample.
    order = np.lexsort((d2, s))
    s, c = s[order], c[order]
    first = np.ones(len(s), dtype=bool)
    first[1:] = s[1:] != s[:-1]
    nearest = np.full(len(samples), -1, dtype=np.int64)
    nearest[s[first]] = c[first]
    missing = np.flatnonzero(nearest < 0)
    for chunk in np.array_split(missing, max(1, len(missing) // 2048 + 1)):
        if len(chunk):
            dist = np.linalg.norm(samples[chunk, None, :] - centers[None, :, :], axis=2)
            nearest[chunk] = dist.argmin(axis=1)
    return nearest


# FD. nearest_distances
# Signature: POINTS (samples), POINTS (centres), float -> DISTANCES
def nearest_distances(samples: np.ndarray, centers: np.ndarray, radius: float) -> np.ndarray:
    """Distance from each sample to its nearest centre (inf if there are no centres)."""
    if len(centers) == 0:
        return np.full(len(samples), np.inf)
    nearest = assign_nearest(samples, centers, radius)
    return np.linalg.norm(np.asarray(samples, dtype=float) - centers[nearest], axis=1)
