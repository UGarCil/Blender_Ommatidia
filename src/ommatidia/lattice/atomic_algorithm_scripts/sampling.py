"""Atomic: random points on a SURFACE, and dart throwing to fill the seams between patches."""

from __future__ import annotations

import numpy as np

from .neighbors import _candidate_pairs, neighbor_pairs


# FD. sample_candidates
# Signature: SURFACE, int, RNG -> POINTS
def sample_candidates(surface, n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw n points uniformly per unit area.

    Algorithm:
        delegate to surface.sample_area_weighted (the SURFACE owns its geometry)
    """
    return np.asarray(surface.sample_area_weighted(int(n), rng), dtype=float).reshape(-1, 3)


# FD. poisson_filter
# Signature: POINTS, float -> (K,) bool
def poisson_filter(candidates: np.ndarray, radius: float) -> np.ndarray:
    """Greedy dart throwing: keep candidates in order unless a kept one is closer than radius.

    Algorithm:
        list all candidate pairs closer than radius (spatial hash)
        build an adjacency list
        for each candidate in order: if still alive, kill all its close neighbours
        return the alive mask
    """
    n = len(candidates)
    pairs = neighbor_pairs(candidates, radius)
    both = np.concatenate([pairs, pairs[:, ::-1]])
    both = both[np.argsort(both[:, 0], kind="stable")]
    starts = np.searchsorted(both[:, 0], np.arange(n + 1))
    others = both[:, 1]
    alive = np.ones(n, dtype=bool)
    for i in range(n):
        if alive[i] and starts[i] != starts[i + 1]:
            alive[others[starts[i]:starts[i + 1]]] = False
    return alive


# FD. fill_mask
# Signature: POINTS (candidates), POINTS (existing), NORMALS, NORMALS, float -> (K,) bool
def fill_mask(candidates: np.ndarray, existing: np.ndarray, candidate_normals: np.ndarray,
              existing_normals: np.ndarray, radius: float) -> np.ndarray:
    """Dart throwing into the gaps of an existing set: keep candidates in order unless
    an existing point or an earlier kept candidate (facing the same side) is closer
    than radius. Existing points are never removed.

    Algorithm:
        drop candidates within radius of an existing point facing the same side
        Poisson-disk filter the rest among themselves (poisson_filter)
    """
    keep = np.ones(len(candidates), dtype=bool)
    i, j = _candidate_pairs(candidates, existing, radius)
    close = (np.sum((candidates[i] - existing[j]) ** 2, axis=1) < radius**2) & (
        np.sum(candidate_normals[i] * existing_normals[j], axis=1) > 0)
    keep[i[close]] = False
    rest = np.flatnonzero(keep)
    keep[rest] = poisson_filter(candidates[rest], radius)
    return keep
