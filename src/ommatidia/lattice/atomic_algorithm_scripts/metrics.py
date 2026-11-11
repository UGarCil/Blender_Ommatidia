"""Atomic: regularity measures of a set of centres."""

from __future__ import annotations

import numpy as np


# FD. nearest_neighbor_distances
# Signature: POINTS, NEIGHBOR_PAIRS -> DISTANCES
def nearest_neighbor_distances(points: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    """Distance from each point to its closest neighbour among the pairs (inf if none).

    Algorithm:
        pair lengths; for each point take the minimum over the pairs it is in
    """
    nn = np.full(len(points), np.inf)
    if len(pairs):
        d = np.linalg.norm(points[pairs[:, 0]] - points[pairs[:, 1]], axis=1)
        np.minimum.at(nn, pairs[:, 0], d)
        np.minimum.at(nn, pairs[:, 1], d)
    return nn


# FD. spacing_cv
# Signature: DISTANCES, (N,) bool -> float
def spacing_cv(nn_distances: np.ndarray, mask: np.ndarray) -> float:
    """Coefficient of variation (std / mean) of the finite, masked distances.

    Returns nan when there are fewer than two such distances.
    """
    d = nn_distances[mask & np.isfinite(nn_distances)]
    if len(d) < 2:
        return float("nan")
    return float(d.std() / d.mean())
