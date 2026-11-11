"""Atomic: one step of each relaxation method, plus boundary trimming and jitter.

None of these project back onto the SURFACE; the compound does that after
each step (one action per function).
"""

from __future__ import annotations

import numpy as np

# Tuned on a 750-point sphere (agent/logbook.md, 2026-11-03). The radius reaches
# the second ring of a hexagonal lattice (at sqrt(3)·spacing); shorter ranges
# (< ~1.7) lose shear stiffness and stay disordered. Momentum anneals the
# pentagon-heptagon defects about 3x faster than plain gradient steps.
REPULSION_RADIUS = 2.0      # x spacing
REPULSION_STEP = 0.6        # x spacing per unit force
REPULSION_MOMENTUM = 0.9    # fraction of the previous move carried over
REPULSION_MAX_MOVE = 0.15   # x spacing, clamps the move of a single step

# DD. VELOCITIES
# velocities = np.ndarray, shape (N, 3), float64
# interp. the previous relaxation move of each centre (specimen metres), carried
#         over as momentum; tangential to the surface
velocities0 = np.zeros((0, 3))  # stub
velocities1 = np.zeros((1, 3))  # at rest


def _tangential(vectors: np.ndarray, normals: np.ndarray) -> np.ndarray:
    return vectors - np.sum(vectors * normals, axis=1, keepdims=True) * normals


# FD. repulsion_step
# Signature: POINTS, VELOCITIES, NORMALS, NEIGHBOR_PAIRS, SPACING -> POINTS, VELOCITIES
def repulsion_step(points: np.ndarray, velocities: np.ndarray, normals: np.ndarray,
                   pairs: np.ndarray, spacing: float) -> tuple[np.ndarray, np.ndarray]:
    """Push every centre away from its neighbours along the surface (Turk 1991
    style repulsion, with heavy-ball momentum).

    Algorithm:
        r = REPULSION_RADIUS · spacing
        for each pair (i, j) at distance d < r: weight w = (1 − d/r)²
            force on i += w · unit(p_i − p_j); force on j −= the same
        move = MOMENTUM · (tangential part of previous move)
               + STEP · spacing · (tangential part of force)
        clamp each move to MAX_MOVE · spacing
        return points + move, move
    """
    force = np.zeros_like(points)
    if len(pairs):
        i, j = pairs[:, 0], pairs[:, 1]
        diff = points[i] - points[j]
        dist = np.linalg.norm(diff, axis=1)
        r = REPULSION_RADIUS * spacing
        w = np.where(dist < r, (1.0 - dist / r) ** 2, 0.0) / np.maximum(dist, 1e-12 * spacing)
        f = diff * w[:, None]
        for axis in range(3):
            force[:, axis] = (np.bincount(i, f[:, axis], len(points))
                              - np.bincount(j, f[:, axis], len(points)))
    move = (REPULSION_MOMENTUM * _tangential(velocities, normals)
            + REPULSION_STEP * spacing * _tangential(force, normals))
    length = np.linalg.norm(move, axis=1, keepdims=True)
    move *= np.minimum(1.0, REPULSION_MAX_MOVE * spacing / np.maximum(length, 1e-300))
    return points + move, move


# FD. trim_boundary
# Signature: DISTANCES, SPACING -> (N,) bool
def trim_boundary(boundary_distance: np.ndarray, spacing: float) -> np.ndarray:
    """Mask of centres far enough from an open edge to keep (>= 0.5 · spacing).

    Relaxation piles centres onto an open rim; those half-cells are dropped and
    the GN rim falloff fades the edge instead.
    """
    return np.asarray(boundary_distance) >= 0.5 * spacing


# FD. jitter_tangent
# Signature: POINTS, NORMALS, IRREGULARITY, SPACING, RNG -> POINTS
def jitter_tangent(points: np.ndarray, normals: np.ndarray, irregularity: float,
                   spacing: float, rng: np.random.Generator) -> np.ndarray:
    """Add Gaussian tangential noise with standard deviation irregularity · spacing.

    Algorithm:
        draw a 3D Gaussian offset per point, keep its tangential part, scale, add
    """
    offset = _tangential(rng.normal(size=points.shape), normals)
    return points + irregularity * spacing * offset / np.sqrt(2.0)  # 2 tangential dims
