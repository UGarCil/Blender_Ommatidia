"""Atomic: interommatidial bristles as tapered spindles standing in the facet corners.

Each bristle is the same small template (rings along its axis, a tip vertex),
placed at a BRISTLE_SITE and tilted from the surface normal toward the comb
direction. Only the template's local coordinates are stored, so length,
thickness and tilt can change without rebuilding the mesh.
"""

from __future__ import annotations

import numpy as np

BRISTLE_SIDES = 6
BRISTLE_SEGMENTS = 6
SINK = 0.5     # x bristle radius; the base starts this far below the corner so no gap shows


# CDD. BRISTLE_TEMPLATE
# template = (axial (K,), angle (K,), faces (Q, 4) int with -1 padding for triangles)
# interp. one bristle in local coordinates: axial 0 = base, 1 = tip; angle around
#         the axis. K = SIDES·SEGMENTS + 1 (the tip). The base ring is open: it
#         sits below the cornea.
def bristle_template(sides: int = BRISTLE_SIDES, segments: int = BRISTLE_SEGMENTS):
    """Local vertex coordinates and faces of one bristle.

    Algorithm:
        rings s = 0..segments−1 at axial s/segments, `sides` vertices each
        tip vertex at axial 1
        quads between consecutive rings; triangles from the last ring to the tip
    """
    axial = np.repeat(np.arange(segments) / segments, sides)
    angle = np.tile(np.arange(sides) * 2 * np.pi / sides, segments)
    axial = np.append(axial, 1.0)
    angle = np.append(angle, 0.0)
    ring = np.arange(segments * sides).reshape(segments, sides)
    nxt = np.roll(np.arange(sides), -1)
    quads = np.stack([ring[:-1], ring[:-1][:, nxt], ring[1:][:, nxt], ring[1:]], axis=-1).reshape(-1, 4)
    tip = segments * sides
    tris = np.stack([ring[-1], ring[-1][nxt], np.full(sides, tip)], axis=-1)
    faces = np.concatenate([quads, np.column_stack([tris, np.full(sides, -1)])])
    return axial, angle, faces


# FD. spindle_radius
# Signature: (K,) axial -> (K,) float in [0, 1]
def spindle_radius(axial: np.ndarray) -> np.ndarray:
    """Relative thickness along a bristle: 0.6 at the base, widest at ~1/4 of the
    length, tapering to a point (SEM panel b)."""
    a = np.asarray(axial)
    return np.minimum(1.0, 0.6 + 1.6 * a) * (1.0 - a) ** 0.8 / 0.75 ** 0.8


# FD. bristle_positions
# Signature: POINTS (sites, per vertex), NORMALS (per vertex), (K,) axial, (K,) angle,
#            length, radius, tilt (degrees), comb direction -> POINTS
def bristle_positions(sites: np.ndarray, normals: np.ndarray, axial: np.ndarray, angle: np.ndarray,
                      length: float, radius: float, tilt_deg: float, comb: np.ndarray) -> np.ndarray:
    """World positions of bristle vertices (one row per vertex).

    Algorithm:
        comb tangent t = comb projected onto the tangent plane (any tangent if parallel)
        axis d = cos(tilt)·n + sin(tilt)·t; (u, v) completes the frame
        p = site − SINK·radius·n + axial·length·d + spindle(axial)·radius·(cos θ u + sin θ v)
    """
    n = normals
    comb = np.asarray(comb, dtype=float)
    t = comb - (n @ comb)[:, None] * n
    weak = np.linalg.norm(t, axis=1) < 1e-6 * np.linalg.norm(comb)
    helper = np.where(np.abs(n[:, :1]) < 0.9, [[1.0, 0, 0]], [[0, 1.0, 0]])
    t[weak] = np.cross(n[weak], helper[weak])
    t /= np.linalg.norm(t, axis=1, keepdims=True)
    a = np.radians(tilt_deg)
    d = np.cos(a) * n + np.sin(a) * t
    u = np.cross(d, t)
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    v = np.cross(d, u)
    r = (spindle_radius(axial) * radius)[:, None]
    return (sites - SINK * radius * n + (axial * length)[:, None] * d
            + r * (np.cos(angle)[:, None] * u + np.sin(angle)[:, None] * v))
