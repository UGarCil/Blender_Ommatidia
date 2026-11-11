"""Atomic: how high each facet vertex sits above the surface (dossier §4.3).

Two profiles, blended by FACET_SHAPE:
- dome: a lens over a round base of the facet's inner radius r. INFLATION sets
  its height as a fraction of r (1 = a semi-dome, i.e. a hemisphere);
  POKEYNESS bends the outline from a spherical cap (0) toward a cone (1).
  Beyond r (the facet's corners) it is 0: the small flat valleys where
  bristles sit in Drosophila.
- flat facet: a plateau at the same height that falls off near the walls
  (many Hymenoptera).
CREASE_SHARPNESS sets how narrow the fall-off into the groove is.
"""

from __future__ import annotations

import numpy as np

MIN_CREASE_WIDTH = 0.04   # in rho units; keeps the plateau's walls from going vertical


# FD. crease_width
# Signature: CREASE_SHARPNESS -> float
def crease_width(crease_sharpness: float) -> float:
    """Fall-off width as a fraction of the centre-to-wall distance: 0.5 at
    sharpness 0 (soft), MIN_CREASE_WIDTH at sharpness 1 (sharp)."""
    return max(0.5 * (1.0 - crease_sharpness), MIN_CREASE_WIDTH)


def _smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


# FD. dome_heights
# Signature: (V,) rho_round, (V,) rho_hex, (V,) radius, INFLATION, POKEYNESS, FACET_SHAPE,
#            CREASE_SHARPNESS -> (V,) float (height, specimen metres)
def dome_heights(rho_round: np.ndarray, rho_hex: np.ndarray, radius: np.ndarray, inflation: float,
                 pokeyness: float, facet_shape: float, crease_sharpness: float) -> np.ndarray:
    """Height of every facet vertex.

    Algorithm:
        H = inflation · r
        cap(u)  = spherical cap of height H over base radius r, at distance u·r
        cone(u) = H · (1 − u)
        dome = mix(cap, cone; pokeyness)(rho_round) · smoothstep((1 − rho_round) / w)
        flat = H · smoothstep((1 − rho_hex) / w)
        return mix(dome, flat; facet_shape)
    """
    w = crease_width(crease_sharpness)
    height = inflation * radius
    u = np.clip(rho_round, 0.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        sphere = (radius**2 + height**2) / (2.0 * height)         # radius of the cap's sphere
        cap = np.sqrt(np.maximum(sphere**2 - (u * radius) ** 2, 0.0)) - (sphere - height)
    cap = np.where(height > 0, cap, 0.0)
    cone = height * (1.0 - u)
    dome = ((1.0 - pokeyness) * cap + pokeyness * cone) * _smoothstep((1.0 - rho_round) / w)
    flat = height * _smoothstep((1.0 - rho_hex) / w)
    return (1.0 - facet_shape) * dome + facet_shape * flat


# FD. rim_fade
# Signature: DISTANCES (boundary distance / SPACING), RIM_WIDTH -> (V,) float in [0, 1]
def rim_fade(rim_distance: np.ndarray, rim_width: float) -> np.ndarray:
    """Facets fade to flat over RIM_WIDTH spacings at an open edge (1 elsewhere).
    rim_distance is measured from the facet's centre, which is already at least
    half a spacing from the edge, so the fade starts there."""
    if rim_width <= 0:
        return np.ones_like(rim_distance)
    return _smoothstep((np.asarray(rim_distance) - 0.5) / rim_width)
