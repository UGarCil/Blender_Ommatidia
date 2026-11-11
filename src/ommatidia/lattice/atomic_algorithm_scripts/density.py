"""Atomic: turn the user's density choice into both SPACING and a count."""

from __future__ import annotations

import math

# The area of one cell of a regular hexagonal lattice with centre spacing a is
# HEX_CELL_FACTOR * a**2.
HEX_CELL_FACTOR = math.sqrt(3.0) / 2.0


# FD. resolve_density
# Signature: float (SURFACE.area), DENSITY_MODE, SPACING, TARGET_COUNT, MAX_COUNT -> SPACING, TARGET_COUNT
def resolve_density(area: float, density_mode: str, spacing: float,
                    target_count: int, max_count: int) -> tuple[float, int]:
    """Derive the missing one of (SPACING, count) from the surface area.

    Algorithm:
        if mode is SPACING: count = max(1, round(area / (sqrt(3)/2 * spacing²)))
        if mode is COUNT:   spacing = sqrt(area / (sqrt(3)/2 * count))
        if count > max_count: raise ValueError (protects Blender from freezing)
        return spacing, count
    """
    if not area > 0:
        raise ValueError(f"Surface area must be > 0, got {area}")
    if density_mode == "SPACING":
        count = max(1, round(area / (HEX_CELL_FACTOR * spacing**2)))
    elif density_mode == "COUNT":
        count = int(target_count)
        spacing = math.sqrt(area / (HEX_CELL_FACTOR * count))
    else:
        raise ValueError(f"Unknown density_mode {density_mode!r}")
    if count > max_count:
        raise ValueError(
            f"Spacing {spacing:.3g} m on a surface of {area:.3g} m² needs {count} ommatidia, "
            f"more than max_count={max_count}. Check model_scale (Blender units per specimen "
            f"metre) or use density_mode='COUNT'."
        )
    return float(spacing), int(count)
