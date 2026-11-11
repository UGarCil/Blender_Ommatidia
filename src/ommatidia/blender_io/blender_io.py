"""Compound: the only bridge between bpy data and the pure modules. See blender_io.md."""

from __future__ import annotations

import bpy
import numpy as np

from ..config import Config
from ..lattice.atomic_algorithm_scripts.density import HEX_CELL_FACTOR
from .atomic_algorithm_scripts.lattice_to_object import lattice_object_name, write_lattice_object
from .atomic_algorithm_scripts.mesh_to_surface import BvhSurface, inner_layer_faces, mesh_arrays

# Name of the Geometry Nodes modifier this add-on puts on the target (from M3).
# It is skipped when reading the target, so domes never feed back into the surface.
MODIFIER_NAME = "Ommatidia"
# Two opposite-facing skins closer than this (x estimated spacing) are one thin solid.
THIN_GAP = 1.0


# FD. surface_from_object
# Signature: bpy Object, CONFIG -> SURFACE
def surface_from_object(obj, config: Config) -> BvhSurface:
    """Build the SURFACE of a mesh object, in specimen metres.

    Algorithm:
        evaluated mesh (without our own modifier), world space / MODEL_SCALE
        estimate the facet spacing from the whole area
        thin solid (two skins closer than THIN_GAP spacings)? drop the inner skin
        -> BVH-backed SURFACE (remembers how many faces were dropped)
    """
    depsgraph = bpy.context.evaluated_depsgraph_get()
    vertices, triangles, normals = mesh_arrays(obj, depsgraph, config.model_scale, MODIFIER_NAME)
    a, b, c = (vertices[triangles[:, k]] for k in range(3))
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()
    spacing = config.spacing if config.density_mode == "SPACING" else         float(np.sqrt(area / (HEX_CELL_FACTOR * config.target_count)))
    inner = inner_layer_faces(vertices, triangles, THIN_GAP * spacing)
    surface = BvhSurface(vertices, triangles[~inner], normals[~inner])
    surface.dropped_inner_faces = int(inner.sum())
    return surface


__all__ = [
    "MODIFIER_NAME",
    "surface_from_object",
    "write_lattice_object",
    "lattice_object_name",
]
