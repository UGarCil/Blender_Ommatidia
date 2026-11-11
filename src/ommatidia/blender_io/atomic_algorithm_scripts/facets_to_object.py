"""Atomic: write FACET_MESH as the visible ommatidia object, and reshape it in place."""

from __future__ import annotations

import bpy
import numpy as np

from ...facets.facets import FacetMesh, shape_positions

FACETS_SUFFIX = "_Ommatidia"

# FACET_MESH field -> (attribute name, data type, foreach field, components)
_ATTRIBUTES = {
    "base": ("omm_base", "FLOAT_VECTOR", "vector", 3),
    "normal": ("omm_normal", "FLOAT_VECTOR", "vector", 3),
    "kind": ("omm_kind", "INT", "value", 1),
    "rho_round": ("omm_rho_round", "FLOAT", "value", 1),
    "rho_hex": ("omm_rho_hex", "FLOAT", "value", 1),
    "radius": ("omm_radius", "FLOAT", "value", 1),
    "rim": ("omm_rim", "FLOAT", "value", 1),
    "axial": ("omm_axial", "FLOAT", "value", 1),
    "angle": ("omm_angle", "FLOAT", "value", 1),
}

# CDD. FACETS_OBJECT
# A mesh object "<target>_Ommatidia", parented to the target (world positions kept):
#   vertices = shape_positions(FACET_MESH, CONFIG) × MODEL_SCALE; faces from FACET_MESH
#   material slots: 0 = facets, 1 = bristles (face material_index = face_kind)
#   point attributes = the FACET_MESH per-vertex fields (specimen metres), so the
#   shape can be recomputed without the surface
#   custom properties: omm_spacing (specimen metres), omm_target (target name)


# FD. facets_object_name
# Signature: str -> str
def facets_object_name(target_name: str) -> str:
    """Name of the visible ommatidia object that belongs to a target object."""
    return target_name + FACETS_SUFFIX


# FD. write_facets_object
# Signature: bpy Object (target), FACET_MESH, CONFIG -> FACETS_OBJECT
def write_facets_object(target, mesh: FacetMesh, config):
    """Create or replace the ommatidia object of target.

    Algorithm:
        build a fresh mesh: vertices at the shaped positions × MODEL_SCALE,
            loops and polygons straight from the arrays, smooth shading,
            material index per face
        store the FACET_MESH fields as point attributes
        reuse the object of that name (swap its mesh, keep its materials) or create
            one in the target's first collection; parent it to the target
    """
    name = facets_object_name(target.name)
    data = bpy.data.meshes.new(name)
    positions = shape_positions(mesh, config) * config.model_scale
    data.vertices.add(len(mesh))
    data.vertices.foreach_set("co", positions.astype(np.float32).ravel())
    data.loops.add(len(mesh.loop_vertices))
    data.loops.foreach_set("vertex_index", mesh.loop_vertices.astype(np.int32))
    data.polygons.add(len(mesh.loop_start))
    data.polygons.foreach_set("loop_start", mesh.loop_start.astype(np.int32))
    data.update(calc_edges=True)
    data.polygons.foreach_set("material_index", mesh.face_kind.astype(np.int32))
    data.polygons.foreach_set("use_smooth", np.ones(len(mesh.loop_start), dtype=bool))
    for field, (attr, data_type, slot, _) in _ATTRIBUTES.items():
        values = getattr(mesh, field)
        values = values.astype(np.int32) if data_type == "INT" else values.astype(np.float32)
        attribute = data.attributes.new(attr, data_type, "POINT")
        if len(mesh):
            attribute.data.foreach_set(slot, values.ravel())
    data.update()

    obj = bpy.data.objects.get(name)
    if obj is None:
        obj = bpy.data.objects.new(name, data)
        collection = target.users_collection[0] if target.users_collection else bpy.context.scene.collection
        collection.objects.link(obj)
    else:
        old = obj.data
        for material in old.materials:
            data.materials.append(material)
        obj.data = data
        if old.users == 0:
            bpy.data.meshes.remove(old)
    obj.parent = target
    obj.matrix_parent_inverse = target.matrix_world.inverted()
    obj.matrix_basis.identity()
    obj["omm_spacing"] = mesh.spacing
    obj["omm_target"] = target.name
    return obj


# FD. read_facet_mesh
# Signature: FACETS_OBJECT -> FACET_MESH (faces left empty: reshaping does not need them)
def read_facet_mesh(obj) -> FacetMesh:
    """Read the per-vertex FACET_MESH fields back from the object's attributes."""
    data = obj.data
    n = len(data.vertices)
    values = {}
    for field, (attr, data_type, slot, size) in _ATTRIBUTES.items():
        buffer = np.empty(n * size, dtype=np.int32 if data_type == "INT" else np.float32)
        data.attributes[attr].data.foreach_get(slot, buffer)
        buffer = buffer.astype(np.int64 if data_type == "INT" else np.float64)
        values[field] = buffer.reshape(n, 3) if size == 3 else buffer
    empty = np.zeros(0, dtype=np.int64)
    return FacetMesh(**values, loop_start=empty, loop_vertices=empty, face_kind=empty,
                     spacing=float(obj["omm_spacing"]))


# FD. reshape_facets_object
# Signature: FACETS_OBJECT, CONFIG -> None
def reshape_facets_object(obj, config) -> None:
    """Move the vertices to the shape of CONFIG (domes and bristles), in place.

    Algorithm:
        FACET_MESH fields from the attributes -> shape_positions × MODEL_SCALE -> vertex coordinates
    """
    mesh = read_facet_mesh(obj)
    positions = shape_positions(mesh, config) * config.model_scale
    obj.data.vertices.foreach_set("co", positions.astype(np.float32).ravel())
    obj.data.update()
