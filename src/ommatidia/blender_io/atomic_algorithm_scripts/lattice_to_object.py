"""Atomic: write LATTICE + TOPOLOGY into one Blender mesh object (LATTICE_OBJECT)."""

from __future__ import annotations

import json

import bpy
import numpy as np

LATTICE_SUFFIX = "_OmmatidiaLattice"
KIND_CENTER, KIND_CELL_VERTEX = 0, 1

# CDD. LATTICE_OBJECT
# A hidden mesh object "<target>_OmmatidiaLattice", parented to the target, world positions:
#   vertices = LATTICE.centers ∪ TOPOLOGY.cell_points; edges = TOPOLOGY.cell_edges
#   point attributes: omm_kind INT (0 centre, 1 cell vertex), omm_normal FLOAT_VECTOR,
#                     omm_bristle BOOLEAN, omm_valence INT (centres; 0 for cell vertices)
#   custom properties: omm_config (CONFIG JSON), omm_metrics (LATTICE_METRICS JSON),
#                      omm_target (name of the covered object)


# FD. lattice_object_name
# Signature: str -> str
def lattice_object_name(target_name: str) -> str:
    """Name of the lattice object that belongs to a target object."""
    return target_name + LATTICE_SUFFIX


# FD. write_lattice_object
# Signature: bpy Object, LATTICE, TOPOLOGY, CONFIG, LATTICE_METRICS -> LATTICE_OBJECT
def write_lattice_object(target, lattice, topology, config, metrics: dict):
    """Create or replace the lattice object of target.

    Algorithm:
        positions = [centres; cell points] × MODEL_SCALE (back to model units)
        edges = cell edges shifted past the centre block
        build a fresh mesh; write the point attributes
        reuse the existing object of that name (swap its mesh) or create one
            in the target's first collection; parent it to the target; hide it
        store CONFIG, metrics and the target's name as custom properties
    """
    n_c = len(lattice.centers)
    positions = np.concatenate([lattice.centers, topology.cell_points]) * config.model_scale
    edges = (topology.cell_edges + n_c).tolist()
    name = lattice_object_name(target.name)

    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(positions.tolist(), edges, [])
    n = len(positions)
    kind = np.concatenate([np.full(n_c, KIND_CENTER), np.full(len(topology.cell_points), KIND_CELL_VERTEX)])
    normals = np.concatenate([lattice.normals, topology.cell_normals]).reshape(-1, 3)
    bristle = np.concatenate([np.zeros(n_c, dtype=bool), topology.bristle_mask])
    valence = np.concatenate([topology.valence, np.zeros(len(topology.cell_points), dtype=np.int64)])
    _write_attribute(mesh, "omm_kind", "INT", "value", kind.astype(np.int32), n)
    _write_attribute(mesh, "omm_normal", "FLOAT_VECTOR", "vector", normals.astype(np.float32).ravel(), n)
    _write_attribute(mesh, "omm_bristle", "BOOLEAN", "value", bristle, n)
    _write_attribute(mesh, "omm_valence", "INT", "value", valence.astype(np.int32), n)
    mesh.update()

    obj = bpy.data.objects.get(name)
    if obj is None:
        obj = bpy.data.objects.new(name, mesh)
        collection = target.users_collection[0] if target.users_collection else bpy.context.scene.collection
        collection.objects.link(obj)
    else:
        old = obj.data
        obj.data = mesh
        if old.users == 0:
            bpy.data.meshes.remove(old)
    obj.parent = target
    obj.matrix_parent_inverse = target.matrix_world.inverted()
    obj.matrix_basis.identity()
    obj.display_type = "WIRE"
    obj.hide_render = True
    try:
        obj.hide_set(True)      # data for the add-on; the visible result is <target>_Ommatidia
    except RuntimeError:        # not in the active view layer
        pass
    obj["omm_config"] = config.to_json()
    obj["omm_metrics"] = json.dumps(metrics, indent=2)
    obj["omm_target"] = target.name
    return obj


# FD. read_lattice_metrics
# Signature: bpy Object (target) -> LATTICE_METRICS | None
def read_lattice_metrics(target) -> dict | None:
    """Metrics stored on the target's lattice object, or None if it has none yet."""
    obj = bpy.data.objects.get(lattice_object_name(target.name))
    if obj is None or "omm_metrics" not in obj:
        return None
    return json.loads(obj["omm_metrics"])


def _write_attribute(mesh, name: str, data_type: str, field: str, values, n: int):
    attribute = mesh.attributes.new(name, data_type, "POINT")
    if n:
        attribute.data.foreach_set(field, values)
