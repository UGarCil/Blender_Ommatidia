"""Compound: the add-on's pipeline for one target object, in the order of main.py.

The Generate button and the live updates call these; main.py spells out the
same steps for batch runs.
"""

from __future__ import annotations

import time

import bpy

from .blender_io.atomic_algorithm_scripts.facets_to_object import (
    facets_object_name,
    reshape_facets_object,
    write_facets_object,
)
from .blender_io.blender_io import lattice_object_name, surface_from_object, write_lattice_object
from .config import Config
from .facets.facets import build_facet_mesh
from .lattice.lattice import build_lattice
from .material.material import assign_materials
from .topology.topology import build_topology


# FD. generate
# Signature: bpy Object (target), CONFIG -> LATTICE_METRICS
def generate(target, config: Config) -> dict:
    """Everything from the mesh to the visible ommatidia.

    Algorithm:
        SURFACE of the target -> LATTICE -> TOPOLOGY -> FACET_MESH
        write the lattice object (data, hidden) and the ommatidia object (visible)
        give the ommatidia object its materials
    """
    surface = surface_from_object(target, config)
    lattice, metrics = build_lattice(surface, config)
    topology, topology_metrics = build_topology(lattice, surface, config)
    metrics.update(topology_metrics)
    t0 = time.perf_counter()
    facet_mesh = build_facet_mesh(lattice, topology, surface, config)
    facets = write_facets_object(target, facet_mesh, config)
    assign_materials(facets)
    metrics["SECONDS"]["facets"] = time.perf_counter() - t0
    metrics["N_FACET_VERTICES"] = len(facet_mesh)
    metrics["DROPPED_INNER_FACES"] = getattr(surface, "dropped_inner_faces", 0)
    write_lattice_object(target, lattice, topology, config, metrics)
    return metrics


# FD. reshape
# Signature: bpy Object (target), CONFIG -> bool
def reshape(target, config: Config) -> bool:
    """Apply the shape settings to the existing ommatidia object (fast). False if
    the target has not been generated yet."""
    facets = bpy.data.objects.get(facets_object_name(target.name))
    if facets is None:
        return False
    reshape_facets_object(facets, config)
    lattice = bpy.data.objects.get(lattice_object_name(target.name))
    if lattice is not None:
        lattice["omm_config"] = config.to_json()     # the .blend keeps saying what made it
    return True
