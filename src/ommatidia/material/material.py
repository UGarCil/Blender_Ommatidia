"""Compound: materials for the ommatidia object. See material.md."""

from __future__ import annotations

from .atomic_algorithm_scripts.sem_look import sem_material


# FD. assign_materials
# Signature: bpy Object -> None
def assign_materials(obj) -> None:
    """Slot 0 = corneal facets, slot 1 = bristles; both start with the SEM look.
    Existing slots are kept, so a user's own materials survive regeneration.

    Algorithm:
        add slots until there are two; fill empty slots with the SEM material
    """
    material = sem_material()
    while len(obj.material_slots) < 2:
        obj.data.materials.append(None)
    for slot in obj.material_slots[:2]:
        if slot.material is None:
            slot.material = material
