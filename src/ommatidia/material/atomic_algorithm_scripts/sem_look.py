"""Atomic: a material that reads like a scanning electron micrograph.

SEM brightness comes mostly from the edge effect: surfaces seen at a grazing
angle emit more secondary electrons, so outlines and slopes glow and faces
seen head-on stay dark. Layer Weight 'Facing' gives exactly that ramp.
Sockets are linked by index: their names changed in Blender 5 ("Fac" -> "Factor").
"""

from __future__ import annotations

import bpy

MATERIAL_NAME = "Ommatidia SEM"
EDGE_BLEND = 0.35          # Layer Weight blend: how fast slopes brighten
DARK, BRIGHT = 0.02, 0.85  # grey levels of faces seen head-on and edge-on
DIFFUSE_SHARE = 0.25       # a little ordinary shading so shapes read in renders


# FD. sem_material
# Signature: None -> bpy Material
def sem_material():
    """Get the add-on's SEM material, building it the first time.

    Algorithm:
        reuse the material of that name if it exists (users may edit it)
        else: Layer Weight (Facing) -> Color Ramp (dark -> bright) -> Emission,
              mixed with a light grey Diffuse -> Material Output
    """
    material = bpy.data.materials.get(MATERIAL_NAME)
    if material is not None:
        return material
    material = bpy.data.materials.new(MATERIAL_NAME)
    if bpy.app.version < (5, 0, 0):
        material.use_nodes = True
    material.diffuse_color = (0.55, 0.55, 0.55, 1.0)      # solid viewport colour
    nodes, links = material.node_tree.nodes, material.node_tree.links
    nodes.clear()
    weight = nodes.new("ShaderNodeLayerWeight")
    weight.inputs[0].default_value = EDGE_BLEND
    ramp = nodes.new("ShaderNodeValToRGB")
    ramp.color_ramp.elements[0].color = (DARK, DARK, DARK, 1.0)
    ramp.color_ramp.elements[1].color = (BRIGHT, BRIGHT, BRIGHT, 1.0)
    emission = nodes.new("ShaderNodeEmission")
    diffuse = nodes.new("ShaderNodeBsdfDiffuse")
    diffuse.inputs[0].default_value = (0.5, 0.5, 0.5, 1.0)
    mix = nodes.new("ShaderNodeMixShader")
    mix.inputs[0].default_value = DIFFUSE_SHARE
    output = nodes.new("ShaderNodeOutputMaterial")
    links.new(weight.outputs[1], ramp.inputs[0])          # Facing -> ramp
    links.new(ramp.outputs[0], emission.inputs[0])
    links.new(emission.outputs[0], mix.inputs[1])
    links.new(diffuse.outputs[0], mix.inputs[2])
    links.new(mix.outputs[0], output.inputs[0])
    for x, node in enumerate((weight, ramp, emission, mix, output)):
        node.location = (x * 250, 0)
    diffuse.location = (500, -200)
    return material
