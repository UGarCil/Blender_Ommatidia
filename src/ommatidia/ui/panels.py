"""Sidebar panels (View3D > N > Ommatidia). Corneal-nipple controls arrive with M6."""

from __future__ import annotations

import bpy

from ..blender_io.atomic_algorithm_scripts.lattice_to_object import read_lattice_metrics
from .operators import is_coverable


class _Sub:
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Ommatidia"
    bl_parent_id = "OMMATIDIA_PT_main"


class OMMATIDIA_PT_main(bpy.types.Panel):
    bl_label = "Ommatidia"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Ommatidia"

    @classmethod
    def poll(cls, context):
        return is_coverable(context.active_object)

    def draw(self, context):
        props = context.active_object.ommatidia
        layout = self.layout
        layout.use_property_split = True
        col = layout.column(heading="Density")
        col.prop(props, "density_mode")
        col.prop(props, "spacing" if props.density_mode == "SPACING" else "target_count")
        col.prop(props, "model_scale")
        row = layout.row(align=True)
        row.scale_y = 1.4
        row.operator("ommatidia.generate", icon="MESH_ICOSPHERE")
        layout.prop(props, "live_update")
        if props.status:
            layout.label(text=props.status, icon="ERROR")


class OMMATIDIA_PT_lattice(_Sub, bpy.types.Panel):
    bl_label = "Lattice"

    def draw(self, context):
        props = context.active_object.ommatidia
        layout = self.layout
        layout.use_property_split = True
        col = layout.column()
        col.prop(props, "origin_source")
        if props.origin_source == "OBJECT":
            col.prop(props, "origin_object")
        elif props.origin_source == "POINT":
            col.prop(props, "origin_point")
        col.operator("ommatidia.use_cursor", icon="PIVOT_CURSOR")
        col.prop(props, "lattice_rotation_deg")
        col.separator()
        col.prop(props, "seed")
        col.prop(props, "relax_iterations")
        col.prop(props, "irregularity")


class OMMATIDIA_PT_facets(_Sub, bpy.types.Panel):
    bl_label = "Facets"

    def draw(self, context):
        props = context.active_object.ommatidia
        layout = self.layout
        layout.use_property_split = True
        col = layout.column()
        col.prop(props, "inflation", slider=True)
        col.prop(props, "facet_shape", slider=True)
        col.prop(props, "pokeyness", slider=True)
        col.prop(props, "crease_sharpness", slider=True)
        col.prop(props, "rim_width")
        col.prop(props, "subdivision_level")


class OMMATIDIA_PT_bristles(_Sub, bpy.types.Panel):
    bl_label = "Bristles"

    def draw(self, context):
        props = context.active_object.ommatidia
        layout = self.layout
        layout.use_property_split = True
        col = layout.column()
        col.prop(props, "bristle_pattern")
        row = col.row()
        row.enabled = props.bristle_pattern == "RANDOM"
        row.prop(props, "bristle_frequency", slider=True)
        sub = col.column()
        sub.enabled = props.bristle_pattern != "NONE"
        sub.prop(props, "bristle_length")
        sub.prop(props, "bristle_radius")
        sub.prop(props, "bristle_tilt_deg")
        sub.prop(props, "comb_direction")


class OMMATIDIA_PT_diagnostics(_Sub, bpy.types.Panel):
    bl_label = "Diagnostics"
    bl_options = {"DEFAULT_CLOSED"}

    def draw(self, context):
        metrics = read_lattice_metrics(context.active_object)
        col = self.layout.column(align=True)
        if metrics is None:
            col.label(text="Not generated yet")
            return
        col.label(text=f"Ommatidia: {metrics['N']}  (target {metrics['N_TARGET']})")
        col.label(text=f"Spacing: {metrics['SPACING'] * 1e6:.2f} µm   CV {metrics['SPACING_CV']:.3f}")
        col.label(text=f"Hexagonal: {metrics['HEX_FRACTION']:.1%}")
        col.label(text=f"Pentagons / heptagons: {metrics['N_PENTAGONS']} / {metrics['N_HEPTAGONS']}")
        col.label(text=f"Regular patches: {metrics.get('N_PATCHES', 0)}   seam fill: {metrics.get('N_SEAM', 0)}")
        col.label(text=f"Bristles: {metrics['N_BRISTLES']}")
        if metrics.get("DROPPED_INNER_FACES", 0):
            col.label(text="Thin solid: covered its outer skin only", icon="INFO")
        broken = metrics.get("HOLES_ON_CLOSED_SURFACE", 0) + metrics.get("N_OVERLAP_EDGES", 0)
        if broken:
            col.label(text=f"{broken} broken facet edges: sharp creases?", icon="ERROR")
            col.label(text="Add a Subdivision Surface or Bevel modifier")
        seconds = metrics.get("SECONDS", {})
        col.label(text=f"Time: {sum(seconds.values()):.1f} s")


CLASSES = (OMMATIDIA_PT_main, OMMATIDIA_PT_lattice, OMMATIDIA_PT_facets,
           OMMATIDIA_PT_bristles, OMMATIDIA_PT_diagnostics)
