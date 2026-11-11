"""Operators: thin callers of the pipeline. No algorithms here."""

from __future__ import annotations

import bpy

from .live import regenerate


def is_coverable(obj) -> bool:
    """A mesh object that is not itself made by this add-on."""
    return obj is not None and obj.type == "MESH" and "omm_target" not in obj


class OMMATIDIA_OT_generate(bpy.types.Operator):
    """Cover the active mesh with ommatidia (regular lattice, domed facets and bristles)"""

    bl_idname = "ommatidia.generate"
    bl_label = "Generate Ommatidia"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return is_coverable(context.active_object)

    def execute(self, context):
        target = context.active_object
        metrics = regenerate(target)
        if metrics is None:
            self.report({"ERROR"}, target.ommatidia.status)
            return {"CANCELLED"}
        self.report({"INFO"}, f"{metrics['N']} ommatidia, {metrics['HEX_FRACTION']:.0%} hexagonal, "
                              f"{metrics['N_BRISTLES']} bristles")
        return {"FINISHED"}


class OMMATIDIA_OT_cursor_to_centre(bpy.types.Operator):
    """Use the 3D cursor as the centre of the regular patch (place it with Shift + right-click)"""

    bl_idname = "ommatidia.use_cursor"
    bl_label = "Use 3D Cursor"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return is_coverable(context.active_object)

    def execute(self, context):
        props = context.active_object.ommatidia
        props.origin_point = context.scene.cursor.location      # live-regenerates if POINT
        props.origin_source = "POINT"
        return {"FINISHED"}


CLASSES = (OMMATIDIA_OT_generate, OMMATIDIA_OT_cursor_to_centre)
