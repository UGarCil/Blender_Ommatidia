"""Ommatidia: cover any mesh with a regular lattice of insect ommatidia.

This file only registers the Blender parts. Blender modules are imported inside
register()/unregister() so the pure packages (lattice/, topology/) and config.py
can be imported and tested without bpy.
"""

_registered = []


def register():
    import bpy

    from .blender_io.atomic_algorithm_scripts.config_props import OmmatidiaSettings
    from .ui import operators, panels

    for cls in (OmmatidiaSettings, *operators.CLASSES, *panels.CLASSES):
        bpy.utils.register_class(cls)
        _registered.append(cls)
    bpy.types.Object.ommatidia = bpy.props.PointerProperty(type=OmmatidiaSettings)


def unregister():
    import bpy

    from .ui import live

    live.cancel_pending()
    del bpy.types.Object.ommatidia
    for cls in reversed(_registered):
        bpy.utils.unregister_class(cls)
    _registered.clear()
