"""Atomic: the UI PropertyGroup is a *view* of CONFIG (AGENTS §5.3).

This is the only place where UI values and CONFIG are converted. The UI shows
µm where CONFIG stores metres; the lattice centre is chosen in the UI as
"automatic / 3D cursor / an object / a fixed point" and becomes
CONFIG.lattice_origin; everything else maps one to one.

Changing a setting updates the result live (ui/live.py): shape settings move the
existing vertices at once; lattice settings regenerate after a short pause.
"""

from __future__ import annotations

from dataclasses import fields

import bpy
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    FloatVectorProperty,
    IntProperty,
    PointerProperty,
    StringProperty,
)

from ...config import BRISTLE_PATTERNS, DENSITY_MODES, NIPPLE_MODES, Config

# CONFIG field -> factor from the UI value to CONFIG (µm -> m). Data only.
UI_TO_CONFIG = {
    "spacing": 1e-6,
    "bristle_length": 1e-6,
    "bristle_radius": 1e-6,
    "socket_radius": 1e-6,
    "socket_depth": 1e-6,
    "nipple_spacing": 1e-6,
    "nipple_height": 1e-6,
}
# CONFIG fields set in the UI through origin_source / origin_object / origin_point.
DERIVED_FIELDS = ("lattice_origin",)
ORIGIN_SOURCES = ("AUTO", "CURSOR", "OBJECT", "POINT")

_D = Config()  # defaults shown in the UI come from CONFIG, the single source of truth


def _enum(values, descriptions=None):
    return [(v, v.title().replace("_", " "), (descriptions or {}).get(v, "")) for v in values]


def _um(field):
    return getattr(_D, field) / UI_TO_CONFIG[field]


def _shape_changed(self, context):
    from ...ui import live
    live.shape_changed(self.id_data)


def _lattice_changed(self, context):
    from ...ui import live
    live.lattice_changed(self.id_data)


# CDD. OMMATIDIA_SETTINGS
# settings = OmmatidiaSettings (bpy PropertyGroup on each mesh Object)
# interp. the UI view of CONFIG for that object; same field names as CONFIG,
#         plus UI-only fields: live_update, origin_source/object/point, status
class OmmatidiaSettings(bpy.types.PropertyGroup):
    # density
    density_mode: EnumProperty(name="Density", items=_enum(DENSITY_MODES, {
        "SPACING": "Set the facet spacing; the count follows from the surface area",
        "COUNT": "Set the number of ommatidia; the spacing follows from the surface area"}),
        default=_D.density_mode, update=_lattice_changed)
    spacing: FloatProperty(name="Spacing (µm)", min=1e-6, default=_um("spacing"), update=_lattice_changed,
                           description="Centre-to-centre distance between neighbouring facets, in real µm")
    target_count: IntProperty(name="Count", min=1, default=_D.target_count, update=_lattice_changed,
                              description="Number of ommatidia (open sheets lose ~10% to the rim)")
    max_count: IntProperty(name="Max count", min=1, default=_D.max_count,
                           description="Safety cap against accidentally huge lattices")
    model_scale: FloatProperty(name="Model scale", min=1e-12, default=_D.model_scale, update=_lattice_changed,
                               description="Blender units per real metre of specimen "
                                           "(e.g. a 400 µm eye modelled 2 units wide: 5000)")
    # lattice
    seed: IntProperty(name="Seed", min=0, default=_D.seed, update=_lattice_changed,
                      description="Random choices in the seams (the regular patch does not change)")
    origin_source: EnumProperty(name="Centre", items=_enum(ORIGIN_SOURCES, {
        "AUTO": "Open sheet: the point farthest from its edge. Closed mesh: the side facing -Y (front view)",
        "CURSOR": "The 3D cursor (Shift + right-click on the mesh to place it)",
        "OBJECT": "The location of an object, e.g. an Empty",
        "POINT": "A fixed point (world coordinates)"}),
        default="AUTO", update=_lattice_changed,
        description="Where the perfectly regular patch starts; defects collect far from it")
    origin_object: PointerProperty(name="Centre object", type=bpy.types.Object, update=_lattice_changed)
    origin_point: FloatVectorProperty(name="Centre point", size=3, subtype="TRANSLATION",
                                      update=_lattice_changed)
    lattice_rotation_deg: FloatProperty(name="Rotation (°)", default=_D.lattice_rotation_deg,
                                        soft_min=-180.0, soft_max=180.0, update=_lattice_changed,
                                        description="Turns the rows of facets around the centre")
    relax_iterations: IntProperty(name="Seam relaxation", min=0, default=_D.relax_iterations,
                                  update=_lattice_changed,
                                  description="Smoothing iterations for the irregular seams")
    irregularity: FloatProperty(name="Irregularity", min=0.0, max=0.3, default=_D.irregularity,
                                update=_lattice_changed,
                                description="Random sideways jitter, as a fraction of the spacing")
    # facet shape (live)
    inflation: FloatProperty(name="Inflation", min=0.0, max=1.0, default=_D.inflation, update=_shape_changed,
                             description="Dome height as a fraction of the facet radius: "
                                         "0 = flat, 1 = semi-dome (globular)")
    facet_shape: FloatProperty(name="Facet shape", min=0.0, max=1.0, default=_D.facet_shape,
                               update=_shape_changed,
                               description="0 = round dome, 1 = flat-topped hexagonal facet")
    pokeyness: FloatProperty(name="Pokeyness", min=0.0, max=1.0, default=_D.pokeyness, update=_shape_changed,
                             description="Dome outline: 0 = round lens, 1 = pointed cone")
    crease_sharpness: FloatProperty(name="Crease sharpness", min=0.0, max=1.0, default=_D.crease_sharpness,
                                    update=_shape_changed,
                                    description="How narrow the groove between facets is")
    rim_width: FloatProperty(name="Rim fade", min=0.0, default=_D.rim_width, update=_shape_changed,
                             description="Facets flatten over this many spacings at an open edge")
    subdivision_level: IntProperty(name="Resolution", min=-1, max=6, default=_D.subdivision_level,
                                   update=_lattice_changed,
                                   description="Mesh detail per facet; -1 = automatic (2)")
    # bristles
    bristle_pattern: EnumProperty(name="Bristles", items=_enum(BRISTLE_PATTERNS, {
        "ALTERNATING": "Every other facet corner: about one bristle per ommatidium (Drosophila)",
        "RANDOM": "Each corner with probability 'Frequency'",
        "NONE": "No bristles"}), default=_D.bristle_pattern, update=_lattice_changed)
    bristle_frequency: FloatProperty(name="Frequency", min=0.0, max=1.0, default=_D.bristle_frequency,
                                     update=_lattice_changed)
    bristle_length: FloatProperty(name="Length (µm)", min=1e-6, default=_um("bristle_length"),
                                  update=_shape_changed)
    bristle_radius: FloatProperty(name="Radius (µm)", min=1e-6, default=_um("bristle_radius"),
                                  update=_shape_changed)
    bristle_tilt_deg: FloatProperty(name="Tilt (°)", min=0.0, max=90.0, default=_D.bristle_tilt_deg,
                                    update=_shape_changed,
                                    description="Angle between a bristle and the surface normal")
    comb_direction: FloatVectorProperty(name="Comb direction", size=3, default=_D.comb_direction,
                                        update=_shape_changed,
                                        description="Bristles lean toward this direction (world axes)")
    socket_radius: FloatProperty(name="Socket radius (µm)", min=0.0, default=_um("socket_radius"))
    socket_depth: FloatProperty(name="Socket depth (µm)", min=0.0, default=_um("socket_depth"))
    # corneal nipples
    nipple_mode: EnumProperty(name="Corneal nipples", items=_enum(NIPPLE_MODES), default=_D.nipple_mode)
    nipple_spacing: FloatProperty(name="Nipple spacing (µm)", min=1e-6, default=_um("nipple_spacing"))
    nipple_height: FloatProperty(name="Nipple height (µm)", min=0.0, default=_um("nipple_height"))
    # UI only
    live_update: BoolProperty(name="Live update", default=True,
                              description="Regenerate automatically when a lattice setting changes "
                                          "(shape settings always update at once)")
    status: StringProperty(name="Status", default="")


# FD. lattice_origin_from_props
# Signature: OMMATIDIA_SETTINGS, bpy Scene, MODEL_SCALE -> LATTICE_ORIGIN
def lattice_origin_from_props(props, scene, model_scale: float):
    """The lattice centre in specimen metres (world axes), or None for automatic."""
    source = props.origin_source
    if source == "CURSOR" and scene is not None:
        location = scene.cursor.location
    elif source == "OBJECT" and props.origin_object is not None:
        location = props.origin_object.matrix_world.translation
    elif source == "POINT":
        location = props.origin_point
    else:
        return None
    return tuple(float(x) / model_scale for x in location)


# FD. props_to_config
# Signature: OMMATIDIA_SETTINGS, bpy Scene | None -> CONFIG
def props_to_config(props, scene=None) -> Config:
    """Read every CONFIG field from the UI, converting µm to metres, and validate.

    Algorithm:
        for each CONFIG field: value = UI value × UI_TO_CONFIG factor (1 if none)
        tuples for vector fields; lattice_origin from the centre settings
        build Config; validate
    """
    values = {}
    for f in fields(Config):
        if f.name in DERIVED_FIELDS:
            continue
        value = getattr(props, f.name)
        if f.name in UI_TO_CONFIG:
            value = value * UI_TO_CONFIG[f.name]
        elif f.name == "comb_direction":
            value = tuple(float(x) for x in value)
        values[f.name] = value
    values["lattice_origin"] = lattice_origin_from_props(props, scene, props.model_scale)
    config = Config(**values)
    config.validate()
    return config


# FD. config_to_props
# Signature: CONFIG, OMMATIDIA_SETTINGS -> None
def config_to_props(config: Config, props) -> None:
    """Write CONFIG into the UI, converting metres to µm (without triggering live updates)."""
    from ...ui import live

    with live.paused():
        for f in fields(Config):
            if f.name in DERIVED_FIELDS:
                continue
            value = getattr(config, f.name)
            if f.name in UI_TO_CONFIG:
                value = value / UI_TO_CONFIG[f.name]
            setattr(props, f.name, value)
        if config.lattice_origin is None:
            props.origin_source = "AUTO"
        else:
            props.origin_source = "POINT"
            props.origin_point = tuple(x * config.model_scale for x in config.lattice_origin)
