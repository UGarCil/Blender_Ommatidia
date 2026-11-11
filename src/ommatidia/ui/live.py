"""Live updates: what happens when a setting changes in the panel.

- Shape settings (inflation, facet shape, bristle length, ...) move the existing
  vertices immediately; that takes milliseconds.
- Lattice settings (count, seed, centre, ...) need a new lattice, which takes
  about a second, so they regenerate after a short pause: dragging a slider
  regenerates once, when you let go, not at every step.
Nothing happens before the first Generate, or with Live update switched off.
"""

from __future__ import annotations

from contextlib import contextmanager

import bpy

DEBOUNCE_SECONDS = 0.4

_pending: set = set()     # names of targets waiting for regeneration
_paused = [False]          # set while the add-on itself writes settings


@contextmanager
def paused():
    """Suspend live updates while settings are written programmatically."""
    before = _paused[0]
    _paused[0] = True
    try:
        yield
    finally:
        _paused[0] = before


def _generated(target) -> bool:
    from ..blender_io.atomic_algorithm_scripts.facets_to_object import facets_object_name
    return bpy.data.objects.get(facets_object_name(target.name)) is not None


# FD. shape_changed
# Signature: bpy Object (target) -> None
def shape_changed(target) -> None:
    """Reshape the target's ommatidia object now (no-op before the first Generate)."""
    if _paused[0] or not _generated(target):
        return
    from ..blender_io.atomic_algorithm_scripts.config_props import props_to_config
    from ..pipeline import reshape

    props = target.ommatidia
    try:
        reshape(target, props_to_config(props, bpy.context.scene))
        props.status = ""
    except (ValueError, TypeError) as error:
        props.status = str(error).splitlines()[-1].strip(" -")


# FD. lattice_changed
# Signature: bpy Object (target) -> None
def lattice_changed(target) -> None:
    """Schedule a regeneration DEBOUNCE_SECONDS after the last change."""
    if _paused[0] or not target.ommatidia.live_update or not _generated(target):
        return
    _pending.add(target.name)
    if bpy.app.timers.is_registered(_regenerate_pending):
        bpy.app.timers.unregister(_regenerate_pending)
    bpy.app.timers.register(_regenerate_pending, first_interval=DEBOUNCE_SECONDS)


# FD. regenerate
# Signature: bpy Object (target) -> LATTICE_METRICS | None
def regenerate(target):
    """Run the whole pipeline for target; report problems in its status line."""
    from ..blender_io.atomic_algorithm_scripts.config_props import props_to_config
    from ..pipeline import generate

    props = target.ommatidia
    try:
        metrics = generate(target, props_to_config(props, bpy.context.scene))
    except (ValueError, TypeError) as error:
        props.status = str(error).splitlines()[-1].strip(" -")
        return None
    props.status = ""
    return metrics


def _regenerate_pending():
    names = sorted(_pending)
    _pending.clear()
    for name in names:
        target = bpy.data.objects.get(name)
        if target is not None:
            regenerate(target)
    return None   # run once


def cancel_pending() -> None:
    """Drop scheduled regenerations (on unregister)."""
    _pending.clear()
    if bpy.app.timers.is_registered(_regenerate_pending):
        bpy.app.timers.unregister(_regenerate_pending)
