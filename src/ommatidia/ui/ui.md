# module name: ui

**purpose:** the sidebar panels, operators and live updates. They hold no
algorithms; they call [pipeline.py](../pipeline.py), which runs the same steps
as [main.py](../../../main.py).

## Introduction

The panel lives in *View3D ▸ Sidebar (N) ▸ Ommatidia* and appears when the
active object is a mesh that can be covered (not an object made by the add-on).
Each mesh object keeps its own settings (`object.ommatidia`, see
[blender_io.md](../blender_io/blender_io.md)), so two eyes in one scene can
differ.

**Generate Ommatidia** (`ommatidia.generate`) runs the slow part: reading the
mesh, growing the lattice, outlines and bristle sites, building the facets, and
writing `<target>_Ommatidia` (visible) and `<target>_OmmatidiaLattice` (hidden
data). Invalid settings, or a density that would need more than `max_count`
ommatidia (usually a wrong `model_scale`), are reported instead of freezing
Blender.

**Live updates** ([live.py](live.py)). After the first Generate:
- shape settings (inflation, facet shape, pokeyness, crease, rim fade, bristle
  length/radius/tilt/comb) move the existing vertices at once (~20–30 ms);
- lattice settings (density, seed, centre, rotation, seam relaxation,
  irregularity, resolution, bristle pattern) regenerate 0.4 s after the last
  change, so dragging a slider regenerates once. *Live update* switches this
  off for very large lattices. Problems show as a status line in the panel.

**Centre** of the regular patch: automatic, the 3D cursor (read at each generation), an object
(e.g. an Empty), or a fixed point. *Use 3D Cursor* copies the cursor into a
fixed point.

Sub-panels: Lattice, Facets, Bristles, Diagnostics (count, spacing, hexagon
fraction, 5/7 defects, patches and seam fill, bristles, time; a warning when
hard creases broke facet edges).

## Atomic pieces

- [operators.py](operators.py): `OMMATIDIA_OT_generate`, `OMMATIDIA_OT_cursor_to_centre`, `is_coverable`.
- [live.py](live.py): `shape_changed`, `lattice_changed`, `regenerate`, `paused`.
- [panels.py](panels.py): `OMMATIDIA_PT_main` and its sub-panels.
