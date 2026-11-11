# module name: blender_io

**purpose:** the only bridge between Blender data and the pure modules, and
the single place where Blender units become specimen metres and back.

Code: [blender_io.py](blender_io.py) · uses: [lattice.md](../lattice/lattice.md),
[topology.md](../topology/topology.md) · axioms:
[data_definitions.md](../../../docs/agents/data_definitions.md)

## Introduction

[lattice](../lattice/lattice.md) and [topology](../topology/topology.md) never
touch Blender. They need a **SURFACE**, and their results must end up in the
scene as objects. This module does both conversions,
and nothing else in the add-on reads or writes mesh data.

**Reading.** The target object is evaluated (so the user's own modifiers, such
as Subdivision Surface, count), triangulated, and moved to world space. The
add-on's own "Ommatidia" modifier is switched off while reading, so the domes
never feed back into the surface they are placed on. Coordinates are divided
by `model_scale` (Blender units per real metre), so a 2-unit-wide sculpt of a
400 µm eye becomes a 400 µm surface. A BVH tree answers "closest point" queries.
**Thin solids** (e.g. Sculpt ▸ Mask Extract with "Extract as Solid"): when at
least a quarter of the area has an opposite-facing skin closer than about one
facet spacing, only the outer skin (faces pointing away from the shape's
centre) is kept. Two such skins cannot both carry ommatidia: the lattice would
hop between them. The number of dropped faces is reported in Diagnostics.
Normals come from the mesh's *corner* (split) normals interpolated across each
triangle: smooth on smooth-shaded meshes, exact face normals on flat faces. Open
edges are sampled into a KD-tree for "distance to the rim".

**Writing the ommatidia.** `<target>_Ommatidia` is the visible result: the
FACET_MESH of [facets.md](../facets/facets.md) as a real mesh (vertices, loops and
polygons set straight from numpy arrays), smooth-shaded, with material slot 0
on facets and 1 on bristles. Its per-vertex base data is stored as point
attributes, so a shape setting only reads them back, recomputes positions and
writes the vertex coordinates (`reshape_facets_object`). Regenerating swaps the
mesh but keeps the material slots.

**Writing the lattice.** A hidden mesh object, `<target>_OmmatidiaLattice`, holds the ommatidium
centres and the facet corners as vertices, and the facet walls as edges. Point
attributes tell which is which (`omm_kind`), where bristles go
(`omm_bristle`), each point's normal (`omm_normal`) and each centre's neighbour
count (`omm_valence`, handy for viewing defects). The object also stores the
CONFIG and metrics as JSON, so a .blend file records what made it.
Regenerating swaps the mesh in place; nothing piles up. Both objects are
parented to the target (world positions kept), so they follow it when it moves.

**Settings.** The UI panel edits `OmmatidiaSettings`, a property group on each
mesh object that mirrors CONFIG field by field. It shows µm where CONFIG stores
metres, and `props_to_config` / `config_to_props` are the only conversions.

## Data definitions

| DD | Type | Meaning |
|---|---|---|
| MODEL_SCALE | float | Blender units per specimen metre |
| MESH_ARRAYS | (V,3), (T,3), (T,3,3) | world-space vertices / model_scale, triangles, corner normals |
| SURFACE (`BvhSurface`) | object | BVH-backed SURFACE in specimen metres |
| LATTICE_OBJECT | bpy Object | centres + facet corners + walls, with attributes and JSON properties (hidden) |
| FACETS_OBJECT | bpy Object | the visible domed facets and bristles, with FACET_MESH attributes |
| OMMATIDIA_SETTINGS | PropertyGroup | the UI view of CONFIG on each object |

## Atomic algorithms

- **mesh_arrays** ([mesh_to_surface.py](atomic_algorithm_scripts/mesh_to_surface.py)):
  evaluated, triangulated mesh in world space / model_scale, with corner
  normals; skips our own modifier.
- **BvhSurface** (same file): `closest_point` (BVH + interpolated corner
  normals), `sample_area_weighted` (triangle-area CDF), `boundary_distance`
  (KD-tree of rim samples, refined exactly on the nearest rim segment), `area`.
- **write_lattice_object / read_lattice_metrics / lattice_object_name**
  ([lattice_to_object.py](atomic_algorithm_scripts/lattice_to_object.py)).
- **write_facets_object / read_facet_mesh / reshape_facets_object / facets_object_name**
  ([facets_to_object.py](atomic_algorithm_scripts/facets_to_object.py)).
- **props_to_config / config_to_props / lattice_origin_from_props**
  ([config_props.py](atomic_algorithm_scripts/config_props.py)). Settings carry
  update callbacks into [ui/live.py](../ui/live.py).

## Compound algorithm

```
target object ──mesh_arrays (skip "Ommatidia" modifier, ÷ model_scale)──► BvhSurface (SURFACE)
                                         │
                       lattice / topology (pure)
                                         ▼
LATTICE + TOPOLOGY + CONFIG + metrics ──write_lattice_object (× model_scale)──► <target>_OmmatidiaLattice
```

`surface_from_object(obj, config)` is the reading half; `write_lattice_object`
is the writing half. `main.py` and the Generate operator call them in that order.

Measured: reading a 164k-vertex mesh takes ~1 s; 800 ommatidia on it ~4 s in
total (relaxation dominates: ~160k BVH queries).
