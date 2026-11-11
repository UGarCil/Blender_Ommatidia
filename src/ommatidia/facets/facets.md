# module name: facets

**purpose:** turn the lattice into the visible eye: one connected mesh of domed
facets (corneal lenses) and interommatidial bristles, reshaped instantly when
a shape setting changes.

Code: [facets.py](facets.py) · inputs: [lattice.md](../lattice/lattice.md),
[topology.md](../topology/topology.md) · axioms:
[data_definitions.md](../../../docs/agents/data_definitions.md)

## Introduction

`topology/` gives every ommatidium its outline: the corners where three facets
meet and the walls between them. This module builds the surface the SEM shows.

**Composite mesh.** Each facet is a fan of rings from its centre out to its
walls (by default 18 points around, 6 rings). Corners and wall points are
shared by neighbouring facets, so the whole eye is one watertight sheet with no
cracks between lenses. Only facets whose ring of corners closes get a dome:
the rim row of an open sheet and cells next to a hole have none.

**Two stages.** Building the layout needs the surface (every vertex is
projected onto it) and happens once per Generate. Each vertex also remembers
how far it is from its facet's centre, both as a fraction of the way to the
wall (`rho_hex`, hexagonal) and in units of the facet's inner radius
(`rho_round`, circular). From that alone, `shape_positions` lifts every vertex
along its normal; no surface queries, so 100k vertices take ~20–30 ms and the
sliders feel immediate.

**Shape** (dossier §4.3): a dome over a round base of the facet's inner radius,
blended with a flat-topped hexagonal plateau.
- *Inflation*: dome height / facet radius; 0 flat, 1 a semi-dome (globular).
- *Pokeyness*: dome outline, round lens (0) to cone (1).
- *Facet shape*: 0 round dome (Drosophila), 1 flat hexagonal facet (wasps).
- *Crease sharpness*: width of the fall-off into the groove.
- *Rim fade*: facets flatten within a few spacings of an open edge.

**Bristles** stand in the facet corners chosen by `topology/` (one per
ommatidium in Drosophila). Each is a small spindle: thickest at a quarter of
its length, pointed, its base sunk slightly into the cornea. It tilts from the
normal toward the comb direction. Length and radius are in real µm; in COUNT
mode one lattice spacing counts as `spacing` µm, so bristles keep their size
relative to the facets on a mesh of any scale.

## Data definitions

| DD | Type | Meaning |
|---|---|---|
| CELL_LOOPS | list of (m,) int | corner triangles around each centre, CCW; empty = no facet |
| FACET_RESOLUTION | (sides, rings) | from SUBDIVISION_LEVEL L: (L+1, 2(L+1)) |
| FACET_LAYOUT | dataclass | unprojected vertices, owner, rho_hex, rho_round, faces |
| FACET_MESH | dataclass | per-vertex base data + faces; input of `shape_positions` |
| BRISTLE_TEMPLATE | arrays | one bristle in local coordinates (axial, angle, faces) |
| INFLATION, POKEYNESS, FACET_SHAPE, CREASE_SHARPNESS, RIM_WIDTH | floats | shape controls |

## Atomic algorithms

- **cell_loops** ([cell_loops.py](atomic_algorithm_scripts/cell_loops.py)): walk
  from triangle to triangle around each centre.
- **facet_resolution / facet_layout** ([dome_mesh.py](atomic_algorithm_scripts/dome_mesh.py)):
  vertices (shared corners and walls, then each facet's centre and rings) and
  faces (centre fan, quads between rings), vectorised per facet size.
- **dome_heights / rim_fade / crease_width** ([dome_profile.py](atomic_algorithm_scripts/dome_profile.py)).
- **bristle_template / spindle_radius / bristle_positions** ([bristles.py](atomic_algorithm_scripts/bristles.py)).

## Compound algorithm

```
build_facet_mesh(LATTICE, TOPOLOGY, SURFACE, CONFIG):
    cell_loops(TRIANGLES) ──facet_layout(sides, rings)──► FACET_LAYOUT
        ──SURFACE.closest_point──► base positions + normals
    SURFACE.boundary_distance(centres) / SPACING ──► rim distance per vertex
    bristle sites (TOPOLOGY.cell_points[bristle_mask]) × BRISTLE_TEMPLATE
    ──► FACET_MESH

shape_positions(FACET_MESH, CONFIG):
    facets:   base + normal · (0.01·SPACING + dome_heights(...) · rim_fade(...))
    bristles: bristle_positions(sites, tilt, comb, length, radius × real_length_scale)
    ──► POINTS
```

Sockets (a ring around each bristle base) are not built yet.
