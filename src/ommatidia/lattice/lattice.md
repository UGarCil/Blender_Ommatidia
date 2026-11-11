# module name: lattice

**purpose:** place the ommatidium centres on any surface as a regular,
unstretched hexagonal lattice, with the unavoidable defects pushed to the edges.

Code: [lattice.py](lattice.py) · next module: [topology.md](../topology/topology.md) ·
axioms: [data_definitions.md](../../../docs/agents/data_definitions.md)

## Introduction

A compound eye is a curved sheet tiled by facets of almost equal size. To
imitate it on an arbitrary mesh we first need the *centre* of each facet, evenly
spaced **measured on the surface itself**, not in UV space.

The module never sees a Blender mesh. It receives a **SURFACE**: any object that
can project points onto itself, draw random points uniformly per unit area, say
how far a point is from an open edge, and report its area. Blender meshes and
analytic test shapes (sphere, ellipsoid, open cap) both provide one, so the
whole module is tested with plain `pytest`.

**Regularity first (since M3).** A hexagonal lattice cannot cover a curved
surface perfectly (Euler's formula; on a sphere at least twelve 5-sided cells).
Relaxing random points spreads the resulting 5/7 defects everywhere (about 10%
of cells at 750). Instead, a **perfect patch grows** from a centre, ring by
ring: each new centre sits at the apex of an equilateral triangle on two
placed neighbours, so every facet starts exactly one spacing wide. Curvature
compresses each new ring a little; a few **spring steps** along the lattice
edges after every ring spread that strain over the last rings instead of
letting it become defects. Where the patch would overlap itself (two fronts
meeting on the far side of a sphere) or reaches an open edge, the new centre is
not placed. The gaps left behind are the **seams**: they are filled by dart
throwing and smoothed by repulsion while the patch stays fixed. All defects
therefore sit in the seams, far from the centre, where a user can hide them.

The **centre** is chosen by the user (a point, the 3D cursor, an object) or
automatically: the point farthest from the rim on an open sheet (an eye cut out
of a head), the side facing −Y on a closed mesh. **Rotation** turns the rows.
Any large uncovered region (a second eye in the same object) gets its own patch.

Measured at 750 centres: a sphere is perfect within 62° of the centre and 95%
hexagonal overall (relaxation alone: 90%); an eye-like 60° cap is perfect up to
the last ~2 rows (96%). The price is a little size variation: cells shrink by
up to ~10% toward the seams (spacing CV ≈ 0.055 instead of 0.032).

On an **open sheet** the patch stops half a spacing from the edge (the rim
half-facets are dropped, ~10% of the count); the facets fade there later
(`facets/`, rim fade). Optional **irregularity** adds small random sideways
jitter at the end. The seed only changes the seams.

## Data definitions

| DD | Type | Meaning |
|---|---|---|
| SURFACE | object | the 2-manifold the lattice lives on |
| SPACING | float, m | target centre-to-centre distance (= facet flat-to-flat width) |
| LATTICE_ORIGIN, LATTICE_ROTATION | point or None, degrees | where the patch starts; how its rows turn |
| HEX_INDEX | dict (q, r) → row or DEAD | one patch in axial coordinates |
| PATCH_EDGES | (E,2) int | lattice-neighbour rows of a patch |
| LATTICE_LINKS | (L,2) int | new-node pairs that are lattice neighbours (may be closer) |
| NEIGHBOR_PAIRS | (M,2) int | all pairs closer than a radius |
| VELOCITIES | (N,3) float | previous move (momentum) in springs and repulsion |
| LATTICE | dataclass | centres, normals, spacing, patch id (−1 = seam) |

## Atomic algorithms

- **resolve_density** ([density.py](atomic_algorithm_scripts/density.py)): area
  + density choice → spacing and count; refuses counts above `max_count`.
- **growth** ([growth.py](atomic_algorithm_scripts/growth.py)):
  `auto_origin`, `tangent_frame`, `first_ring`, `ring_candidates`,
  `apex_triples`, `apex_points`, `mean_by_owner`, `node_links`,
  `placement_mask` (rim, off-surface, collision: 0.7·spacing, or 0.55 for
  lattice neighbours), `patch_edges`, `spring_step`.
- **sample_candidates / poisson_filter / fill_mask**
  ([sampling.py](atomic_algorithm_scripts/sampling.py)): surface samples; dart
  throwing into the gaps of an existing set.
- **neighbor_pairs / assign_nearest / nearest_distances**
  ([neighbors.py](atomic_algorithm_scripts/neighbors.py)): numpy spatial hash.
- **repulsion_step / trim_boundary / jitter_tangent**
  ([relax.py](atomic_algorithm_scripts/relax.py)): seam relaxation, rim filter,
  irregularity.
- **nearest_neighbor_distances / spacing_cv**
  ([metrics.py](atomic_algorithm_scripts/metrics.py)).

## Compound algorithm: `build_lattice(surface, config) -> (LATTICE, metrics)`

```
SURFACE.area, CONFIG ──resolve_density──► SPACING, N
SURFACE ──sample_candidates(30·N)──► samples (coverage, seams)
origin = LATTICE_ORIGIN or auto_origin(samples)
repeat (≤ 8 patches):
    grow_patch(origin): ring 1, then per ring k
        candidates ─► apex proposals ─► project ─► placement_mask ─► place / DEAD
        spring_step × 10 on rings k−4 … k
    nearest_distances(samples) ─► stop, or origin = deepest gap (≥ 2·SPACING)
uncovered samples ──fill_mask(0.85·SPACING)──► seam centres
relax_seams: repulsion on seam centres + patch edges + crowded centres (rest fixed)
trim_boundary (≥ 0.5·SPACING from an open edge) ──► jitter_tangent ──► LATTICE
──► metrics {N, N_TARGET, N_INTERIOR, SPACING, MEAN_NN, MIN_NN, SPACING_CV,
             N_PATCHES, N_SEAM, SECONDS}
```

Lloyd relaxation (2026-11-03) and global random relaxation (M3) were
benchmarked and removed; see the logbook.
