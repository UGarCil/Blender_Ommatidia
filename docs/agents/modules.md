# Modular structure (axioms)

Status: **approved by Uriel 2026-11-03.** Data definitions referenced in
CAPITALS are in [data_definitions.md](data_definitions.md). Items marked
**[PROPOSAL]** were accepted (data_definitions.md §9); M1 additions are listed
in data_definitions.md §10. **[M3] changes (data_definitions.md §11) await
review:** `geometry_nodes/` is replaced by the pure `facets/` module, and the
lattice grows a perfect patch instead of relaxing random points. Module-level docs:
[lattice.md](../../src/ommatidia/lattice/lattice.md),
[topology.md](../../src/ommatidia/topology/topology.md),
[facets.md](../../src/ommatidia/facets/facets.md),
[material.md](../../src/ommatidia/material/material.md).

---

## 1. Overview

```
                 main.py  /  pipeline.py (Generate button, live updates)   (orchestrators)
                                   │
      ┌───────────────┬────────────┼──────────────┬──────────────────┐
      ▼               ▼            ▼              ▼                  ▼
  blender_io      lattice      topology        facets            material
 (bpy ⇄ numpy)  (pure numpy) (pure numpy)   (pure numpy)     (bpy shader tree)
      │               ▲            ▲              ▲
      └── SURFACE ────┴────────────┴──────────────┘
                      └─ LATTICE ──┴─ TOPOLOGY ───┘

  config.py (CONFIG) is imported by every module; it imports nothing from the package.
```

Data flow of one "Generate" run:

```
target mesh ──blender_io──► SURFACE ──lattice──► LATTICE ──topology──► TOPOLOGY
                                                    │                     │
                                                    └──────┬──────────────┘
                                                           ▼
                                  facets: FACET_MESH (domed facets + bristles)
                                                           ▼
                    blender_io: write FACETS_OBJECT (visible) and LATTICE_OBJECT (hidden data)
                                                           ▼
                                  material: SEM material on the facets object

A shape setting changes -> facets.shape_positions on the stored FACET_MESH
attributes -> blender_io moves the vertices (no SURFACE, ~20-30 ms).
```

Dependency rules

- `lattice/` and `topology/` import only `numpy`, the standard library and
  `config.py`. Never `bpy` or `mathutils`.
- `topology/` may import `lattice/` data definitions (LATTICE); `lattice/` never
  imports `topology/`.
- Only `blender_io/` reads or writes mesh data and object properties.
- `facets/` is pure numpy too (it may read LATTICE and TOPOLOGY).
- `material/` builds shader trees; it never touches mesh data.
- `ui/` calls the same compound functions as `main.py`; it holds no algorithm.

---

## 2. `lattice/`: centres on a surface

**Compound** `lattice.py`

```
# FD. build_lattice
# Signature: SURFACE, CONFIG -> LATTICE, LATTICE_METRICS (partial)
```

```
SURFACE, CONFIG
   │
   ├─ resolve_density ─────────► SPACING, N
   ├─ sample_candidates ───────► surface samples (30·N), for coverage and seams
   ├─ origin = LATTICE_ORIGIN or auto_origin
   ├─ repeat (≤ 8 patches):
   │     grow_patch: first_ring, then per ring
   │         ring_candidates ─► apex_triples ─► apex_points ─► mean_by_owner
   │         ─► closest_point ─► node_links ─► placement_mask (place or DEAD)
   │         ─► spring_step × 10 on the last 4 rings (patch_edges)
   │     nearest_distances(samples) ─► stop, or seed a new patch in the deepest gap
   ├─ fill_mask (seams) ───────► dart-thrown seam centres
   ├─ relax_seams: repulsion_step on seam centres, patch-edge and crowded centres
   ├─ trim_boundary ───────────► drop centres with boundary_distance < 0.5·SPACING
   ├─ jitter_tangent ──────────► IRREGULARITY, then re-project
   └─ lattice_metrics ─────────► MEAN_NN, SPACING_CV, N_PATCHES, N_SEAM
   ▼
LATTICE
```

**Atomic algorithms** (`atomic_algorithm_scripts/`)

| File | FD | Signature | Action |
|---|---|---|---|
| `density.py` | `resolve_density` | `SURFACE.area, CONFIG -> SPACING, TARGET_COUNT` | map: derive the missing one; raise if N > MAX_COUNT |
| `sampling.py` | `sample_candidates` | `SURFACE, int, RNG -> POINTS` | map: delegate to `surface.sample_area_weighted` |
| `sampling.py` | `poisson_filter` | `POINTS, float -> (K,) bool` | filter: greedy dart throwing on a hash grid |
| `sampling.py` | `fill_mask` [M3] | `POINTS, POINTS, NORMALS, NORMALS, float -> (K,) bool` | filter: dart throwing into the gaps of an existing set |
| `growth.py` [M3] | `auto_origin` | `POINTS, DISTANCES -> POINT` | evaluate: farthest from the rim, or toward -Y on closed meshes |
| `growth.py` | `tangent_frame`, `first_ring`, `hex_ring` | | map |
| `growth.py` | `ring_candidates` | `fronts, HEX_INDEX, ring -> nodes` | filter |
| `growth.py` | `apex_triples` | `nodes, HEX_INDEX -> owners, a, b, c` | filter: placed neighbour pairs |
| `growth.py` | `apex_points` | `POINTS, NORMALS, a, b, c, SPACING -> POINTS` | map: equilateral apex away from c |
| `growth.py` | `mean_by_owner` | | evaluate |
| `growth.py` | `node_links` | `nodes, HEX_INDEX, offset -> LATTICE_LINKS ×2` | map |
| `growth.py` | `placement_mask` | `new, proposals, …, LATTICE_LINKS -> (K,) bool` | filter: rim, off-surface, collisions (0.7·SPACING; 0.55 for lattice neighbours) |
| `growth.py` | `patch_edges` | `HEX_INDEX -> PATCH_EDGES` | map |
| `growth.py` | `spring_step` | `POINTS, VELOCITIES, NORMALS, PATCH_EDGES, mask, SPACING -> POINTS, VELOCITIES` | map: edges toward SPACING, tangential, momentum 0.8 |
| `neighbors.py` | `neighbor_pairs` | `POINTS, float -> NEIGHBOR_PAIRS` | filter: all pairs within radius via a numpy spatial hash (27 neighbouring cells) |
| `neighbors.py` | `assign_nearest` | `POINTS, POINTS, float -> (K,) int` | map: index of the nearest centre for each sample |
| `neighbors.py` | `nearest_distances` [M3] | `POINTS, POINTS, float -> DISTANCES` | map |
| `relax.py` | `repulsion_step` | `POINTS, VELOCITIES, NORMALS, NEIGHBOR_PAIRS, SPACING -> POINTS, VELOCITIES` | map: kernel force `(1 − d/r)²` with r = 2·SPACING, tangential only, momentum 0.9, move clamped to 0.15·SPACING (seams only since M3) |
| `relax.py` | `trim_boundary` | `POINTS, DISTANCES, SPACING -> mask` | filter |
| `relax.py` | `jitter_tangent` | `POINTS, NORMALS, IRREGULARITY, SPACING, RNG -> POINTS` | map |
| `metrics.py` | `nearest_neighbor_distances` | `POINTS, NEIGHBOR_PAIRS -> DISTANCES` | evaluate: min per point |
| `metrics.py` | `spacing_cv` | `DISTANCES, mask -> float` | count/evaluate |

The relaxation loop lives in the compound; each step is atomic. If profiling in
M1 shows the split costs real time, we stop and discuss (AGENTS §5.2).

---

## 3. `topology/`: facet outlines and bristle sites

**Compound** `topology.py`

```
# FD. build_topology
# Signature: LATTICE, SURFACE, CONFIG -> TOPOLOGY, LATTICE_METRICS (partial)
```

```
LATTICE
   ├─ neighbor_pairs (r = 2·SPACING) ───► candidate neighbours
   ├─ local_delaunay_fans ─────────────► per-centre triangle fans (tangent-plane Delaunay)
   ├─ agree_triangles ─────────────────► TRIANGLES (kept if all 3 vertices list it), oriented
   ├─ fill_small_holes [M3] ───────────► close 3- and 4-edge holes
   ├─ valence ─────────────────────────► VALENCE, HEX_FRACTION, N_PENTAGONS, N_HEPTAGONS
   ├─ circumcenters + surface.closest_point ─► CELL_VERTICES, MAX_HOLE
   ├─ dual_edges ──────────────────────► CELL_EDGES
   ├─ seam_cuts (pair 5/7 facets) ─────► SEAM_CUTS
   ├─ triangle_parity (BFS 2-colouring) ► TRIANGLE_PARITY, PARITY_CONFLICTS
   └─ select_bristles (BRISTLE_PATTERN) ► BRISTLE_MASK
   ▼
TOPOLOGY
```

| File | FD | Signature | Action |
|---|---|---|---|
| `triangulate.py` | `local_delaunay_fans` | `POINTS, NORMALS, NEIGHBOR_PAIRS, float -> (F,3)` | filter: empty-circumcircle triangles around each centre, tested in the plane ⟂ the mean corner normal; circumradius ≤ SPACING |
| `triangulate.py` | `agree_triangles` | `fans -> TRIANGLES` | filter: keep triangles proposed by all three vertices; orient CCW about the normal |
| `triangulate.py` | `fill_small_holes` [M3] | `TRIANGLES, POINTS, float -> TRIANGLES` | append: triangles for 3- and 4-edge hole loops |
| `triangulate.py` | `valence` | `TRIANGLES, N -> VALENCE` | count |
| `triangulate.py` | `open_edges` | `TRIANGLES -> (K,2)` | filter: edges used by one triangle (rim or hole) |
| `triangulate.py` | `overused_edges` | `TRIANGLES -> (K,2)` | filter: edges used by > 2 triangles (hard creases) |
| `dual.py` | `circumcenters` | `POINTS, TRIANGLES -> POINTS` | map |
| `dual.py` | `circumradii` | `POINTS, TRIANGLES, POINTS -> DISTANCES` | map |
| `dual.py` | `dual_edges` | `TRIANGLES -> CELL_EDGES` | map: triangles sharing an edge |
| `dual.py` | `crossed_edges` | `TRIANGLES -> (E,2)` | map: the two centres each cell edge separates |
| `bristle_sites.py` | `seam_cuts` | `(E,2), odd mask, rim mask -> SEAM_CUTS` | map: greedy shortest-path pairing of odd facets |
| `bristle_sites.py` | `triangle_parity` | `int, CELL_EDGES, SEAM_CUTS -> TRIANGLE_PARITY` | map: BFS 2-colouring around the seams + repair sweeps |
| `bristle_sites.py` | `count_conflicts` | `TRIANGLE_PARITY, CELL_EDGES -> int` | count |
| `bristle_sites.py` | `select_bristles` | `TRIANGLE_PARITY, BRISTLE_PATTERN, BRISTLE_FREQUENCY, RNG -> BRISTLE_MASK` | filter |

---

## 4. `blender_io/`: the only bridge to bpy data (and the unit boundary)

**Compound** `blender_io.py`: `surface_from_object`, `write_lattice_object`,
`read_config` / `write_config`.

| File | FD | Signature | Action |
|---|---|---|---|
| `mesh_to_surface.py` | `mesh_arrays` | `bpy Object, depsgraph, MODEL_SCALE, str -> MESH_ARRAYS` | map: evaluated mesh (own modifier skipped) → world space → ÷ MODEL_SCALE, with corner normals |
| `mesh_to_surface.py` | `BvhSurface(vertices, triangles, corner_normals)` | `MESH_ARRAYS -> SURFACE` | BVHTree for `closest_point` (interpolated corner normals); triangle-area CDF for sampling; KD-tree of rim samples + exact segment distance for `boundary_distance` |
| `lattice_to_object.py` | `write_lattice_object` | `bpy Object, LATTICE, TOPOLOGY, CONFIG, metrics -> LATTICE_OBJECT` | map: × MODEL_SCALE, `mesh.from_pydata`, write attributes and custom properties; reuse the object if it exists |
| `lattice_to_object.py` | `read_lattice_metrics` | `bpy Object -> LATTICE_METRICS \| None` | map: metrics JSON of the target's lattice object (for the UI) |
| `facets_to_object.py` [M3] | `write_facets_object` | `bpy Object, FACET_MESH, CONFIG -> FACETS_OBJECT` | map: arrays → mesh via `foreach_set`, attributes, parent to the target |
| `facets_to_object.py` | `read_facet_mesh` / `reshape_facets_object` | `FACETS_OBJECT (, CONFIG)` | map: attributes → FACET_MESH → new vertex positions |
| `config_props.py` | `props_to_config` / `config_to_props` | `PropertyGroup ⇄ CONFIG` | map: µm in the UI ⇄ metres in CONFIG |

---

## 5. `facets/`: the visible ommatidia [M3, replaces `geometry_nodes/`]

**Compound** `facets.py`: `build_facet_mesh(LATTICE, TOPOLOGY, SURFACE, CONFIG)
-> FACET_MESH` (once per Generate) and `shape_positions(FACET_MESH, CONFIG) ->
POINTS` (every shape change; no SURFACE).

| File | FD | Signature | Action |
|---|---|---|---|
| `cell_loops.py` | `cell_loops` | `TRIANGLES, N -> CELL_LOOPS` | map: corners around each centre, CCW |
| `dome_mesh.py` | `facet_resolution` | `SUBDIVISION_LEVEL -> (sides, rings)` | map |
| `dome_mesh.py` | `facet_layout` | `POINTS, POINTS, CELL_LOOPS, CELL_EDGES, sides, rings -> FACET_LAYOUT, radii` | map: rings from centre to wall, walls and corners shared |
| `dome_profile.py` | `dome_heights` | `rho_round, rho_hex, radius, INFLATION, POKEYNESS, FACET_SHAPE, CREASE_SHARPNESS -> heights` | evaluate (dossier §4.3) |
| `dome_profile.py` | `rim_fade` | `DISTANCES, RIM_WIDTH -> weights` | evaluate |
| `bristles.py` | `bristle_template`, `spindle_radius`, `bristle_positions` | | map |

Why not Geometry Nodes (Q4 in data_definitions.md §11): the same live response
(~20–30 ms for 100k vertices), tested with pytest, and no dependence on node
names that change between Blender versions.

---

## 6. `material/`: SEM look and corneal nipples

| File | Builds | Notes |
|---|---|---|
| `sem_look.py` | "Ommatidia SEM" material [M3] | Layer Weight facing → ramp → emission, with a little diffuse: the SEM edge effect. Compound `material.py`: `assign_materials` (two slots, keeps the user's) |
| `corneal_nipples.py` | "Omm Nipples" shader group | Voronoi (scale = 1/NIPPLE_SPACING in model units) → Bump with NIPPLE_HEIGHT; GEOMETRY mode deferred |

---

## 7. `ui/` and `pipeline.py`

`pipeline.py` [M3]: `generate(target, CONFIG)` runs SURFACE → LATTICE → TOPOLOGY →
FACET_MESH → objects → material; `reshape(target, CONFIG)` applies shape settings.
`ui/operators.py`: `ommatidia.generate`, `ommatidia.use_cursor`.
`ui/live.py`: property update callbacks: shape settings reshape at once,
lattice settings regenerate 0.4 s after the last change (`bpy.app.timers`).
`ui/panels.py`: Ommatidia (density, Generate, live update) with sub-panels
Lattice, Facets, Bristles, Diagnostics.

---

## 8. `main.py` (headless orchestrator)

```
blender -b input.blend --python main.py -- path/to/config.json
```

Statement order is the pipeline:

1. load CONFIG (`config.json` path after `--`, or the defaults below `if __name__ == "__main__":`), `validate()`
2. pick the target object (active, or `CONFIG`-named) → `surface_from_object`
3. `build_lattice` → `build_topology`
4. `build_facet_mesh` → `write_facets_object` → `assign_materials`
5. `write_lattice_object`
6. save `00_results/<timestamp>_<seed>/` with `scene.blend`, `config.json`, `metrics.json`

---

## 9. Tests

```
test/
    fixtures/analytic_surfaces.py   # SphereSurface, EllipsoidSurface, CapSurface (pure)
    pure/                           # pytest, any Python with numpy
        test_config.py              # invariants, JSON round trip, mutations
        test_surfaces.py            # fixtures obey the SURFACE invariants
        test_density.py, test_sampling.py, test_neighbors.py, test_relax.py,
        test_lattice.py             # regularity, coverage, determinism
        test_triangulate.py, test_dual.py, test_bristle_sites.py,
        test_topology.py            # hex fraction, no holes on closed meshes
    blender/
        run_blender_tests.py        # blender -b --factory-startup --python ...
        test_bvh_surface.py, test_lattice_object.py, test_node_groups.py
```

Order: DDs → atomics → compounds, as in AGENTS §8. Performance (800 ommatidia on
a 200k-vertex mesh) is logged, not asserted. Development setup: add `numpy` and
`pytest` to `.venv` and a `pyproject.toml` with the `src/` layout.

---

## 10. Milestone mapping

| Milestone | Modules touched |
|---|---|
| M1 | `config.py`, `lattice/`, test fixtures, pure tests; repulsion benchmark (Lloyd removed) |
| M1.5 [PROPOSAL] | `topology/` (pure, so it belongs with M1 before any Blender work) |
| M2 | `blender_io/`, minimal `ui/` operator, `main.py` |
| M3–M5 | `facets/` domes, facet shape, rim fade, bristles (sockets pending); live UI; lattice growth [M3] |
| M6 | `material/` (SEM look done; corneal nipples pending) |
| M7 | `ui/` polish, presets, `blender_manifest.toml`, `docs/user/` |
