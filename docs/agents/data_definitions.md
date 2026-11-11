# Data definitions (axioms)

Status: **approved by Uriel 2026-11-03** (P1–P9 accepted). Later changes are
marked *[M1]*/*[M2]* and listed in §10. **[M3] changes (§11) follow Uriel's
feedback of 2026-11-03 and await review.**

This file lists every piece of information the program handles, from atomic
values up to the compound objects that the modules pass to each other. Code
copies these blocks verbatim (`# DD.` / `# CDD.`) above the first place each
definition is used. Items marked **[PROPOSAL]** differ from, or add to,
`AGENTS.md`; they are summarised in §9.

Conventions

- `N` = number of ommatidia, `T` = number of triangles, `V` = number of cell
  vertices, `E` = number of cell edges, `B` = number of bristles.
- Arrays are `numpy.ndarray`, C-contiguous, `float64` for coordinates and
  `int64` for indices, unless stated.
- `#stub` marks the base/degenerate example; it is also the first test case.

---

## 1. Units [PROPOSAL]

Users sculpt at any scale: an eye mesh 2 Blender units (BU) wide may stand for a
400 µm eye. AGENTS.md says "SI in code"; this makes it precise.

1. **Pure modules (`lattice/`, `topology/`) work in specimen SI metres.** A
   16 µm facet is `16e-6`.
2. **`blender_io/` is the single unit boundary.** It divides mesh coordinates
   by `MODEL_SCALE` when it builds a SURFACE, and multiplies by it when it
   writes arrays back to Blender or sets GN modifier inputs.
3. **Lengths that belong to the facet are stored relative to it**: rim width
   in multiples of `SPACING`, dome height as a fraction of the facet's inner
   radius (`INFLATION`), crease width as a fraction of the centre-to-wall
   distance. They then stay correct when density changes. Lengths of independent
   structures (bristles, sockets, nipples) are stored in real SI metres. *[M3]*
   In `COUNT` mode the specimen's real size is unknown, so one lattice spacing is
   taken to be `CONFIG.spacing` real metres when real lengths are applied
   (`real_length_scale`); otherwise bristles on a unit-sized mesh at
   `model_scale` 1 would be microscopic.

```
# DD. MODEL_SCALE
# model_scale = float
# interp. Blender units per specimen metre: how much bigger the mesh is than
#         the real eye. 1.0 means the mesh is modelled at real size.
# Invariant: model_scale > 0
model_scale0 = 1.0 #stub: mesh at real size (a 400 µm eye is 0.0004 BU wide)
model_scale1 = 5000.0  # a 400 µm eye sculpted 2 BU wide
```

---

## 2. Atomic data

### 2.1 Geometry primitives

```
# DD. POINTS
# points = np.ndarray, shape (K, 3), float64
# interp. K positions in specimen metres, world orientation
# Invariant: finite values; K >= 0
points0 = np.zeros((0, 3)) #stub: empty
points1 = np.array([[0.0, 0.0, 1e-4]])

# DD. NORMALS
# normals = np.ndarray, shape (K, 3), float64
# interp. outward unit surface normals, row k belongs to POINTS row k
# Invariant: same K as its POINTS; every row has norm 1 (±1e-9)
normals0 = np.zeros((0, 3)) #stub
normals1 = np.array([[0.0, 0.0, 1.0]])

# DD. DISTANCES
# distances = np.ndarray, shape (K,), float64
# interp. a non-negative length per point, in specimen metres
# Invariant: >= 0; may contain np.inf (e.g. boundary distance on a closed mesh)
distances0 = np.zeros(0) #stub
distances1 = np.array([np.inf, 3e-6])
```

### 2.2 Density

```
# DD. SPACING
# spacing = float
# interp. target centre-to-centre distance between neighbouring ommatidia, in
#         specimen metres. For a regular hexagonal facet this equals the facet's
#         flat-to-flat width ("facet diameter" in dossier §3).
# Invariant: spacing > 0
spacing0 = 1.0 #stub: unit spacing, used with unit test surfaces
spacing1 = 16e-6  # Drosophila default (dossier §3, verify)

# DD. TARGET_COUNT
# target_count = int
# interp. desired number of ommatidia covering the surface
# Invariant: target_count >= 1
target_count0 = 1 #stub
target_count1 = 750  # Drosophila: ~700 (male) to ~750 (female) per eye,
                     # "An anatomical atlas of Drosophila melanogaster - the wild-type"

# DD. DENSITY_MODE
# density_mode = "SPACING" | "COUNT"
# interp. which of SPACING / TARGET_COUNT the user set; the other is derived:
#         count = round(area / (sqrt(3)/2 * spacing**2))
#         spacing = sqrt(2 * area / (sqrt(3) * count))
# (sqrt(3)/2 * spacing**2 is the area of one hexagonal cell.)
density_mode0 = "SPACING" #stub
density_mode1 = "COUNT"

# DD. MAX_COUNT
# max_count = int
# interp. safety cap; a derived count above it raises a clear error instead of
#         freezing Blender (e.g. a real-size spacing on a mesh with model_scale 1
#         that is actually 2 m wide)
# Invariant: max_count >= 1
max_count0 = 1 #stub
max_count1 = 50_000
```

Sanity check of the defaults: 750 cells of 16 µm cover
`750 · 0.866 · (16e-6)² ≈ 1.7e-7 m²`, a hemisphere of radius ≈ 160 µm. That is a
plausible *Drosophila* eye.

### 2.3 Randomness

```
# DD. SEED
# seed = int
# interp. the single source of randomness for a run
# Invariant: seed >= 0. Same SEED + same CONFIG + same SURFACE -> identical output
seed0 = 0 #stub
seed1 = 12345

# DD. RNG
# rng = np.random.Generator
# interp. a generator for one pipeline stage, spawned from SEED with
#         np.random.SeedSequence(seed).spawn(n) in a fixed stage order
#         (sampling, relaxation samples, jitter, random bristles), so adding
#         randomness to one stage never changes the others
rng0 = np.random.default_rng(0) #stub
```

### 2.4 Shape controls (dimensionless)

```
# DD. FACET_SHAPE
# facet_shape = float
# interp. blend between profiles: 0 = rounded dome (Drosophila), 1 = flat
#         hexagonal facet (many Hymenoptera)
# Invariant: 0 <= facet_shape <= 1
facet_shape0 = 0.0 #stub: pure dome
facet_shape1 = 1.0  # pure flat facet

# DD. INFLATION  [M3, PROPOSAL]
# inflation = float
# interp. dome height as a fraction of the facet's inner radius: 0 = flat,
#         1 = a semi-dome (hemisphere; highly globular lenses)
# Invariant: 0 <= inflation <= 1
inflation0 = 0.0 #stub: flat
inflation1 = 0.4   # Drosophila default, a guess from SEM (dossier §8)

# DD. POKEYNESS  [M3, PROPOSAL: was "height / SPACING"; height is now INFLATION]
# pokeyness = float
# interp. outline of the dome: 0 = round lens (spherical cap), 1 = pointed cone
# Invariant: 0 <= pokeyness <= 1
pokeyness0 = 0.0 #stub: round (Drosophila)
pokeyness1 = 1.0   # conical

# DD. LATTICE_ORIGIN  [M3]
# lattice_origin = (x, y, z) | None
# interp. where the perfect hexagonal patch starts, in specimen metres (world
#         axes); projected onto the surface. None = automatic: the point
#         farthest from the rim on an open sheet, the surface point facing -Y
#         (Blender's front view) on a closed mesh
lattice_origin0 = None #stub: automatic
lattice_origin1 = (0.0, -170e-6, 0.0)

# DD. LATTICE_ROTATION  [M3]
# lattice_rotation_deg = float
# interp. turns the rows of facets about the normal at the origin; 0 = one
#         lattice axis along world +Z projected onto the surface (+X if the
#         normal is vertical)
lattice_rotation0 = 0.0 #stub

# DD. CREASE_SHARPNESS
# crease_sharpness = float
# interp. how narrow the groove between facets is. The flat-facet fall-off width
#         is w = (1 - crease_sharpness) * SPACING / 2 (plus a tiny epsilon)
# Invariant: 0 <= crease_sharpness <= 1
crease_sharpness0 = 0.0 #stub: widest, softest groove
crease_sharpness1 = 0.8

# DD. IRREGULARITY
# irregularity = float
# interp. tangential jitter applied after relaxation, as a fraction of SPACING
# Invariant: 0 <= irregularity <= 0.3 (above that the lattice stops reading as hexagonal)
irregularity0 = 0.0 #stub: perfect lattice
irregularity1 = 0.03

# DD. RIM_WIDTH
# rim_width = float
# interp. width of the fade at an open boundary, in multiples of SPACING
# Invariant: rim_width >= 0 (0 = hard edge)
rim_width0 = 0.0 #stub
rim_width1 = 1.5
```

### 2.5 Bristles

```
# DD. BRISTLE_PATTERN
# bristle_pattern = "ALTERNATING" | "RANDOM" | "NONE"
# interp. which cell vertices (where 3 facets meet) carry a bristle:
#   ALTERNATING: one of the two triangle orientations, i.e. every other corner of
#                each hexagon -> about one bristle per ommatidium (Drosophila:
#                "each accompanied by a single mechanosensory bristle", atlas above)
#   RANDOM:      each cell vertex independently, probability BRISTLE_FREQUENCY
#   NONE:        no bristles
bristle_pattern0 = "NONE" #stub
bristle_pattern1 = "ALTERNATING"
bristle_pattern2 = "RANDOM"

# DD. BRISTLE_FREQUENCY
# bristle_frequency = float
# interp. fraction of cell vertices with a bristle. Used by RANDOM only;
#         ALTERNATING is 0.5 by construction (UI greys the slider out).
# Invariant: 0 <= bristle_frequency <= 1
bristle_frequency0 = 0.0 #stub
bristle_frequency1 = 0.5

# DD. BRISTLE_GEOMETRY  (SI metres / degrees; all biology defaults are guesses
#                        from SEM panel b, listed in dossier §8)
# bristle_length    = float > 0      # 10e-6
# bristle_radius    = float > 0      # base radius, 0.4e-6
# bristle_tilt_deg  = float in [0, 90]  # angle from the surface normal, 35
# comb_direction    = (3,) float, non-zero  # world vector; projected onto each
#                                           # tangent plane to set tilt direction
# socket_radius     = float >= 0     # 1.0e-6 (0 = no socket)
# socket_depth      = float >= 0     # 0.3e-6
```

### 2.6 Corneal nipples

```
# DD. NIPPLE_MODE
# nipple_mode = "OFF" | "SHADER" | "GEOMETRY"
# interp. how corneal nipples are made: none, shader bump (default, cheap), or
#         GN micro-displacement (expensive, close-ups only)
nipple_mode0 = "OFF" #stub
nipple_mode1 = "SHADER"

# DD. NIPPLE_GEOMETRY  (SI metres; values need citations, dossier §8)
# nipple_spacing = float > 0   # 0.2e-6
# nipple_height  = float >= 0  # 0.05e-6
```

---

## 3. SURFACE (the purity boundary)

```
# CDD. SURFACE
# surface = Surface(...)  any object exposing:
#   surface.closest_point(points: POINTS) -> (POINTS, NORMALS)
#       projection of each point onto the surface and the normal there
#   surface.sample_area_weighted(n: int, rng: RNG) -> POINTS
#       n points, uniform per unit area
#   surface.boundary_distance(points: POINTS) -> DISTANCES
#       geodesic-or-Euclidean distance to the nearest open edge; inf on closed surfaces
#   surface.area() -> float                                  [PROPOSAL: added]
#       total area in m^2; needed to turn SPACING into a count and back
# interp. any 2-manifold the lattice lives on, real (bpy mesh) or analytic.
#         Coordinates are specimen metres (see §1).
# Invariant: closest_point is idempotent (projecting twice = projecting once);
#            normals are unit length and point outward.
# Implementations:
#   BvhSurface      (blender_io; mathutils BVHTree on the evaluated mesh, world space / MODEL_SCALE)
#   SphereSurface   (test fixture; radius r)
#   EllipsoidSurface(test fixture; semi-axes a, b, c; Newton projection)
#   CapSurface      (test fixture; sphere radius r, polar angle limit theta_max; open)
surface0 = SphereSurface(radius=1.0) #stub
surface1 = EllipsoidSurface(1.0, 0.8, 0.6)
surface2 = CapSurface(radius=1.0, theta_max=np.pi / 3)
```

Template for SURFACE ("pass the object, not its attributes"):

```
... (surface)
... surface.closest_point(points)
... surface.sample_area_weighted(n, rng)
... surface.boundary_distance(points)
... surface.area()
```

---

## 4. Lattice data (output of `lattice/`)

```
# DD. NEIGHBOR_PAIRS                                         [PROPOSAL: new]
# neighbor_pairs = np.ndarray, shape (M, 2), int64
# interp. every unordered pair (i, j) of POINTS closer than a cut-off radius
# Invariant: i < j in every row; no duplicate rows; rows sorted lexicographically
neighbor_pairs0 = np.zeros((0, 2), dtype=np.int64) #stub
neighbor_pairs1 = np.array([[0, 1], [0, 2], [1, 2]])

# DD. RELAX_METHOD  [M1: REMOVED 2026-11-03, Uriel. Only repulsion remains; kept here as history]
# relax_method = "REPULSION" | "LLOYD"
# interp. REPULSION: Turk-style short-range point repulsion, projected to SURFACE.
#         LLOYD: sample-based centroidal Voronoi iteration (dense area-weighted
#         samples assigned to the nearest centre; centre moves to the mean, then
#         is projected back). [M1] Benchmark: REPULSION (with momentum) wins clearly
#         (CV 0.032 vs 0.106, hexagons 90% vs 52%, 10x faster). LLOYD is kept only
#         until Uriel decides whether to remove it.
relax_method0 = "REPULSION" #stub
relax_method1 = "LLOYD"

# DD. RELAX_ITERATIONS
# relax_iterations = int
# Invariant: relax_iterations >= 0 (0 = keep the initial sample)
relax_iterations0 = 0 #stub
relax_iterations1 = 200  # default [M1]

# DD. VELOCITIES  [M1]
# velocities = np.ndarray, shape (N, 3), float64
# interp. the previous relaxation move of each centre, carried over as momentum
velocities0 = np.zeros((0, 3)) #stub

# DD. HEX_INDEX, PATCH_EDGES, LATTICE_LINKS  [M3] (see lattice/growth.py)
# hex_index = dict[(q, r) -> row | DEAD]: one patch's nodes in axial coordinates
# patch_edges = (E, 2) int64: rows that are lattice neighbours in a patch
# lattice_links = (L, 2) int64: (new node, other point) pairs that are lattice
#         neighbours, so they may be closer than the collision distance

# DD. PATCH_IDS  [M3]
# patch = np.ndarray, shape (N,), int64
# interp. the patch each centre grew in (0, 1, ...), or SEAM (-1) for a centre
#         dart-thrown into a seam between patches

# CDD. LATTICE
# lattice = Lattice(centers: POINTS, normals: NORMALS, spacing: SPACING, patch: PATCH_IDS)
# interp. the ommatidium centres on a SURFACE: perfect hexagonal patches grown
#         from LATTICE_ORIGIN [M3], seams filled and relaxed, boundary trimmed,
#         jittered
# Invariant: centres lie on the surface (closest_point(centers) == centers ± 1e-9·spacing);
#            no two centres closer than 0.5·spacing;
#            on open surfaces every centre has boundary_distance >= 0.5·spacing  [PROPOSAL]
lattice0 = Lattice(np.zeros((0, 3)), np.zeros((0, 3)), 1.0) #stub: empty
lattice1 = Lattice(centers=<~800 points on a 170 µm hemisphere>, normals=..., spacing=16e-6)
```

Template for LATTICE:

```
... (lattice)
... lattice.centers
... lattice.normals
... lattice.spacing
```

---

## 5. Topology data (output of `topology/`)

```
# DD. TRIANGLES
# triangles = np.ndarray, shape (T, 3), int64
# interp. triangulation of the LATTICE centres on the surface; each row indexes
#         LATTICE.centers. Vertex order is counter-clockwise seen from outside
#         (along the normal).
# Invariant: every edge is shared by at most 2 triangles; on a closed surface by
#            exactly 2 (no holes); no triangle has a circumradius above SPACING
triangles0 = np.zeros((0, 3), dtype=np.int64) #stub
triangles1 = np.array([[0, 1, 2], [0, 2, 3]])

# DD. VALENCE
# valence = np.ndarray, shape (N,), int64
# interp. number of triangulation neighbours of each centre = number of sides
#         of its cell. 6 = hexagon; 5 and 7 are lattice defects.
valence0 = np.zeros(0, dtype=np.int64) #stub
valence1 = np.array([6, 6, 5, 7, 6])

# DD. CELL_VERTICES
# cell_vertices = (POINTS, NORMALS) with K = T
# interp. one point per triangle: its circumcentre projected onto SURFACE.
#         These are the corners of the hexagonal facets, where 3 facets meet.

# DD. CELL_EDGES
# cell_edges = np.ndarray, shape (E, 2), int64
# interp. the facet boundaries: an edge joins the cell vertices of two triangles
#         that share a triangulation edge (the dual graph)
# Invariant: i < j; indices < T; no duplicates
cell_edges0 = np.zeros((0, 2), dtype=np.int64) #stub
cell_edges1 = np.array([[0, 1]])

# DD. TRIANGLE_PARITY
# triangle_parity = np.ndarray, shape (T,), int8, values in {0, 1}
# interp. 2-colouring of triangles so neighbours differ (the "up" and "down"
#         triangles of a hex lattice). Exact on a perfect lattice; around 5/7
#         defects an odd cycle forces some same-colour neighbours, which are
#         counted as PARITY_CONFLICTS.

# DD. SEAM_CUTS  [M1]
# seam_cuts = np.ndarray, shape (E,), bool, aligned with CELL_EDGES
# interp. facet walls allowed to break the alternation: short seams that pair
#         each 5-/7-sided facet with a nearby one (or with the open rim), so the
#         rest of the pattern alternates perfectly

# DD. BRISTLE_MASK
# bristle_mask = np.ndarray, shape (T,), bool
# interp. True where the cell vertex of that triangle carries a bristle

# CDD. TOPOLOGY
# topology = Topology(triangles: TRIANGLES, valence: VALENCE,
#                     cell_points: POINTS, cell_normals: NORMALS,
#                     cell_edges: CELL_EDGES, parity: TRIANGLE_PARITY,
#                     bristle_mask: BRISTLE_MASK)
# interp. everything derived from the centres' connectivity: facet outlines and
#         bristle sites
topology0 = Topology(<all empty>) #stub
```

---

## 6. Metrics

```
# DD. LATTICE_METRICS
# lattice_metrics = {
#   "N":                  int,    # number of centres
#   "N_TARGET":           int,    # count before rim trimming  [M1]
#   "N_INTERIOR":         int,    # centres > SPACING from an open edge  [M1]
#   "MIN_NN":             float,  # smallest nearest-neighbour distance  [M1]
#   "SPACING":            float,  # target spacing used (m)
#   "MEAN_NN":            float,  # mean nearest-neighbour distance (m)
#   "SPACING_CV":         float,  # std/mean of nearest-neighbour distances, interior centres only
#   "HEX_FRACTION":       float,  # share of interior centres with VALENCE == 6
#   "N_PENTAGONS":        int,    # interior centres with valence 5
#   "N_HEPTAGONS":        int,    # interior centres with valence 7
#   "N_OTHER_DEFECTS":    int,    # interior valence < 5 or > 7  [M1]
#   "N_OPEN_EDGES":       int,    # triangulation rim (open sheets) or holes  [M1]
#   "HOLES_ON_CLOSED_SURFACE": int,  # open edges when the SURFACE is closed  [M1]
#   "N_OVERLAP_EDGES":    int,    # edges used by > 2 triangles (hard creases only)  [M2]
#   "MAX_HOLE":           float,  # largest triangle circumradius / SPACING
#                                 #   (perfect lattice = 1/sqrt(3) ≈ 0.577; > 1 = a hole)
#   "PARITY_CONFLICTS":   int,
#   "N_BRISTLES":         int,
#   "SECONDS":            {"sampling": float, "relax": float, "topology": float}
# }
# interp. quality report written to 00_results/<run>/metrics.json and to the
#         lattice object's custom properties
# "Interior" = boundary_distance > SPACING (all centres on a closed surface), so
# the rim does not pollute the regularity numbers. [M1] For the topology metrics,
# interior also excludes centres touching an open triangulation edge (cut-open cells).
# Acceptance (dossier §5): SPACING_CV <= 0.05 sphere, <= 0.08 deformed ellipsoid;
# HEX_FRACTION >= 0.90; MAX_HOLE < 1.0.
# [M1] Tests assert HEX_FRACTION >= 0.88: a 750-centre sphere measures 0.896 and
# that is physics, not a bug (extra 5-7 scars, Bausch et al. 2003). 5000 centres: 0.92.
```

Dictionaries hold data only, never functions (house rule).

---

## 7. CONFIG

```
# CDD. CONFIG
# config = Config(...)  a dataclass in src/ommatidia/config.py
# interp. every user-settable parameter of a run. Serialised to config.json next
#         to every result and stored as JSON on the generated lattice object.
#         The UI PropertyGroup is a view of it (converted only in config_props.py).
@dataclass
class Config:
    # density
    density_mode: str = "COUNT"          # DENSITY_MODE
    spacing: float = 16e-6               # SPACING (m)
    target_count: int = 750              # TARGET_COUNT
    max_count: int = 50_000              # MAX_COUNT
    model_scale: float = 1.0             # MODEL_SCALE (BU per specimen metre)
    # lattice
    seed: int = 0                        # SEED (only the seams are random) [M3]
    lattice_origin: tuple | None = None  # LATTICE_ORIGIN [M3]
    lattice_rotation_deg: float = 0.0    # LATTICE_ROTATION [M3]
    relax_iterations: int = 200          # RELAX_ITERATIONS (seams only) [M3]
    irregularity: float = 0.0            # IRREGULARITY (× spacing)
    # facet shape
    facet_shape: float = 0.0             # FACET_SHAPE
    inflation: float = 0.4               # INFLATION (× facet inner radius) [M3]
    pokeyness: float = 0.0               # POKEYNESS (0 round, 1 cone) [M3]
    crease_sharpness: float = 0.5        # CREASE_SHARPNESS
    rim_width: float = 1.5               # RIM_WIDTH (× spacing)
    subdivision_level: int = -1          # facet resolution 0..6; -1 = auto (2) [M3]
    # bristles
    bristle_pattern: str = "ALTERNATING" # BRISTLE_PATTERN
    bristle_frequency: float = 0.5       # BRISTLE_FREQUENCY
    bristle_length: float = 10e-6
    bristle_radius: float = 0.4e-6
    bristle_tilt_deg: float = 20.0       # [M3] was 35; SEM bristles fan out near the normal
    comb_direction: tuple = (0.0, 0.0, -1.0)
    socket_radius: float = 1.0e-6
    socket_depth: float = 0.3e-6
    # corneal nipples
    nipple_mode: str = "SHADER"          # NIPPLE_MODE
    nipple_spacing: float = 0.2e-6
    nipple_height: float = 0.05e-6

config0 = Config() #stub: Drosophila defaults
config1 = Config(density_mode="SPACING", spacing=25e-6, facet_shape=1.0,
                 inflation=0.1, bristle_pattern="NONE")  # wasp-like
```

`Config` gets a `validate()` method that checks every invariant above and a
`to_json()` / `from_json(path)` pair. `DENSITY_MODE = "COUNT"` is the default
because it is scale-free: it works on any mesh without first setting
`MODEL_SCALE`.

---

## 8. Blender-side data (used by `blender_io/`, `geometry_nodes/`, `ui/`)

```
# CDD. LATTICE_OBJECT                                         [PROPOSAL]
# A single Blender mesh object "<target>_OmmatidiaLattice", world space, model units:
#   vertices  = LATTICE.centers  ∪  TOPOLOGY.cell_points
#   edges     = TOPOLOGY.cell_edges (offset into the cell-point block)
#   point attributes:
#     "omm_kind"    INT     0 = ommatidium centre, 1 = cell vertex
#     "omm_normal"  FLOAT_VECTOR
#     "omm_bristle" BOOLEAN True on cell vertices that carry a bristle
#     "omm_valence" INT     valence for centres (debug view of defects), 0 otherwise
#   custom properties:
#     "omm_config"  str     CONFIG as JSON
#     "omm_metrics" str     LATTICE_METRICS as JSON
#     "omm_target"  str     name of the covered object  [M2]
# interp. the hand-off from the one-off Python computation to the live GN tree.
#         GN separates the parts by "omm_kind". One object keeps the scene tidy
#         and gives one place to store what made it.

# CDD. OMMATIDIA_MODIFIER   [REPLACED in M3 by FACETS_OBJECT; proposal §11]

# CDD. FACETS_OBJECT  [M3]
# A mesh object "<target>_Ommatidia", parented to the target (world positions kept):
#   vertices = shape_positions(FACET_MESH, CONFIG) × MODEL_SCALE; faces from FACET_MESH
#   material slots: 0 = facets, 1 = bristles (face material_index = face kind)
#   point attributes = the per-vertex FACET_MESH fields ("omm_base", "omm_normal",
#     "omm_kind", "omm_rho_round", "omm_rho_hex", "omm_radius", "omm_rim",
#     "omm_axial", "omm_angle"; specimen metres), so shape settings are applied
#     without the surface
#   custom properties: "omm_spacing" (specimen metres), "omm_target"
# The LATTICE_OBJECT is now hidden and parented to the target too.
```

```
# DD. CELL_LOOPS  [M3]   facets/cell_loops.py
# cell_loops = list of (m,) int64, one per centre: the corner triangles around it,
#   counter-clockwise; empty where the loop does not close (rim cells, holes)

# DD. FACET_RESOLUTION  [M3]   facets/dome_mesh.py
# (sides, rings) = (L + 1, 2·(L + 1)) for SUBDIVISION_LEVEL L (-1 -> 2)

# CDD. FACET_LAYOUT  [M3]   facets/dome_mesh.py
# points, owner, rho_hex, rho_round, loop_start, loop_vertices: the unprojected
#   composite. Corners and wall points are shared by neighbouring facets.
#   rho_hex = fraction from centre to wall along the ray; rho_round = distance
#   from the centre / the facet's inner radius (mean centre-to-wall-midpoint)

# CDD. FACET_MESH  [M3]   facets/facets.py
# per vertex: base, normal, kind (facet/bristle), rho_round, rho_hex, radius, rim
#   (centre's boundary distance / SPACING), axial, angle; faces: loop_start,
#   loop_vertices, face_kind; spacing
# Invariant: shape_positions(FACET_MESH, CONFIG) is defined for every valid CONFIG

# CDD. BRISTLE_TEMPLATE  [M3]   facets/bristles.py
# (axial (K,), angle (K,), faces): one bristle in local coordinates, base 0, tip 1
```

---

## 9. Proposals that differ from or extend AGENTS.md (please review)

| # | Proposal | Why |
|---|---|---|
| P1 | Units: pure code in specimen SI metres; one `model_scale` applied only in `blender_io/`; facet-tied lengths as multiples of SPACING | Sculpts are rarely real-size; one unit boundary is easy to test |
| P2 | `SURFACE.area()` added to the protocol | Needed to convert spacing ↔ count |
| P3 | New atomic `lattice/.../neighbors.py` (numpy spatial-hash grid) and DD `NEIGHBOR_PAIRS` | No SciPy KD-tree, and `mathutils.kdtree` is off-limits in pure code |
| P4 | `HEX_FRACTION` / defects computed in `topology/` (it needs the triangulation); `lattice/metrics.py` keeps spacing CV and nearest-neighbour stats | Keeps the dependency arrow lattice → topology one-way |
| P5 | Open boundaries: centres closer than 0.5·spacing to the rim are dropped after relaxation; metrics use interior centres only; the visual fade is done live in GN from the input mesh's boundary edges | Repulsion piles points onto the rim; trimming is simple and testable. A mirror-force wall is the fallback if the rim looks ragged |
| P6 | Triangulation = local tangent-plane Delaunay per centre, keeping triangles all three vertices agree on | No SciPy; robust on a relaxed point set (a triangular lattice has no co-circular quadruples) |
| P7 | Both REPULSION and LLOYD implemented behind `relax_method`; benchmark in M1, then pick the default | Dossier §8 open question |
| P8 | One LATTICE_OBJECT with an `omm_kind` attribute, not three objects; GN modifier lives on the user's target object | Tidy scene; one place for CONFIG JSON; non-destructive |
| P9 | `DENSITY_MODE` default = COUNT (800) | Works on a mesh of any scale without setting `model_scale` first |

---

## 10. Changes during M1 and M2 (atomic additions only; no approved axiom was replaced)

| Change | Why |
|---|---|
| `VELOCITIES` + momentum in `repulsion_step` | Plain repulsion left ~27% non-hexagonal cells after 60 steps; momentum 0.9 reaches 90% in ~100 steps |
| `SEAM_CUTS` + `seam_cuts`, `crossed_edges`, `open_edges`, `rim_cells` | A plain BFS colouring broke the bristle alternation on 8% of facet walls; pairing defects by short seams brings it to ≈ one break per defect (3–4%) |
| Triangulation empty-circle test in the triangle's own plane (not each centre's tangent plane) | Tangent-plane tests disagreed near ties and left 4–16 holes; the triangle-plane test is watertight on all fixtures |
| Defaults: `target_count` 750, `relax_iterations` 200 | Atlas value; convergence benchmark (logbook) |
| LLOYD and `CONFIG.relax_method` removed (Uriel, 2026-11-03) | Benchmark: worse on every measure |
| Open sheets: `target_count` includes the trimmed rim row (~12% fewer centres) | Accepted by Uriel (662 vs ~700 is fine) |
| Hexagon-fraction acceptance = 0.88 at 750 centres | Accepted by Uriel: a limit of the problem (5-7 scars on curved surfaces), not of the method |
| [M2] `BvhSurface` normals = interpolated *corner* (split) normals, not vertex normals | Vertex normals on a flat-shaded low-poly mesh point diagonally across whole faces; corner normals follow the user's shading |
| [M2] Triangulation tests each triangle in the plane ⟂ the mean normal of its corners; neighbours with normals facing away are ignored | The triangle-plane test produced 288 overlapping edges on a cube's 90° creases (now 54-68); smooth meshes unchanged (0) |
| [M2] Metric `N_OVERLAP_EDGES`; UI warns on broken facet edges | Hard creases remain a known limitation (dossier §8) |
| [M2] DD. MESH_ARRAYS, CDD. OMMATIDIA_SETTINGS (see blender_io.md); lattice object property `omm_target` | Blender side of the hand-off |

---

## 11. [M3] Changes after Uriel's feedback of 2026-11-03 (please review)

Feedback: (1) the lattice should be regular across the eye; carry the defects to
the edges where they can be hidden; (2) the result should look like the SEM:
a composite of multi-face ommatidia, with an inflation control from flat to
globular; (3) changes should apply at once; (4) bristles were missing.

| # | Change | Why |
|---|---|---|
| Q1 | `lattice/`: a perfect hexagonal patch grows ring by ring from `LATTICE_ORIGIN` (equilateral apex placement, spring steps on the patch edges); new patches start in any big uncovered region (second component); the seams are dart-thrown and relaxed with the patches fixed. Random sampling + global repulsion is gone. New DDs HEX_INDEX, PATCH_EDGES, LATTICE_LINKS, PATCH_IDS; new CONFIG `lattice_origin`, `lattice_rotation_deg` | Feedback 1. 750-sphere: hexagons 0.896 → 0.948, perfect to 62° from the centre; cap 0.884 → 0.960, defects only in the last ~2 rows |
| Q2 | Spacing-CV threshold on the sphere 0.05 → 0.07; hexagon threshold 0.88 → 0.93; count within ±7% of the target (curvature compresses the patch) | Regularity of the pattern traded for a few % of cell-size variation, as asked |
| Q3 | `topology/`: `fill_small_holes` closes 3- and 4-edge holes | Very regular lattices have near-co-circular quadruples whose candidate triangles were rejected in slightly different planes |
| Q4 | New compound `facets/` (pure numpy) builds the composite mesh: every facet a fan of rings, walls and corners shared, plus bristle spindles; `shape_positions` applies the shape settings in ~20–30 ms | Feedback 2 and 3. Replaces the planned `geometry_nodes/` (OMMATIDIA_MODIFIER): Python-side shaping is testable with pytest, version-proof (no node API), and still instant |
| Q5 | `INFLATION` (new) = dome height / facet radius; `POKEYNESS` redefined as the dome outline (round → cone) | Feedback 2: "inflation" from flat to semi-dome. Height and outline are separate controls |
| Q6 | Real lengths (bristles) in COUNT mode scale with the lattice (`real_length_scale`) | Feedback 4: with COUNT and model_scale 1 a 10 µm bristle on a 2 m mesh is invisible |
| Q7 | UI: settings update live (shape at once; lattice after 0.4 s); `live_update`, `origin_source`/`origin_object`/`origin_point` and `status` are UI-only properties | Feedback 3 |
| Q8 | `material/`: an SEM-look material (Layer Weight facing → emission) on both slots | First step of M6 |
| Q9 | Default `bristle_tilt_deg` 35 → 20 | In the reference SEM the bristles fan out radially, i.e. they stand close to the normal |
| Q10 | `blender_io`: a thin solid (≥ 25% of its area has an opposite-facing skin closer than one estimated spacing) keeps only its outer skin (faces pointing away from the area-weighted centroid); metric `DROPPED_INNER_FACES` | Mask Extract "as Solid" meshes: the lattice hopped between skins (holes, stray bristles) and overlapping triangles crashed `facet_layout` (KeyError) |
| Q11 | Guards: growth refuses a node whose normal turns > ~95° from its parents' (`MIN_TURN_COSINE` −0.1); triangulation drops neighbour pairs facing opposite ways (cos < −0.5); `facet_layout` skips a facet whose corner loop lacks a wall | No crash and no mixing of skins even where the detection does not apply (thin parts of a larger mesh) |
