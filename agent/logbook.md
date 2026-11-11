# Logbook

## 2026-11-03 — M0: axioms drafted

- Wrote `docs/agents/data_definitions.md` and `docs/agents/modules.md` from the
  dossier and AGENTS.md. Waiting for Uriel's review; no implementation code yet.
- Environment: Blender 4.5 and 5.1 installed. `.venv` is Python 3.12.10 with
  no numpy or pytest yet (needed for M1).
- Main proposals (P1–P9 in data_definitions.md §9):
  - Units: pure modules in specimen metres; `model_scale` converts in
    `blender_io/` only. Facet-tied lengths are multiples of spacing.
  - `SURFACE.area()` added; numpy spatial hash replaces a KD-tree.
  - Hexagon fraction moves to `topology/`.
  - Open boundaries: trim centres within 0.5·spacing of the rim; fade in GN.
  - Triangulation by local tangent-plane Delaunay (no SciPy).
  - One lattice object with an `omm_kind` attribute; GN modifier on the target.
  - Default density mode = COUNT (scale-free).
- Biology values still needing citations: bristle placement (Cagan & Ready
  1989?), ommatidium count (Ready, Hanson & Benzer 1976?), nipple spacing and
  height, bristle tilt and comb direction.

## 2026-11-03 — M0 approved; M1 (lattice) and M1.5 (topology) done

**Review outcome.** Uriel approved P1–P9 and `model_scale`. Source for counts
and bristles: "An anatomical atlas of *Drosophila melanogaster*—the wild-type":
~700 ommatidia (male) to ~750 (female), "each accompanied by a single
mechanosensory bristle". Default count is now 750. Uriel's note: exact counts
matter less than they seem, since users eyeball density on their own mesh, and
anyone needing precision should work from SEM/µCT data.

**Setup.** `.venv` now has numpy 2.5.4 and pytest 9.2.0. Blender 4.5 bundles
numpy 1.26.4, so code avoids numpy-2-only APIs. `pyproject.toml` puts `src/`
and `test/` on the pytest path. 133 pure tests pass in ~13 s.

**Finding 1: spacing CV is a weak regularity measure.** A lattice can have
CV ≈ 0.03 and still only 73% hexagonal cells (it is "glassy": grains and
dislocations). From now on, the hexagon fraction from the triangulation is the
quality number to watch; CV is secondary.

**Finding 2: relaxation benchmark (sphere, 750 centres).**

| Method | Iter. | Spacing CV | Hexagons | Time |
|---|---|---|---|---|
| Repulsion, step 0.1 | 60 | 0.051 | (not measured) | 0.4 s |
| Repulsion, step 0.6 | 60 | 0.033 | 0.73 | 0.5 s |
| Repulsion, step 0.6 | 300 | 0.031 | 0.88 | 1.4 s |
| Repulsion, step 0.6, momentum 0.9 | 100 | 0.033 | 0.90 | 0.6 s |
| Repulsion, step 0.6, momentum 0.9 | 300 | 0.032 | 0.92 | 1.5 s |
| Lloyd (sampled, 40/centre) | 100 | 0.106 | 0.52 | 12 s |
| Exact local Lloyd (Voronoi cell centroid) | 150 | 0.043 | 0.84 | 7 s |

Repulsion radii below ~1.7·spacing fail (no second-ring interaction, so no
shear stiffness: hexagons ≈ 0.43). Chosen: radius 2.0, step 0.6, momentum 0.9,
200 iterations.

**Finding 3: ~90% hexagons is close to the physical limit at 750 on a sphere.**
Pentagons − heptagons is always exactly 12 (Euler). For N ≳ 400, ground states
on a sphere carry extra 5-7 "scars" (Bausch et al. 2003, Science). Measured:
sphere 0.896, ellipsoid 0.907, cap 0.884. At 5000 centres it is 0.92. Tests
use ≥ 0.88; the dossier's 0.90 is borderline at 750.

**Finding 4: triangulation holes.** The first version tested empty circles in
each centre's tangent plane. Neighbouring centres disagreed near ties, leaving
4–16 open edges on the sphere. Testing in the triangle's own plane makes the
verdict identical for all three corners. Result: 0 holes, Euler characteristic
2, T = 2N − 4 on every closed fixture, even on unrelaxed random points. Exact
ties now reject both diagonals (a counted hole, never overlapping triangles).

**Finding 5: bristle alternation.** A plain BFS 2-colouring left long seams
behind every defect, breaking 8% of facet walls (184 of 2244). Pairing odd
(5/7) facets by greedy shortest paths first brings this to 81 (3.6%), about one
per defect. Bristles per ommatidium ≈ 1.00, matching the atlas.

**Finding 6: open rim.** Repulsion pins ~12% of centres onto the rim at
0.88·spacing. They are trimmed as planned (P5), so `target_count` 750 on the
test cap yields ~662. Open question for Uriel (dossier §8).

**Bug found and fixed:** `assign_nearest` trusted grid candidates farther than
the search radius, so a nearer centre just outside the searched cells could be
missed. Fixed; brute-force comparison test added.

**Performance (pure numpy, analytic surfaces).** 750 centres: lattice ~1 s,
topology ~0.1 s. 5000 centres: lattice 7–12 s, topology ~1 s. The ellipsoid is
slower only because its test projection uses bisection. Real-mesh timing
(BVH) comes in M2.

## 2026-11-03 — Lloyd removed; M2 (Blender side) done

**Decisions (Uriel).** Lloyd relaxation and `CONFIG.relax_method` removed. On
open sheets, `target_count` includes the trimmed rim row (662 vs ~700 is
fine). Hexagon fraction 0.88 accepted as a domain limit. AGENTS.md layout
comment updated ("repulsion relaxation").

**Built.** `blender_io/` (BVH SURFACE, lattice object, settings ⇄ CONFIG),
`ui/` (Generate operator, sidebar panel with diagnostics), registration in
`__init__.py`, `blender_manifest.toml`, `main.py`, `test/run_blender_tests.py`
with 16 headless tests (pass on 4.5.3 and 5.1.1), user guide in `docs/user/`.
The built zip installs as `bl_ext.user_default.ommatidia` in an isolated
profile and generates on the default cube.

**API notes.** `BVHTree.FromPolygons(..., all_triangles=...)` is keyword-only
in 5.1, so it is always passed by name. `mesh.corner_normals` and
`loop_triangles.loops` give per-corner shading normals (4.1+). Blender's
icosphere `subdivisions=1` is the bare icosahedron: level 8 is 163,842 vertices.

**Finding: vertex normals were wrong for flat-shaded meshes.** On the default
cube, interpolated vertex normals point diagonally across each face, which
distorted the tangential relaxation (54% hexagons). Corner normals follow the
user's shading (face normals on flat faces); the cube went to 79%, and smooth
meshes are unchanged.

**Finding: hard creases break the triangulation locally.** Testing in the
triangle's own plane ignored points across a 90° fold (288 overlapping edges
on a cube). The test now uses the plane perpendicular to the mean corner
normal and ignores only neighbours whose normal faces away. Results:

| Mesh (750) | Hexagons | Open edges | Overlaps |
|---|---|---|---|
| smooth icosphere | 0.90 | 0 | 0 |
| subdivided cube | 0.89 | 0 | 0 |
| flat icosahedron (42° creases) | 0.88 | 4 | 0 |
| flat icosphere L3 (11°) | 0.86 | 5 | 2 |
| bevelled cube | 0.82 | 23 | 1 |
| 16-sided cylinder (90° caps) | 0.87 | 112 | 24 |
| raw cube | 0.81–0.83 | 30–36 | 54–68 |

Counted as `HOLES_ON_CLOSED_SURFACE` / `N_OVERLAP_EDGES`; the UI suggests a
Subdivision Surface or Bevel modifier. Open question in dossier §8: is a repair
pass worth it, given that eyes are smooth?

**Performance (AGENTS §8).** 800 ommatidia on a 163,842-vertex ellipsoid:
reading the mesh ~1 s, sampling 0.1 s, relaxation 2.6–3.8 s (≈160k BVH
queries in a Python loop), topology 0.1–0.2 s. Hexagons 0.905, CV 0.032.

## 2026-11-03 — M3: regularity first, domed facets, bristles, live updates

**Feedback (Uriel).** (1) The pattern looked inconsistent; prioritise
regularity and carry the artifacts to the edges where they can be hidden.
(2) The output should look like the SEM: multi-face ommatidia forming a
composite, with an "inflation" control from flat to globular. (3) Changes
should apply without pressing Generate. (4) Bristles were missing.

**Lattice: grow a perfect patch.** New `lattice/growth.py`. From a centre,
each ring is placed at the apex of equilateral triangles on two placed
neighbours; 10 spring steps on the last 4 rings spread the strain; a node dies
if it would overlap a non-neighbour (< 0.7·s), a lattice neighbour (< 0.55·s),
leave the surface or reach the rim. Leftover gaps are dart-thrown and relaxed
by repulsion with the patch fixed (patch-edge and crowded nodes may move).
Prototype sweep (750, unit sphere unless noted):

| Variant | Hex (sphere) | Perfect to | CV | Notes |
|---|---|---|---|---|
| Old: Poisson + 200 repulsion | 0.896 | — | 0.032 | defects everywhere |
| Growth, no springs | 0.831 | 49° | — | compression 0.71·s on the 6 sector lines |
| + free repulsion afterwards | 0.92 | — | — | destroys the patch (hex < 60° drops to 0.94) |
| + springs per ring (band 4) | 0.917 | 49° | 0.058 | min NN 0.19·s where fronts meet |
| + seam relaxation of patch edges | 0.917 | 49° | 0.058 | min NN 0.80·s |
| + neighbour collision 0.55·s | 0.943 | 62° | 0.057 | 5000: min NN 0.07·s (fronts drift into each other) |
| + crowded nodes relax too (chosen) | 0.948 | 62° | 0.055 | 5000: 0.958, 65°, CV 0.068 |

Neighbour threshold 0.4 reached 68° but overcrowded the 5000 sphere (N 1.06×,
CV 0.10). A bigger spring band or more steps did not help: the strain of a
perfect patch is geometric (cells must shrink ~cos²(θ/2) toward angle θ).
Final, 750: sphere 0.948, ellipsoid 0.949, cap 0.960 (perfect up to ~2 rows
from the rim); 800 on the 164k-vertex mesh 0.947 in 2.5 s (was 0.905 in 3–4 s).

**Finding: near-co-circular holes.** A very regular lattice has quadruples
where both splits were rejected (each tested in its own mean-normal plane):
6 open edges on the sphere. `fill_small_holes` closes 3- and 4-edge hole loops.

**Facets (new `facets/`, pure numpy).** Each facet is a fan of rings to its
walls; walls and corners are shared, so the composite is watertight and
consistently oriented (tested). Build ~0.3 s; reshaping 104k vertices
20–30 ms in Blender. Shape = dossier §4.3: round cap (inflation × inner radius;
1 = hemisphere; pokeyness bends it toward a cone) blended with a flat plateau
(facet shape); crease width; rim fade. Decided against Geometry Nodes for this
(same responsiveness, pytest-testable, no node-API drift); flagged for review.

**Bristles were invisible, not missing.** Lengths are real µm; with COUNT mode
on a unit mesh (model_scale 1) a 10 µm bristle is 1e-5 of the eye. In COUNT
mode real lengths now scale so that one spacing = `config.spacing` (16 µm).
Default tilt 35° → 20°: in the reference SEM bristles fan out radially.

**Live updates.** Shape settings reshape in their update callback; lattice
settings regenerate 0.4 s after the last change (bpy.app.timers; timers do not
fire in `-b`, so the test runs the scheduled call directly).

**API notes.** Mesh from arrays: `vertices.add`, `loops.add`,
`polygons.add` + `foreach_set("loop_start")` (loop_total derived) works on 4.5
and 5.1. Blender 5 renamed shader sockets ("Fac" → "Factor"): link by index.
`Material.use_nodes` is deprecated in 5.x; only set it before 5.0.

**Tests.** 193 pure (pytest), 26 Blender (4.5.3 and 5.1.1). Raw cube: 78
holes, 88 overlaps (was 30–36 / 54–68): hard creases remain a limitation.

## 2026-11-03 — M3 fix: eyes made with Mask Extract (thin solids)

**Report (Uriel).** On a mask-extracted eye: sometimes `KeyError` in
`facet_layout`; otherwise an incomplete mesh (holes, bristles everywhere).
Planes and single-layer meshes were fine; Uriel suspected two close layers.

**Cause.** Mask Extract's "Extract as Solid" makes two skins a hair apart.
Reproduced with a cap + Solidify (thickness 0.002–0.1 at spacing 0.115): the
patch hopped between skins (closest point picks either), so the outer skin got
holes and inner-skin bristles poked out; neighbour search joined centres of both
skins into overlapping triangles, whose corner loops then lacked a wall
(KeyError). COUNT also spread the count over both skins (~410 visible of 750).

**Fix.** (1) `blender_io.inner_layer_faces`: ray from every face backwards up
to one estimated spacing; if ≥ 25% of the area is doubled, drop the doubled
faces facing toward the centroid. Solid caps now match the single-layer cap
(745 centres, 96% hexagons, 1 patch) at every thickness. (2) Guards, which
alone already prevent the crash: growth refuses nodes turning > ~95° from their
parents (75° dropped the raw cube from 0.854 to 0.720 hexagons, so −0.1 is the
cosine); triangulation ignores opposite-facing neighbour pairs; facets whose
loop lacks a wall are skipped.

**Cost.** One BVH ray per face: reading the 164k-vertex mesh went from ~1 s to
~2–2.6 s. Tests: 197 pure, 28 Blender (4.5.3 and 5.1.1).
