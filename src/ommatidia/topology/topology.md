# module name: topology

**purpose:** from the ommatidium centres, derive the facet outlines (which
facets touch, where their corners and walls are) and choose which corners
carry an interommatidial bristle.

Code: [topology.py](topology.py) · previous module: [lattice.md](../lattice/lattice.md) ·
axioms: [data_definitions.md](../../../docs/agents/data_definitions.md)

## Introduction

[lattice](../lattice/lattice.md) gives us *where* each ommatidium is. To draw
flat hexagonal facets and to place bristles, we also need *how they connect*:
which centres are neighbours, where three facets meet (the facet corners), and
where the walls between facets run.

Connecting each centre to its true neighbours is a **triangulation**. Its
*dual* is the facet pattern: every triangle of three neighbouring centres has
one facet corner (the point equidistant from the three, its circumcentre), and
two triangles sharing a side are joined by one facet wall. On a hexagonal
lattice each centre has six neighbours, so each facet is a hexagon. Centres
with five or seven neighbours are the lattice defects (pentagonal or heptagonal
facets).

Without SciPy we triangulate locally. Each centre looks at the centres within
two spacings and proposes every triangle around itself whose circumcircle
(taken in the triangle's own plane) contains no other centre: the Delaunay
rule. Because the test depends only on the triangle and its surroundings, the
three corners of a triangle always reach the same verdict. We keep the
triangles all three proposed. On closed test surfaces this produces a
watertight triangulation (Euler characteristic 2, no holes). Exact ties (four
centres on one circle) are rejected, so they appear as a counted hole rather
than as overlapping triangles.

The plane used for a triangle's test is perpendicular to the mean surface
normal of its three corners. On smooth surfaces that is almost the triangle's
own plane. Across a sharp crease it is the bisecting plane, which "unfolds"
the crease symmetrically. **Limitation:** hard creases (≳ 45°, e.g. a cube or
cylinder caps) still leave a few holes or overlapping triangles near the
crease, and flat-shaded low-poly meshes leave a handful (≈ 0.2% of edges).
They are counted (`HOLES_ON_CLOSED_SURFACE`, `N_OVERLAP_EDGES`) and the UI
suggests a Subdivision Surface or Bevel modifier. Smooth meshes, including a
subdivided cube, have none.

**Bristles.** In *Drosophila* each ommatidium carries about one bristle, at
every other corner of its hexagon. The triangles of a hexagonal lattice come
in two orientations ("up" and "down"), and neighbours always differ, so taking
all triangles of one orientation gives exactly every other corner. That is
about one bristle per ommatidium, because there are about two triangles per
centre. Around a 5- or 7-sided facet, however, the alternation cannot close
(an odd cycle). We therefore draw short **seams** pairing nearby defects (or
running to an open rim) and alternate everywhere else. That leaves roughly one
irregular spot per defect, as in real eyes.

## Data definitions

| DD | Type | Meaning |
|---|---|---|
| TRIANGLES | (T,3) int | centre triangles, counter-clockwise about the outward normal |
| VALENCE | (N,) int | neighbours per centre = sides of its facet (6 = hexagon) |
| CELL_VERTICES | POINTS, NORMALS (T rows) | facet corners: projected circumcentres |
| CELL_EDGES | (E,2) int | facet walls: pairs of adjacent triangles |
| SEAM_CUTS | (E,) bool | facet walls allowed to break the alternation |
| TRIANGLE_PARITY | (T,) int8 | up/down colouring of triangles |
| BRISTLE_MASK | (T,) bool | facet corners that carry a bristle |
| TOPOLOGY | dataclass | all of the above |

## Atomic algorithms

- **local_delaunay_fans** ([triangulate.py](atomic_algorithm_scripts/triangulate.py)):
  each centre proposes its empty-circumcircle triangles (circumradius ≤ spacing).
- **agree_triangles** (same file): keep triangles proposed by all three corners.
- **fill_small_holes** (same file) [M3]: close holes bounded by 3 or 4 edges. A
  very regular lattice has near-co-circular quadruples whose two candidate
  splits are tested in slightly different planes, so both can be rejected.
- **open_edges / valence** (same file): rim/hole edges; facet side counts.
- **circumcenters / circumradii** ([dual.py](atomic_algorithm_scripts/dual.py)):
  facet corners and the size of each triangle.
- **dual_edges / crossed_edges** (same file): facet walls, and which two
  centres each wall separates.
- **seam_cuts** ([bristle_sites.py](atomic_algorithm_scripts/bristle_sites.py)):
  greedily pair odd-sided facets by shortest paths (approximate minimum T-join).
- **triangle_parity / count_conflicts** (same file): 2-colouring around the
  seams, then local repair sweeps.
- **select_bristles** (same file): ALTERNATING, RANDOM or NONE.

## Compound algorithm: `build_topology(lattice, surface, config) -> (TOPOLOGY, metrics)`

```
LATTICE ──neighbor_pairs(2·SPACING)──► local_delaunay_fans ──► agree_triangles ──► fill_small_holes ──► TRIANGLES
TRIANGLES ──valence──► VALENCE
TRIANGLES ──circumcenters──► SURFACE.closest_point ──► CELL_VERTICES
TRIANGLES ──dual_edges──► CELL_EDGES
TRIANGLES ──open_edges──► rim cells ┐
VALENCE (odd, not rim) ─────────────┴─► seam_cuts ──► triangle_parity ──► select_bristles
──► TOPOLOGY, {HEX_FRACTION, N_PENTAGONS, N_HEPTAGONS, N_OTHER_DEFECTS, N_OPEN_EDGES,
               HOLES_ON_CLOSED_SURFACE, MAX_HOLE, PARITY_CONFLICTS, N_BRISTLES, ...}
```

Interior for these metrics = more than one SPACING from the surface's open edge
**and** not touching an open edge of the triangulation.

Measured (750 centres, defaults): hexagon fraction 0.90 (sphere), 0.91
(ellipsoid), 0.88 (open cap); no holes on closed surfaces; one bristle per
ommatidium; colour conflicts ≈ number of defects (3–4% of facet walls). At
5000 centres the hexagon fraction is ≈ 0.92.
