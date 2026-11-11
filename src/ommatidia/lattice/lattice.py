"""Compound: place a regular lattice of ommatidium centres on a SURFACE.

Pure numpy: never imports bpy or mathutils. See lattice.md.

Regularity first: a perfect hexagonal patch grows from a centre until curvature
or the rim stops it; the leftover gaps (seams) are filled with relaxed random
points. Defects therefore collect in the seams, away from the centre.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from ..config import Config
from .atomic_algorithm_scripts.density import HEX_CELL_FACTOR, resolve_density
from .atomic_algorithm_scripts.growth import (
    DEAD,
    HEX_DIRECTIONS,
    RIM_DISTANCE,
    apex_points,
    apex_triples,
    auto_origin,
    first_ring,
    mean_by_owner,
    node_links,
    patch_edges,
    placement_mask,
    ring_candidates,
    spring_step,
    tangent_frame,
)
from .atomic_algorithm_scripts.metrics import nearest_neighbor_distances, spacing_cv
from .atomic_algorithm_scripts.neighbors import nearest_distances, neighbor_pairs
from .atomic_algorithm_scripts.relax import (
    REPULSION_RADIUS,
    jitter_tangent,
    repulsion_step,
    trim_boundary,
)
from .atomic_algorithm_scripts.sampling import fill_mask, sample_candidates

SAMPLES_PER_CENTER = 30          # surface samples per ommatidium, for coverage and seam filling
NN_SEARCH_RADIUS = 1.5           # x spacing; enough to find the nearest neighbour
# A sample farther than this from every centre is uncovered (a perfect lattice
# has at most spacing/sqrt(3) = 0.577; strain stretches that a little).
COVER_DISTANCE = 0.75            # x spacing
SEED_GAP = 2.0                   # x spacing; a gap this deep gets its own patch (e.g. a second eye)
MAX_PATCHES = 8
FILL_RADIUS = 0.85               # x spacing; dart-throwing radius in the seams
SPRING_ITERATIONS = 10           # spring steps after each ring
SPRING_BAND = 4                  # rings behind the front that the springs may move
SEAM_RIM_MARGIN = 1.5            # x spacing; patch edges closer to the rim than this stay put
CROWDED = 0.75                   # x spacing; closer centres join the seam relaxation

# Order of the RNG streams spawned from SEED. Append new stages; never reorder.
RNG_STAGES = ("sampling", "relax", "jitter", "bristles")
SEAM = -1                        # patch id of a seam-fill centre


# CDD. LATTICE
# lattice = Lattice(centers: POINTS, normals: NORMALS, spacing: SPACING, patch: PATCH_IDS)
# interp. the ommatidium centres on a SURFACE: perfect hexagonal patches plus
#         seam-fill centres, after boundary trimming and jitter
# PATCH_IDS: (N,) int; patch number of each centre, SEAM (-1) for seam fill
# Invariant: centres lie on the surface; no two centres closer than 0.5·spacing;
#            on open surfaces every centre has boundary_distance >= 0.5·spacing
@dataclass
class Lattice:
    centers: np.ndarray
    normals: np.ndarray
    spacing: float
    patch: np.ndarray = field(default=None)

    def __post_init__(self):
        if self.patch is None:
            self.patch = np.full(len(self.centers), SEAM, dtype=np.int64)

    def __len__(self) -> int:
        return len(self.centers)


lattice0 = Lattice(np.zeros((0, 3)), np.zeros((0, 3)), 1.0)  # stub: empty


# FD. stage_rngs
# Signature: SEED -> {stage name: RNG}
def stage_rngs(seed: int) -> dict[str, np.random.Generator]:
    """One independent generator per pipeline stage, so randomness in one stage
    never shifts another (DD. RNG)."""
    children = np.random.SeedSequence(seed).spawn(len(RNG_STAGES))
    return {name: np.random.default_rng(s) for name, s in zip(RNG_STAGES, children)}


# FD. build_lattice
# Signature: SURFACE, CONFIG -> LATTICE, LATTICE_METRICS (partial)
def build_lattice(surface, config: Config) -> tuple[Lattice, dict]:
    """Cover SURFACE with a regular, unstretched lattice of centres.

    Algorithm:
        resolve SPACING and N from the surface area and CONFIG
        draw area-weighted surface samples (for coverage tests)
        origin = CONFIG.lattice_origin, or the automatic centre
        repeat (at most MAX_PATCHES):
            grow a perfect patch from origin (grow_patch)
            uncovered = samples farther than COVER_DISTANCE from every centre
            stop unless some uncovered sample is SEED_GAP deep; else origin = deepest one
        seam fill: dart-throw uncovered samples (FILL_RADIUS)
        relax the seams (relax_seams): seam-fill centres, and patch centres with a
            missing neighbour or a too-close one (away from the rim) move; the rest
            stays fixed
        drop centres closer than 0.5·SPACING to an open edge
        if IRREGULARITY > 0: tangential jitter, project
        measure nearest-neighbour spacing on interior centres
    """
    config.validate()
    rngs = stage_rngs(config.seed)
    seconds = {}

    t0 = time.perf_counter()
    spacing, n = resolve_density(surface.area(), config.density_mode, config.spacing,
                                 config.target_count, config.max_count)
    samples, sample_normals = surface.closest_point(
        sample_candidates(surface, SAMPLES_PER_CENTER * n, rngs["sampling"]))
    sample_rim = surface.boundary_distance(samples)
    fillable = sample_rim >= RIM_DISTANCE * spacing
    if config.lattice_origin is not None:
        origin = np.asarray(config.lattice_origin, dtype=float)
    else:
        origin = auto_origin(samples, sample_rim)
    seconds["sampling"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    points, normals = np.zeros((0, 3)), np.zeros((0, 3))
    patch = np.zeros(0, dtype=np.int64)
    open_side = np.zeros(0, dtype=bool)
    for patch_id in range(MAX_PATCHES):
        grown, grown_normals, grown_open = grow_patch(surface, origin, config.lattice_rotation_deg,
                                                      spacing, points, normals)
        points = np.concatenate([points, grown])
        normals = np.concatenate([normals, grown_normals])
        patch = np.concatenate([patch, np.full(len(grown), patch_id, dtype=np.int64)])
        open_side = np.concatenate([open_side, grown_open])
        gap = nearest_distances(samples, points, 2 * SEED_GAP * spacing)
        deep = np.where(fillable, gap, 0.0)
        if deep.max(initial=0.0) < SEED_GAP * spacing:
            break
        origin = samples[int(np.argmax(deep))]
    seconds["growth"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    uncovered = fillable & (nearest_distances(samples, points, 2 * spacing) > COVER_DISTANCE * spacing)
    candidates, candidate_normals = samples[uncovered], sample_normals[uncovered]
    order = rngs["relax"].permutation(len(candidates))
    candidates, candidate_normals = candidates[order], candidate_normals[order]
    keep = fill_mask(candidates, points, candidate_normals, normals, FILL_RADIUS * spacing)
    points = np.concatenate([points, candidates[keep]])
    normals = np.concatenate([normals, candidate_normals[keep]])
    patch = np.concatenate([patch, np.full(int(keep.sum()), SEAM, dtype=np.int64)])
    open_side = np.concatenate([open_side, np.ones(int(keep.sum()), dtype=bool)])
    # Patch edges facing a seam move too (not those facing the rim: they would
    # be pushed off it), and so does any centre crowded by another front.
    away_from_rim = surface.boundary_distance(points) > SEAM_RIM_MARGIN * spacing
    crowded = nearest_neighbor_distances(points, neighbor_pairs(points, CROWDED * spacing)) < CROWDED * spacing
    movable = (patch == SEAM) | ((open_side | crowded) & away_from_rim)
    points, normals = relax_seams(surface, points, normals, movable, spacing,
                                  config.relax_iterations)
    seconds["seams"] = time.perf_counter() - t0

    keep = trim_boundary(surface.boundary_distance(points), spacing)
    points, normals, patch = points[keep], normals[keep], patch[keep]

    if config.irregularity > 0 and len(points):
        points = jitter_tangent(points, normals, config.irregularity, spacing, rngs["jitter"])
        points, normals = surface.closest_point(points)

    lattice = Lattice(points, normals, spacing, patch)
    metrics = lattice_metrics(lattice, surface)
    metrics["N_TARGET"] = n
    metrics["SECONDS"] = seconds
    return lattice, metrics


# FD. grow_patch
# Signature: SURFACE, POINT, float (degrees), SPACING, POINTS, NORMALS -> POINTS, NORMALS, (K,) bool
def grow_patch(surface, origin: np.ndarray, rotation_deg: float, spacing: float,
               existing: np.ndarray, existing_normals: np.ndarray
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Grow one perfect hexagonal patch from origin, avoiding existing points.
    Also returns which nodes have fewer than 6 lattice neighbours (the patch edge).

    Algorithm:
        centre = origin projected onto SURFACE; stop if it may not be placed
        ring 1: six points at SPACING around the centre (rotated by rotation_deg)
        for ring k = 2, 3, ... while the previous ring placed anything:
            candidates = unvisited ring-k neighbours of ring k−1
            repeat until no candidate gains a proposal:
                proposals = apex of equilateral triangles on placed neighbour pairs
                average per node (and its parents' normals); project;
                placement_mask -> place or mark dead
            spring steps on the patch edges, moving only the last SPRING_BAND rings
        patch edge = nodes with fewer than 6 placed lattice neighbours
    """
    center, center_normal = surface.closest_point(np.asarray(origin, dtype=float)[None])
    ok = placement_mask(center, center, center_normal, surface.boundary_distance(center),
                        existing, existing_normals, spacing)
    if not ok[0]:
        return np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0, dtype=bool)
    e1, e2 = tangent_frame(center_normal[0], rotation_deg)
    points, normals = center.copy(), center_normal.copy()
    index = {(0, 0): 0}
    ring_of = [0]

    def place(nodes, proposals, parent_normals, ring):
        nonlocal points, normals
        x, nx = surface.closest_point(proposals)
        links, batch_links = node_links(nodes, index, len(existing))
        alive = placement_mask(x, proposals, nx, surface.boundary_distance(x),
                               np.concatenate([existing, points]),
                               np.concatenate([existing_normals, normals]), spacing,
                               links, batch_links, parent_normals)
        rows = len(points) + np.cumsum(alive) - 1
        for node, a, row in zip(nodes, alive, rows):
            index[node] = int(row) if a else DEAD
        points = np.concatenate([points, x[alive]])
        normals = np.concatenate([normals, nx[alive]])
        ring_of.extend([ring] * int(alive.sum()))
        return [node for node, a in zip(nodes, alive) if a]

    front = place([tuple(d) for d in HEX_DIRECTIONS], first_ring(center[0], e1, e2, spacing),
                  np.repeat(center_normal, 6, axis=0), 1)
    ring = 1
    while front:
        ring += 1
        pending = ring_candidates(front, index, ring)
        front = []
        while pending:
            owner, a, b, c = apex_triples(pending, index)
            if len(owner) == 0:
                break
            means, has = mean_by_owner(apex_points(points, normals, a, b, c, spacing), owner, len(pending))
            parents, _ = mean_by_owner(normals[a] + normals[b], owner, len(pending))
            front += place([pending[k] for k in has], means, parents, ring)
            pending = [node for node in pending if node not in index]
        if front:
            points, normals = _relax_front(surface, points, normals, index, np.asarray(ring_of),
                                           ring, spacing)
    degree = np.bincount(patch_edges(index).ravel(), minlength=len(points))
    return points, normals, degree < 6


def _relax_front(surface, points, normals, index, ring_of, ring, spacing):
    """Spring steps on the last SPRING_BAND rings (one ring more as fixed anchors);
    a step that would cross the rim is undone for that point."""
    sub = np.flatnonzero(ring_of >= ring - SPRING_BAND - 1)
    local = np.full(len(points), -1)
    local[sub] = np.arange(len(sub))
    edges = local[patch_edges(index)]
    edges = edges[(edges >= 0).all(axis=1)]
    movable = ring_of[sub] >= ring - SPRING_BAND
    x, nx = points[sub], normals[sub]
    velocities = np.zeros_like(x)
    for _ in range(SPRING_ITERATIONS):
        moved, velocities = spring_step(x, velocities, nx, edges, movable, spacing)
        px, pn = surface.closest_point(moved[movable])
        bad = surface.boundary_distance(px) < RIM_DISTANCE * spacing
        px[bad], pn[bad] = x[movable][bad], nx[movable][bad]
        v = velocities[movable]
        v[bad] = 0.0
        velocities[movable] = v
        x[movable], nx[movable] = px, pn
    points, normals = points.copy(), normals.copy()
    points[sub], normals[sub] = x, nx
    return points, normals


# FD. relax_seams
# Signature: SURFACE, POINTS, NORMALS, (N,) bool (movable), SPACING, RELAX_ITERATIONS -> POINTS, NORMALS
def relax_seams(surface, points: np.ndarray, normals: np.ndarray, seam: np.ndarray,
                spacing: float, iterations: int) -> tuple[np.ndarray, np.ndarray]:
    """Repulsion relaxation of the seam centres only; every other centre stays fixed.

    Algorithm:
        local set = movable centres + fixed centres within 2·REPULSION_RADIUS of one
        repeat iterations times:
            neighbour pairs -> repulsion step -> undo the move of fixed centres
            -> project the movable centres onto SURFACE
    """
    if not seam.any() or iterations == 0:
        return points, normals
    near = nearest_distances(points, points[seam], 2 * REPULSION_RADIUS * spacing)
    sub = np.flatnonzero(seam | (near < 2 * REPULSION_RADIUS * spacing))
    movable = seam[sub]
    x, nx = points[sub], normals[sub]
    velocities = np.zeros_like(x)
    for _ in range(iterations):
        pairs = neighbor_pairs(x, REPULSION_RADIUS * spacing)
        moved, velocities = repulsion_step(x, velocities, nx, pairs, spacing)
        velocities[~movable] = 0.0
        x[movable], nx[movable] = surface.closest_point(moved[movable])
    points, normals = points.copy(), normals.copy()
    points[sub], normals[sub] = x, nx
    return points, normals


# FD. interior_mask
# Signature: LATTICE, SURFACE -> (N,) bool
def interior_mask(lattice: Lattice, surface) -> np.ndarray:
    """Centres more than one SPACING from an open edge (all of them on a closed surface)."""
    return surface.boundary_distance(lattice.centers) > lattice.spacing


# FD. lattice_metrics
# Signature: LATTICE, SURFACE -> LATTICE_METRICS (partial)
def lattice_metrics(lattice: Lattice, surface) -> dict:
    """Spacing statistics on interior centres, and how much of the lattice is patch.

    Algorithm:
        neighbour pairs within 1.5·SPACING -> nearest-neighbour distances
        mean and coefficient of variation over interior centres
        count patches and seam centres
    """
    pairs = neighbor_pairs(lattice.centers, NN_SEARCH_RADIUS * lattice.spacing)
    nn = nearest_neighbor_distances(lattice.centers, pairs)
    interior = interior_mask(lattice, surface)
    finite = interior & np.isfinite(nn)
    in_patch = lattice.patch != SEAM
    return {
        "N": len(lattice),
        "N_INTERIOR": int(interior.sum()),
        "SPACING": lattice.spacing,
        "MEAN_NN": float(nn[finite].mean()) if finite.any() else float("nan"),
        "MIN_NN": float(nn.min()) if len(nn) else float("nan"),
        "SPACING_CV": spacing_cv(nn, interior),
        "N_PATCHES": int(len(np.unique(lattice.patch[in_patch]))),
        "N_SEAM": int((~in_patch).sum()),
    }


__all__ = ["Lattice", "build_lattice", "grow_patch", "relax_seams", "interior_mask",
           "lattice_metrics", "stage_rngs", "SEAM", "HEX_CELL_FACTOR"]
