"""Compound: facet outlines and bristle sites from the ommatidium centres.

Pure numpy: never imports bpy or mathutils. See topology.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import Config
from ..lattice.atomic_algorithm_scripts.neighbors import neighbor_pairs
from ..lattice.lattice import Lattice, stage_rngs
from .atomic_algorithm_scripts.bristle_sites import (
    count_conflicts,
    seam_cuts,
    select_bristles,
    triangle_parity,
)
from .atomic_algorithm_scripts.dual import circumcenters, circumradii, crossed_edges, dual_edges
from .atomic_algorithm_scripts.triangulate import (
    agree_triangles,
    fill_small_holes,
    local_delaunay_fans,
    open_edges,
    overused_edges,
    valence,
)

# A triangle whose circumradius exceeds this (x spacing) is a hole, not a cell
# corner; a perfect lattice has 1/sqrt(3) ≈ 0.577. The neighbour radius must be
# at least twice it so every corner sees every point near the circle.
MAX_CIRCUMRADIUS = 1.0
TRIANGULATION_RADIUS = 2.0 * MAX_CIRCUMRADIUS
OPPOSITE_NORMALS = -0.5   # cosine; below it two centres are on different skins


# CDD. TOPOLOGY
# topology = Topology(triangles: TRIANGLES, valence: VALENCE,
#                     cell_points: POINTS, cell_normals: NORMALS,
#                     cell_edges: CELL_EDGES, parity: TRIANGLE_PARITY,
#                     bristle_mask: BRISTLE_MASK)
# interp. everything derived from the centres' connectivity: facet outlines
#         (cell_points + cell_edges) and bristle sites (cell_points[bristle_mask])
# Invariant: cell_points[t] is the projected circumcentre of triangles[t]
@dataclass
class Topology:
    triangles: np.ndarray
    valence: np.ndarray
    cell_points: np.ndarray
    cell_normals: np.ndarray
    cell_edges: np.ndarray
    parity: np.ndarray
    bristle_mask: np.ndarray


# FD. build_topology
# Signature: LATTICE, SURFACE, CONFIG -> TOPOLOGY, LATTICE_METRICS (partial)
def build_topology(lattice: Lattice, surface, config: Config) -> tuple[Topology, dict]:
    """Triangulate the centres, build the facet outlines, and place bristles.

    Algorithm:
        neighbour pairs within TRIANGULATION_RADIUS · spacing, minus opposite-facing ones
        local Delaunay fans -> keep the triangles all three corners agree on
        close holes of 3 or 4 edges (near-co-circular quadruples)
        valence of each centre (6 = hexagonal facet)
        circumcentres projected onto SURFACE -> cell vertices (facet corners)
        triangles sharing an edge -> cell edges (facet walls)
        short seams pairing the 5- and 7-sided cells (or running to the rim)
        2-colour the triangles around the seams -> bristle sites by BRISTLE_PATTERN
        metrics on interior centres
    """
    spacing = lattice.spacing
    pairs = neighbor_pairs(lattice.centers, TRIANGULATION_RADIUS * spacing)
    # Centres facing opposite ways (two skins of a thin part) are never neighbours.
    n = lattice.normals
    pairs = pairs[np.sum(n[pairs[:, 0]] * n[pairs[:, 1]], axis=1) > OPPOSITE_NORMALS]
    fans = local_delaunay_fans(lattice.centers, lattice.normals, pairs, MAX_CIRCUMRADIUS * spacing)
    triangles = fill_small_holes(agree_triangles(fans), lattice.centers, MAX_CIRCUMRADIUS * spacing)
    val = valence(triangles, len(lattice))
    centres_3d = circumcenters(lattice.centers, triangles)
    radii = circumradii(lattice.centers, triangles, centres_3d)
    cell_points, cell_normals = surface.closest_point(centres_3d)
    if len(triangles) == 0:
        cell_points, cell_normals = np.zeros((0, 3)), np.zeros((0, 3))
    cell_edges = dual_edges(triangles)
    rim = rim_cells(triangles, len(lattice))
    cuts = seam_cuts(crossed_edges(triangles), ~rim & (val % 2 == 1), rim)
    parity = triangle_parity(len(triangles), cell_edges, cuts)
    bristles = select_bristles(parity, config.bristle_pattern, config.bristle_frequency,
                               stage_rngs(config.seed)["bristles"])
    topology = Topology(triangles, val, cell_points, cell_normals, cell_edges, parity, bristles)
    return topology, topology_metrics(topology, lattice, surface, radii)


# FD. rim_cells
# Signature: TRIANGLES, int -> (N,) bool
def rim_cells(triangles: np.ndarray, n: int) -> np.ndarray:
    """Centres touching an open edge of the triangulation (their cells are cut open)."""
    mask = np.zeros(n, dtype=bool)
    mask[open_edges(triangles).ravel()] = True
    return mask


# FD. topology_metrics
# Signature: TOPOLOGY, LATTICE, SURFACE, DISTANCES (circumradii) -> LATTICE_METRICS (partial)
def topology_metrics(topology: Topology, lattice: Lattice, surface, radii: np.ndarray) -> dict:
    """Hexagon fraction and defect counts on interior centres, plus hole and bristle counts.

    Interior here = more than one SPACING from an open edge of the SURFACE and
    not touching an open edge of the triangulation (whose cells are cut open).
    """
    on_rim = rim_cells(topology.triangles, len(lattice))
    boundary_distance = surface.boundary_distance(lattice.centers)
    interior = ~on_rim & (boundary_distance > lattice.spacing)
    v = topology.valence[interior]
    n_open = len(open_edges(topology.triangles))
    closed = bool(np.all(np.isinf(boundary_distance)))
    return {
        "N_TRIANGLES": len(topology.triangles),
        "N_TOPO_INTERIOR": int(interior.sum()),
        "HEX_FRACTION": float(np.mean(v == 6)) if len(v) else float("nan"),
        "N_PENTAGONS": int(np.sum(v == 5)),
        "N_HEPTAGONS": int(np.sum(v == 7)),
        "N_OTHER_DEFECTS": int(np.sum((v < 5) | (v > 7))),
        "N_OPEN_EDGES": n_open,                      # holes, if the surface is closed
        "HOLES_ON_CLOSED_SURFACE": n_open if closed else 0,
        "N_OVERLAP_EDGES": len(overused_edges(topology.triangles)),  # sharp creases only
        "MAX_HOLE": float(radii.max() / lattice.spacing) if len(radii) else float("nan"),
        "PARITY_CONFLICTS": count_conflicts(topology.parity, topology.cell_edges),
        "N_BRISTLES": int(topology.bristle_mask.sum()),
    }
