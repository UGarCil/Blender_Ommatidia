"""topology/: triangulation, dual (facet outlines), bristle sites, and the compound."""

import math

import numpy as np
import pytest

from fixtures.analytic_surfaces import CapSurface, EllipsoidSurface, SphereSurface
from ommatidia.config import Config
from ommatidia.lattice.atomic_algorithm_scripts.neighbors import neighbor_pairs
from ommatidia.lattice.lattice import Lattice, build_lattice
from ommatidia.topology.atomic_algorithm_scripts.bristle_sites import (
    count_conflicts,
    seam_cuts,
    select_bristles,
    triangle_parity,
)
from ommatidia.topology.atomic_algorithm_scripts.dual import circumcenters, circumradii, dual_edges
from ommatidia.topology.atomic_algorithm_scripts.triangulate import (
    agree_triangles,
    fill_small_holes,
    local_delaunay_fans,
    open_edges,
    overused_edges,
    valence,
)
from ommatidia.topology.topology import build_topology


def rhombus_patch(n=12, a=1.0):
    """A perfect triangular lattice on an n x n rhombus (straight edges, no ragged rim)."""
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    x = (i + 0.5 * j) * a
    y = j * a * math.sqrt(3) / 2
    return np.column_stack([x.ravel(), y.ravel(), np.zeros(n * n)])


def triangulate(points, normals, spacing=1.0):
    pairs = neighbor_pairs(points, 2.0 * spacing)
    return agree_triangles(local_delaunay_fans(points, normals, pairs, spacing))


def edge_use_counts(tri):
    e = np.sort(tri[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1)
    return np.unique(e, axis=0, return_counts=True)[1]


UP = np.array([0.0, 0.0, 1.0])


# ---------- triangulation on a perfect planar lattice ----------

@pytest.fixture(scope="module")
def hexgrid():
    p = rhombus_patch(12)
    n = np.tile(UP, (len(p), 1))
    return p, n, triangulate(p, n)


def test_perfect_lattice_triangulation(hexgrid):
    p, _, tri = hexgrid
    assert len(tri) == 2 * (12 - 1) * (12 - 1)             # every rhombus split in two
    assert edge_use_counts(tri).max() == 2                   # manifold
    v = valence(tri, len(p))
    rim = np.zeros(len(p), bool)
    rim[open_edges(tri).ravel()] = True
    assert np.all(v[~rim] == 6)


def test_triangles_are_counter_clockwise_about_normal(hexgrid):
    p, _, tri = hexgrid
    nrm = np.cross(p[tri[:, 1]] - p[tri[:, 0]], p[tri[:, 2]] - p[tri[:, 0]])
    assert np.all(nrm @ UP > 0)


def test_co_circular_tie_makes_a_hole_not_an_overlap():
    # Unit square: all four corners on one circle, both diagonals are ties.
    p = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    tri = triangulate(p, np.tile(UP, (4, 1)), spacing=1.0)
    assert len(tri) == 0


def test_fans_drop_large_rim_slivers():
    # Nearly collinear trio: huge circumcircle -> no triangle.
    p = np.array([[0.0, 0, 0], [1, 0, 0], [0.5, 0.01, 0]])
    assert len(triangulate(p, np.tile(UP, (3, 1)))) == 0


def test_empty_inputs():
    z = np.zeros((0, 3))
    assert agree_triangles(np.zeros((0, 3), np.int64)).shape == (0, 3)
    assert local_delaunay_fans(z, z, np.zeros((0, 2), np.int64), 1.0).shape == (0, 3)
    assert valence(np.zeros((0, 3), np.int64), 4).tolist() == [0, 0, 0, 0]
    assert open_edges(np.zeros((0, 3), np.int64)).shape == (0, 2)
    assert dual_edges(np.zeros((0, 3), np.int64)).shape == (0, 2)
    assert circumcenters(z, np.zeros((0, 3), np.int64)).shape == (0, 3)


def _square_grid(n=5):
    """n x n unit squares in z=0, each split along its (0,0)-(1,1) diagonal."""
    pts = np.array([(x, y, 0.0) for y in range(n + 1) for x in range(n + 1)], dtype=float)
    v = lambda x, y: y * (n + 1) + x  # noqa: E731
    tris = []
    for y in range(n):
        for x in range(n):
            tris += [[v(x, y), v(x + 1, y), v(x + 1, y + 1)], [v(x, y), v(x + 1, y + 1), v(x, y + 1)]]
    return pts, np.array(tris)


def test_fill_small_holes_closes_triangle_and_quad_holes():
    pts, full = _square_grid()
    rim = len(open_edges(full))
    middle = 2 * (2 * 5 + 2)                       # the square at (2, 2)
    one = np.delete(full, middle, axis=0)          # a 3-edge hole
    assert len(open_edges(one)) == rim + 3
    out = fill_small_holes(one, pts, 2.0)
    assert len(out) == len(full) and len(open_edges(out)) == rim
    two = np.delete(full, [middle, middle + 1], axis=0)   # a 4-edge hole
    out = fill_small_holes(two, pts, 2.0)
    assert len(out) == len(full) and len(open_edges(out)) == rim
    assert len(overused_edges(out)) == 0


def test_fill_small_holes_leaves_big_gaps_and_rims():
    pts = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
    tri = np.array([[0, 1, 2]])
    assert len(fill_small_holes(tri, pts, 2.0)) == 1   # the rim of a lone triangle is not a hole
    far = np.array([[0.0, 0, 0], [10, 0, 0], [0, 10, 0], [10, 10, 0]])
    two = np.array([[0, 1, 3]])
    assert len(fill_small_holes(two, far, 1.0)) == 1


def test_agree_triangles_needs_all_three_corners():
    proposals = np.array([[0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 3], [2, 3, 0]])
    assert agree_triangles(proposals).tolist() == [[0, 1, 2]]


# ---------- dual ----------

def test_circumcentre_of_right_triangle_is_hypotenuse_midpoint():
    p = np.array([[0.0, 0, 0], [2, 0, 0], [0, 2, 0]])
    tri = np.array([[0, 1, 2]])
    c = circumcenters(p, tri)
    assert np.allclose(c, [[1, 1, 0]])
    assert circumradii(p, tri, c)[0] == pytest.approx(math.sqrt(2))


def test_perfect_lattice_cells_are_regular_hexagons(hexgrid):
    p, _, tri = hexgrid
    c = circumcenters(p, tri)
    assert np.allclose(circumradii(p, tri, c), 1 / math.sqrt(3))
    e = dual_edges(tri)
    assert np.allclose(np.linalg.norm(c[e[:, 0]] - c[e[:, 1]], axis=1), 1 / math.sqrt(3))
    # One cell edge per interior triangulation edge.
    assert len(e) == int(np.sum(edge_use_counts(tri) == 2))
    assert np.all(e[:, 0] < e[:, 1])


# ---------- bristle sites ----------

def test_perfect_lattice_alternates_without_conflict(hexgrid):
    _, _, tri = hexgrid
    e = dual_edges(tri)
    parity = triangle_parity(len(tri), e)
    assert count_conflicts(parity, e) == 0
    mask = select_bristles(parity, "ALTERNATING", 0.0, np.random.default_rng(0))
    assert mask.sum() == len(tri) // 2
    # No two bristles on adjacent corners.
    assert not np.any(mask[e[:, 0]] & mask[e[:, 1]])


def test_odd_cycle_has_minimal_conflicts():
    # A 5-cycle of triangles (around a pentagonal cell) cannot be 2-coloured.
    e = np.array([[0, 1], [1, 2], [2, 3], [3, 4], [0, 4]])
    assert count_conflicts(triangle_parity(5, e), e) == 1


def test_seam_pairs_two_odd_cells_by_the_shortest_path():
    # Centres 0-1-2-3 in a chain plus a detour 0-4-5-3; cells 0 and 3 are odd.
    crossed = np.array([[0, 1], [1, 2], [2, 3], [0, 4], [4, 5], [5, 6], [6, 3]])
    odd = np.zeros(7, bool)
    odd[[0, 3]] = True
    cut = seam_cuts(crossed, odd, np.zeros(7, bool))
    assert cut.tolist() == [True, True, True, False, False, False, False]


def test_seam_runs_to_rim_when_alone():
    crossed = np.array([[0, 1], [1, 2]])
    odd = np.array([True, False, False])
    rim = np.array([False, False, True])
    assert seam_cuts(crossed, odd, rim).tolist() == [True, True]


def test_parity_respects_cuts_on_odd_cycle():
    e = np.array([[0, 1], [1, 2], [2, 3], [3, 4], [0, 4]])
    cuts = np.array([False, False, True, False, False])
    parity = triangle_parity(5, e, cuts)
    same = parity[e[:, 0]] == parity[e[:, 1]]
    assert same.tolist() == cuts.tolist()


def test_select_bristles_patterns():
    rng = np.random.default_rng(5)
    parity = np.array([0, 1] * 5000, dtype=np.int8)
    assert not select_bristles(parity, "NONE", 1.0, rng).any()
    assert select_bristles(parity, "RANDOM", 0.0, rng).sum() == 0
    assert select_bristles(parity, "RANDOM", 0.3, rng).mean() == pytest.approx(0.3, abs=0.02)
    with pytest.raises(ValueError):
        select_bristles(parity, "SOME", 0.3, rng)


# ---------- compound on the analytic surfaces ----------

SURFACES = {
    "sphere": SphereSurface(1.0),
    "ellipsoid": EllipsoidSurface(1.0, 0.7, 0.5),
    "cap": CapSurface(1.0, np.pi / 3),
}


@pytest.fixture(scope="module", params=list(SURFACES))
def built(request):
    surface = SURFACES[request.param]
    lattice, _ = build_lattice(surface, Config())
    topology, metrics = build_topology(lattice, surface, Config())
    return request.param, surface, lattice, topology, metrics


def test_hexagon_fraction(built):
    """Growing a perfect patch and pushing defects into the seams gives 0.95-0.96
    at 750 (relaxation alone reached 0.88-0.91; logbook 2026-11-03)."""
    _, _, _, _, m = built
    assert m["HEX_FRACTION"] >= 0.93
    assert m["N_OTHER_DEFECTS"] <= 0.01 * m["N_TOPO_INTERIOR"]  # 4- or 8-sided cells are rare


def test_closed_surfaces_are_watertight_with_euler_two(built):
    name, _, lattice, topo, m = built
    if name == "cap":
        assert m["N_OPEN_EDGES"] > 0
        return
    assert m["HOLES_ON_CLOSED_SURFACE"] == 0 and m["N_OVERLAP_EDGES"] == 0
    assert edge_use_counts(topo.triangles).max() == 2
    n_edges = len(edge_use_counts(topo.triangles))
    assert len(lattice) - n_edges + len(topo.triangles) == 2
    # Euler on a sphere: pentagons − heptagons = 12 (every centre is interior)
    assert m["N_PENTAGONS"] - m["N_HEPTAGONS"] == 12


def test_no_hole_larger_than_one_cell(built):
    _, _, _, _, m = built
    assert m["MAX_HOLE"] < 1.0


def test_cell_points_on_surface(built):
    _, surface, _, topo, _ = built
    x, _ = surface.closest_point(topo.cell_points)
    assert np.allclose(x, topo.cell_points, atol=1e-9)
    assert topo.cell_edges.max() < len(topo.triangles)


def test_about_one_bristle_per_ommatidium(built):
    """Drosophila: each ommatidium is accompanied by a single bristle (atlas)."""
    name, _, lattice, topo, m = built
    ratio = m["N_BRISTLES"] / len(lattice)
    assert 0.85 <= ratio <= 1.1 if name != "cap" else 0.75 <= ratio <= 1.1
    # Seams pair the 5- and 7-sided cells, so conflicts stay near one per defect.
    assert m["PARITY_CONFLICTS"] < 0.05 * len(topo.cell_edges)
    assert m["PARITY_CONFLICTS"] <= 1.5 * (m["N_PENTAGONS"] + m["N_HEPTAGONS"]) + 20


def test_topology_of_empty_lattice():
    empty = Lattice(np.zeros((0, 3)), np.zeros((0, 3)), 1.0)
    topo, m = build_topology(empty, SphereSurface(1.0), Config())
    assert len(topo.triangles) == 0 and m["N_BRISTLES"] == 0


def test_two_skins_are_never_joined():
    """Two parallel sheets 0.01 apart facing opposite ways (a thin solid): no
    triangle or wall mixes them."""
    from ommatidia.lattice.lattice import Lattice
    from fixtures.analytic_surfaces import CapSurface

    g = np.array([(x + 0.5 * (y % 2), y * np.sqrt(3) / 2) for y in range(8) for x in range(8)])
    top = np.column_stack([g, np.full(len(g), 0.005)])
    bottom = np.column_stack([g + 0.3, np.full(len(g), -0.005)])
    centers = np.concatenate([top, bottom])
    normals = np.concatenate([np.tile([0.0, 0, 1], (len(g), 1)), np.tile([0.0, 0, -1], (len(g), 1))])

    class Flat:
        def closest_point(self, p):
            p = np.asarray(p, float).reshape(-1, 3)
            return p, np.where(p[:, 2:] >= 0, [[0.0, 0, 1]], [[0.0, 0, -1]])

        def boundary_distance(self, p):
            return np.full(len(p), np.inf)

    topo, _ = build_topology(Lattice(centers, normals, 1.0), Flat(), Config())
    layer = (np.arange(len(centers)) >= len(g)).astype(int)
    assert len(topo.triangles) > 0
    assert np.all(layer[topo.triangles].min(axis=1) == layer[topo.triangles].max(axis=1))
