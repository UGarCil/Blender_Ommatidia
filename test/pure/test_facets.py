"""facets/: corner loops, facet layout, dome profiles, bristles, and the composite."""

import numpy as np
import pytest

from fixtures.analytic_surfaces import CapSurface, SphereSurface
from ommatidia.config import Config
from ommatidia.facets.atomic_algorithm_scripts.bristles import (
    SINK,
    bristle_positions,
    bristle_template,
    spindle_radius,
)
from ommatidia.facets.atomic_algorithm_scripts.cell_loops import cell_loops
from ommatidia.facets.atomic_algorithm_scripts.dome_mesh import facet_layout, facet_resolution
from ommatidia.facets.atomic_algorithm_scripts.dome_profile import crease_width, dome_heights, rim_fade
from ommatidia.facets.facets import (
    KIND_BRISTLE,
    KIND_FACET,
    LIFT,
    build_facet_mesh,
    real_length_scale,
    shape_positions,
)
from ommatidia.lattice.lattice import build_lattice
from ommatidia.topology.topology import build_topology


def _faces(mesh_or_layout, kind=None):
    starts = np.append(mesh_or_layout.loop_start, len(mesh_or_layout.loop_vertices))
    rows = range(len(mesh_or_layout.loop_start))
    if kind is not None:
        rows = np.flatnonzero(mesh_or_layout.face_kind == kind)
    return [mesh_or_layout.loop_vertices[starts[p]:starts[p + 1]] for p in rows]


def _edge_counts(faces):
    edges = np.array([(a, b) for f in faces for a, b in zip(f, np.roll(f, -1))])
    undirected, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    _, directed_counts = np.unique(edges, axis=0, return_counts=True)
    return counts, directed_counts


@pytest.fixture(scope="module", params=["sphere", "cap"])
def built(request):
    surface = SphereSurface(1.0) if request.param == "sphere" else CapSurface(1.0, np.pi / 3)
    config = Config(target_count=300)
    lattice, _ = build_lattice(surface, config)
    topology, _ = build_topology(lattice, surface, config)
    return request.param, surface, lattice, topology, build_facet_mesh(lattice, topology, surface, config)


# ---------- cell loops ----------

def test_cell_loops_close_around_interior_centres(built):
    name, _, lattice, topology, _ = built
    loops = cell_loops(topology.triangles, len(lattice))
    has = np.array([len(x) > 0 for x in loops])
    assert np.all([len(x) == v for x, v, h in zip(loops, topology.valence, has) if h])
    if name == "sphere":
        assert has.all()
    else:
        assert 0.7 < has.mean() < 1.0      # the rim row has no facet
    # consecutive corners share a facet wall
    walls = {tuple(e) for e in topology.cell_edges.tolist()}
    for loop in loops[:50]:
        for a, b in zip(loop, np.roll(loop, -1)):
            assert (min(a, b), max(a, b)) in walls


def test_cell_loops_counter_clockwise(built):
    _, _, lattice, topology, _ = built
    loops = cell_loops(topology.triangles, len(lattice))
    for i, loop in enumerate(loops[:100]):
        if len(loop) == 0:
            continue
        c = topology.cell_points[loop] - lattice.centers[i]
        turn = np.cross(c, np.roll(c, -1, axis=0)) @ lattice.normals[i]
        assert np.all(turn > 0)


def test_cell_loops_empty():
    loops = cell_loops(np.zeros((0, 3), dtype=np.int64), 2)
    assert len(loops) == 2 and all(len(x) == 0 for x in loops)


# ---------- layout ----------

@pytest.mark.parametrize("level, expected", [(-1, (3, 6)), (0, (1, 2)), (3, (4, 8))])
def test_facet_resolution(level, expected):
    assert facet_resolution(level) == expected


def test_layout_rho_ranges_and_radius(built):
    _, _, lattice, topology, _ = built
    loops = cell_loops(topology.triangles, len(lattice))
    layout, radius = facet_layout(lattice.centers, topology.cell_points, loops, topology.cell_edges, 3, 6)
    inside = layout.owner >= 0
    assert np.all(layout.rho_hex[~inside] == 1) and np.all(layout.rho_hex[inside] < 1)
    # corners of a regular hexagon are at 2/sqrt(3) = 1.15 inner radii; strained ones a bit more
    assert np.all(layout.rho_round[inside] >= 0) and np.percentile(layout.rho_round[inside], 99.9) < 1.3
    built_cells = np.array([len(x) > 0 for x in loops])
    # a regular hexagon with centre spacing s has inner radius s/2
    assert np.allclose(np.median(radius[built_cells]) / lattice.spacing, 0.5, atol=0.05)
    assert np.all(radius[~built_cells] == 0)


def test_layout_vertex_count(built):
    _, _, lattice, topology, _ = built
    loops = cell_loops(topology.triangles, len(lattice))
    sides, rings = 2, 3
    layout, _ = facet_layout(lattice.centers, topology.cell_points, loops, topology.cell_edges, sides, rings)
    per_cell = sum(1 + (rings - 1) * len(x) * sides for x in loops if len(x))
    assert np.sum(layout.owner >= 0) == per_cell


def test_layout_skips_a_facet_whose_wall_is_missing(built):
    """Overlapping triangles at a hard crease can leave a corner loop without one
    of its walls: that facet is skipped instead of crashing."""
    _, _, lattice, topology, _ = built
    loops = cell_loops(topology.triangles, len(lattice))
    k = next(i for i, x in enumerate(loops) if len(x))
    a, b = sorted(loops[k][:2])
    edges = topology.cell_edges[~((topology.cell_edges[:, 0] == a) & (topology.cell_edges[:, 1] == b))]
    layout, radius = facet_layout(lattice.centers, topology.cell_points, loops, edges, 2, 2)
    assert radius[k] == 0 and np.all(layout.owner != k)


def test_layout_without_facets():
    layout, radius = facet_layout(np.zeros((2, 3)), np.zeros((0, 3)), [np.zeros(0, int)] * 2,
                                  np.zeros((0, 2), int), 3, 6)
    assert len(layout.points) == 0 and np.all(radius == 0)


# ---------- dome profile ----------

R = np.full(5, 1.0)
U = np.array([0.0, 0.3, 0.6, 0.9, 1.0])


def test_dome_zero_at_the_wall_and_peak_at_the_centre():
    h = dome_heights(U, U, R, inflation=0.4, pokeyness=0.0, facet_shape=0.0, crease_sharpness=1.0)
    assert np.isclose(h[0], 0.4) and h[-1] == 0 and np.all(np.diff(h) <= 0)


def test_inflation_one_is_a_semi_dome():
    """Sharp crease, round dome, inflation 1: a hemisphere of the facet radius."""
    u = np.linspace(0, 0.9, 10)
    h = dome_heights(u, u, np.ones(10), 1.0, 0.0, 0.0, 1.0)
    inner = u < 1 - 2 * crease_width(1.0)
    assert np.allclose(h[inner], np.sqrt(1 - u[inner] ** 2))


def test_inflation_zero_is_flat():
    assert np.all(dome_heights(U, U, R, 0.0, 0.5, 0.5, 0.5) == 0)


def test_pokeyness_one_is_a_cone():
    h = dome_heights(U, U, R, 0.5, 1.0, 0.0, 1.0)
    inner = U < 1 - 2 * crease_width(1.0)
    assert np.allclose(h[inner], 0.5 * (1 - U[inner]))


def test_flat_facet_is_a_plateau_with_hexagonal_walls():
    """facet_shape 1 follows rho_hex: corners (rho_round > 1) still reach full height
    inside the hexagon."""
    rho_round = np.array([0.0, 0.5, 1.1])
    rho_hex = np.array([0.0, 0.4, 0.6])
    h = dome_heights(rho_round, rho_hex, np.ones(3), 0.3, 0.0, 1.0, 1.0)
    assert np.allclose(h, 0.3)


def test_softer_crease_lowers_the_dome_near_the_wall():
    sharp = dome_heights(U, U, R, 0.4, 0.0, 0.0, 1.0)
    soft = dome_heights(U, U, R, 0.4, 0.0, 0.0, 0.0)
    assert soft[3] < sharp[3] and np.isclose(soft[0], sharp[0])


def test_rim_fade():
    d = np.array([0.5, 1.0, 2.0, 1e6])
    assert np.allclose(rim_fade(d, 1.5)[[0, 2, 3]], [0, 1, 1]) and 0 < rim_fade(d, 1.5)[1] < 1
    assert np.all(rim_fade(d, 0.0) == 1)


# ---------- bristles ----------

def test_bristle_template():
    axial, angle, faces = bristle_template(6, 4)
    assert len(axial) == 6 * 4 + 1 and axial[-1] == 1 and axial[0] == 0
    assert faces.shape == (6 * 3 + 6, 4) and faces.max() == 24
    assert np.sum(faces[:, 3] == -1) == 6      # tip triangles


def test_spindle_is_widest_near_a_quarter_and_pointed():
    a = np.linspace(0, 1, 101)
    r = spindle_radius(a)
    assert r[-1] == 0 and 0.15 < a[np.argmax(r)] < 0.35 and np.isclose(r.max(), 1.0, atol=0.02)


def test_bristle_axis_tilts_toward_the_comb():
    site, n = np.zeros((2, 3)), np.tile([0.0, 0, 1], (2, 1))
    axial, angle = np.array([0.0, 1.0]), np.zeros(2)
    up = bristle_positions(site, n, axial, angle, 10.0, 1.0, 0.0, np.array([1.0, 0, 0]))
    assert np.allclose(up[1], [0, 0, 10 - SINK])
    tilted = bristle_positions(site, n, axial, angle, 10.0, 1.0, 30.0, np.array([1.0, 0, -1]))
    tip = tilted[1] + SINK * n[1]            # tip relative to the site
    assert np.isclose(np.degrees(np.arctan2(tip[0], tip[2])), 30.0, atol=1e-6) and abs(tip[1]) < 1e-9


def test_bristle_comb_parallel_to_normal_still_works():
    p = bristle_positions(np.zeros((1, 3)), np.array([[0.0, 0, 1]]), np.array([1.0]), np.zeros(1),
                          1.0, 0.1, 45.0, np.array([0.0, 0, -1]))
    assert np.all(np.isfinite(p))


# ---------- composite ----------

def test_composite_is_watertight_and_outward(built):
    name, _, _, _, mesh = built
    faces = _faces(mesh, KIND_FACET)
    counts, directed = _edge_counts(faces)
    assert counts.max() == 2 and directed.max() == 1        # manifold, consistently oriented
    if name == "sphere":
        assert np.all(counts == 2)                           # closed
    positions = shape_positions(mesh, Config())
    for f in faces[::50]:
        p = positions[f]
        assert np.cross(p[1] - p[0], p[2] - p[0]) @ mesh.normal[f[0]] > 0


def test_one_bristle_per_site(built):
    _, _, _, topology, mesh = built
    k = len(bristle_template()[0])
    assert np.sum(mesh.kind == KIND_BRISTLE) == k * topology.bristle_mask.sum()


def test_shape_moves_facets_only_along_normals(built):
    _, _, _, _, mesh = built
    positions = shape_positions(mesh, Config(inflation=0.7))
    f = mesh.kind == KIND_FACET
    offset = positions[f] - mesh.base[f]
    along = np.sum(offset * mesh.normal[f], axis=1)
    assert np.allclose(offset, along[:, None] * mesh.normal[f], atol=1e-12)
    assert along.min() >= LIFT * mesh.spacing * (1 - 1e-9)
    flat = shape_positions(mesh, Config(inflation=0.0))
    assert np.allclose(np.sum((flat[f] - mesh.base[f]) * mesh.normal[f], axis=1), LIFT * mesh.spacing)


def test_open_sheet_fades_at_the_rim(built):
    name, _, _, _, mesh = built
    if name != "cap":
        return
    positions = shape_positions(mesh, Config(rim_width=3.0))
    f = (mesh.kind == KIND_FACET) & (mesh.rho_round < 0.05)
    h = np.sum((positions[f] - mesh.base[f]) * mesh.normal[f], axis=1)
    near, far = mesh.rim[f] < 2.0, mesh.rim[f] > 4.0      # the rim row itself has no facet
    assert near.any() and far.any()
    assert h[near].max() < h[far].min()


def test_real_lengths_follow_the_facets_in_count_mode(built):
    """COUNT mode on a unit sphere: one spacing counts as 16 µm, so a 10 µm bristle
    is ~0.6 spacings long, however big the mesh is."""
    _, _, _, _, mesh = built
    config = Config(target_count=300)
    positions = shape_positions(mesh, config)
    b = mesh.kind == KIND_BRISTLE
    tips = positions[b][mesh.axial[b] == 1]
    bases = mesh.base[b][mesh.axial[b] == 1]
    length = np.linalg.norm(tips - bases, axis=1).mean() / mesh.spacing
    assert np.isclose(length, config.bristle_length / config.spacing, rtol=0.1)
    assert real_length_scale(Config(density_mode="SPACING"), mesh.spacing) == 1.0
