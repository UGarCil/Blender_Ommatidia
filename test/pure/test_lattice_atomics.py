"""Atomic algorithms of lattice/: density, neighbours, sampling, growth, relaxation, metrics."""

import math

import numpy as np
import pytest

from fixtures.analytic_surfaces import SphereSurface
from ommatidia.lattice.atomic_algorithm_scripts.density import HEX_CELL_FACTOR, resolve_density
from ommatidia.lattice.atomic_algorithm_scripts.metrics import nearest_neighbor_distances, spacing_cv
from ommatidia.lattice.atomic_algorithm_scripts.neighbors import assign_nearest, neighbor_pairs
from ommatidia.lattice.atomic_algorithm_scripts.relax import (
    jitter_tangent,
    repulsion_step,
    trim_boundary,
)
from ommatidia.lattice.atomic_algorithm_scripts.growth import (
    DEAD,
    HEX_DIRECTIONS,
    apex_points,
    apex_triples,
    auto_origin,
    first_ring,
    hex_ring,
    mean_by_owner,
    node_links,
    patch_edges,
    placement_mask,
    ring_candidates,
    spring_step,
    tangent_frame,
)
from ommatidia.lattice.atomic_algorithm_scripts.neighbors import nearest_distances
from ommatidia.lattice.atomic_algorithm_scripts.sampling import fill_mask, poisson_filter


def planar_hex_grid(rows=10, cols=10, a=1.0):
    """A perfect triangular lattice of centres in the z = 0 plane."""
    pts = [(c * a + (r % 2) * a / 2, r * a * math.sqrt(3) / 2, 0.0)
           for r in range(rows) for c in range(cols)]
    return np.array(pts)


# ---------- density ----------

def test_density_spacing_mode_counts_hex_cells():
    spacing, n = resolve_density(HEX_CELL_FACTOR * 100.0, "SPACING", 1.0, 1, 10**6)
    assert (spacing, n) == (1.0, 100)


def test_density_count_mode_inverts_spacing_mode():
    spacing, n = resolve_density(4 * math.pi, "COUNT", 99.0, 750, 10**6)
    assert n == 750
    assert resolve_density(4 * math.pi, "SPACING", spacing, 1, 10**6)[1] == 750


def test_density_never_below_one():
    assert resolve_density(1e-12, "SPACING", 1.0, 1, 10)[1] == 1


def test_density_max_count_guard():
    with pytest.raises(ValueError, match="model_scale"):
        resolve_density(1.0, "SPACING", 16e-6, 1, 50_000)


@pytest.mark.parametrize("area", [0.0, -1.0, float("nan")])
def test_density_rejects_bad_area(area):
    with pytest.raises(ValueError):
        resolve_density(area, "COUNT", 1.0, 10, 100)


# ---------- neighbours ----------

def brute_pairs(p, r):
    d = np.linalg.norm(p[:, None] - p[None], axis=2)
    i, j = np.nonzero(np.triu(d < r, 1))
    return np.column_stack([i, j])


@pytest.mark.parametrize("n", [0, 1, 2, 300])
def test_neighbor_pairs_match_brute_force(n):
    p = np.random.default_rng(n).uniform(-1, 1, size=(n, 3))
    got = neighbor_pairs(p, 0.3)
    assert got.shape[1] == 2 and got.dtype == np.int64
    assert np.array_equal(got, brute_pairs(p, 0.3).reshape(-1, 2))


def test_neighbor_pairs_invariants():
    p = np.random.default_rng(1).normal(size=(500, 3))
    pr = neighbor_pairs(p, 0.5)
    assert np.all(pr[:, 0] < pr[:, 1])
    assert len(np.unique(pr, axis=0)) == len(pr)
    assert np.array_equal(pr[np.lexsort((pr[:, 1], pr[:, 0]))], pr)


def test_neighbor_pairs_far_apart_points_and_negative_coordinates():
    p = np.array([[-1e3, 0, 0], [-1e3 + 0.1, 0, 0], [1e3, 0, 0]])
    assert neighbor_pairs(p, 0.2).tolist() == [[0, 1]]


def test_assign_nearest_matches_brute_force_including_far_samples():
    rng = np.random.default_rng(2)
    centers = rng.uniform(0, 1, size=(50, 3))
    samples = np.concatenate([rng.uniform(0, 1, size=(400, 3)), [[10.0, 10.0, 10.0]]])
    got = assign_nearest(samples, centers, 0.2)
    want = np.linalg.norm(samples[:, None] - centers[None], axis=2).argmin(axis=1)
    assert np.array_equal(got, want)


def test_assign_nearest_requires_centres():
    with pytest.raises(ValueError):
        assign_nearest(np.zeros((3, 3)), np.zeros((0, 3)), 1.0)


# ---------- sampling ----------

def test_poisson_filter_respects_radius_and_is_maximal():
    c = np.random.default_rng(3).uniform(0, 1, size=(2000, 3))
    kept = poisson_filter(c, 0.1)
    k = c[kept]
    assert len(neighbor_pairs(k, 0.1)) == 0
    # every rejected candidate is within the radius of a kept one
    d = np.linalg.norm(c[~kept][:, None] - k[None], axis=2).min(axis=1)
    assert np.all(d < 0.1)


def test_poisson_filter_empty():
    assert poisson_filter(np.zeros((0, 3)), 1.0).shape == (0,)


def test_fill_mask_keeps_existing_and_respects_radius():
    existing = np.array([[0.0, 0, 0], [1.0, 0, 0]])
    up = np.tile([0.0, 0, 1], (2, 1))
    cand = np.array([[0.2, 0, 0], [3.0, 0, 0], [3.1, 0, 0], [5.0, 0, 0]])
    keep = fill_mask(cand, existing, np.tile([0.0, 0, 1], (4, 1)), up, 0.5)
    assert keep.tolist() == [False, True, False, True]


def test_fill_mask_ignores_points_facing_the_other_way():
    """A thin sheet: a point on the back face does not block a front-face candidate."""
    keep = fill_mask(np.array([[0.0, 0, 0.01]]), np.array([[0.0, 0, 0]]),
                     np.array([[0.0, 0, 1]]), np.array([[0.0, 0, -1]]), 0.5)
    assert keep.tolist() == [True]


def test_nearest_distances():
    d = nearest_distances(np.array([[0.0, 0, 0], [5, 0, 0]]), np.array([[1.0, 0, 0]]), 2.0)
    assert np.allclose(d, [1, 4])
    assert np.all(np.isinf(nearest_distances(np.zeros((2, 3)), np.zeros((0, 3)), 1.0)))


# ---------- growth ----------

def test_hex_directions_are_ccw_and_consecutive_ones_are_neighbours():
    for t_ in range(6):
        a, b = HEX_DIRECTIONS[t_], HEX_DIRECTIONS[(t_ + 1) % 6]
        assert hex_ring(tuple(a)) == 1 and hex_ring(tuple(b - a)) == 1


def test_hex_ring():
    assert hex_ring((0, 0)) == 0 and hex_ring((2, -1)) == 2 and hex_ring((-3, 3)) == 3


@pytest.mark.parametrize("rotation", [0.0, 37.0])
def test_tangent_frame_is_orthonormal_and_tangent(rotation):
    n = np.array([0.3, -0.5, 0.8]) / np.linalg.norm([0.3, -0.5, 0.8])
    e1, e2 = tangent_frame(n, rotation)
    assert np.allclose([e1 @ e1, e2 @ e2, e1 @ e2, e1 @ n, e2 @ n], [1, 1, 0, 0, 0])
    assert np.allclose(np.cross(e1, e2), n)


def test_tangent_frame_rotation_turns_rows():
    n = np.array([0.0, -1.0, 0.0])
    e1, _ = tangent_frame(n, 0.0)
    r1, _ = tangent_frame(n, 90.0)
    assert np.allclose(e1, [0, 0, 1]) and abs(e1 @ r1) < 1e-12


def test_first_ring_is_a_regular_hexagon():
    ring = first_ring(np.zeros(3), np.array([1.0, 0, 0]), np.array([0.0, 1, 0]), 2.0)
    assert np.allclose(np.linalg.norm(ring, axis=1), 2.0)
    assert np.allclose(np.linalg.norm(ring - np.roll(ring, 1, axis=0), axis=1), 2.0)


def test_ring_candidates_and_apex_triples_on_a_flat_patch():
    index = {(0, 0): 0}
    index.update({tuple(int(x) for x in d): k + 1 for k, d in enumerate(HEX_DIRECTIONS)})
    ring2 = ring_candidates(list(index)[1:], index, 2)
    assert len(ring2) == 12 and all(hex_ring(n) == 2 for n in ring2)
    owner, a, b, c = apex_triples(ring2, index)
    # the 6 edge nodes of ring 2 have one placed pair each; the 6 corners none yet
    assert len(np.unique(owner)) == 6 and len(owner) == 6
    assert np.all(c == 0)


def test_ring_candidates_skip_dead_front_nodes():
    index = {(0, 0): 0, (1, 0): DEAD}
    assert ring_candidates([(1, 0)], index, 2) == []


def test_apex_points_build_equilateral_triangles_away_from_c():
    pts = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.5, -np.sqrt(3) / 2, 0]])
    nrm = np.tile([0.0, 0, 1], (3, 1))
    apex = apex_points(pts, nrm, np.array([0]), np.array([1]), np.array([2]), 1.0)
    assert np.allclose(apex, [[0.5, np.sqrt(3) / 2, 0]])


def test_mean_by_owner():
    m, has = mean_by_owner(np.array([[0.0, 0, 0], [2, 0, 0], [5, 5, 5]]), np.array([0, 0, 2]), 3)
    assert has.tolist() == [0, 2] and np.allclose(m, [[1, 0, 0], [5, 5, 5]])


def _placement(new, existing, links=None, batch=None, bd=None):
    new, existing = np.asarray(new, float), np.asarray(existing, float)
    up = lambda x: np.tile([0.0, 0, 1], (len(x), 1))  # noqa: E731
    bd = np.full(len(new), np.inf) if bd is None else np.asarray(bd)
    kw = {}
    if links is not None:
        kw["links"] = np.asarray(links).reshape(-1, 2)
    if batch is not None:
        kw["batch_links"] = np.asarray(batch).reshape(-1, 2)
    return placement_mask(new, new, up(new), bd, existing, up(existing), 1.0, **kw)


def test_placement_kills_overlap_rim_and_batch_duplicates():
    alive = _placement([[0.65, 0, 0], [5, 0, 0], [5.3, 0, 0], [9, 0, 0]], [[0.0, 0, 0]],
                       bd=[np.inf, np.inf, np.inf, 0.2])
    assert alive.tolist() == [False, True, False, False]


def test_placement_lets_compressed_lattice_neighbours_live():
    assert _placement([[0.6, 0, 0]], [[0.0, 0, 0]], links=[[0, 0]]).tolist() == [True]
    assert _placement([[0.5, 0, 0]], [[0.0, 0, 0]], links=[[0, 0]]).tolist() == [False]
    assert _placement([[5, 0, 0], [5.6, 0, 0]], [[0.0, 0, 0]],
                      batch=[[0, 1], [1, 0]]).tolist() == [True, True]


def test_placement_kills_nodes_that_jump_to_another_skin():
    """A proposal that lands on a surface facing away from its parents (the back of a
    thin solid) dies; a gentle turn is fine."""
    new = np.array([[0.0, 0, 0], [5.0, 0, 0]])
    new_normals = np.array([[0.0, 0, -1], [0.0, 0.5, 0.866]])
    parents = np.tile([0.0, 0, 1], (2, 1))
    alive = placement_mask(new, new, new_normals, np.full(2, np.inf), np.zeros((0, 3)),
                           np.zeros((0, 3)), 1.0, parent_normals=parents)
    assert alive.tolist() == [False, True]


def test_node_links():
    index = {(0, 0): 0, (1, 0): 1, (0, 1): DEAD}
    placed, among = node_links([(1, -1), (2, -1)], index, 10)
    assert {tuple(x) for x in placed} == {(0, 10), (0, 11), (1, 11)}
    assert {tuple(x) for x in among} == {(0, 1), (1, 0)}


def test_patch_edges_lists_each_neighbour_pair_once():
    index = {(0, 0): 0}
    index.update({tuple(int(x) for x in d): k + 1 for k, d in enumerate(HEX_DIRECTIONS)})
    index[(1, 0)] = DEAD
    edges = patch_edges(index)
    assert len(edges) == 12 - 3  # full first ring has 12 edges; the dead node had 3
    assert len({tuple(sorted(e)) for e in edges}) == len(edges)


def test_spring_step_restores_rest_length_and_keeps_fixed_points():
    p = np.array([[0.0, 0, 0], [1.5, 0, 0]])
    nrm = np.tile([0.0, 0, 1], (2, 1))
    v = np.zeros_like(p)
    movable = np.array([False, True])
    for _ in range(100):
        p, v = spring_step(p, v, nrm, np.array([[0, 1]]), movable, 1.0)
    assert np.allclose(p[0], 0) and np.isclose(p[1, 0], 1.0, atol=1e-3) and p[1, 2] == 0


def test_auto_origin_open_and_closed():
    samples = np.array([[0.0, 0, 1], [0, -1, 0], [0, 1, 0]])
    assert np.allclose(auto_origin(samples, np.array([0.1, 0.5, 0.2])), [0, -1, 0])
    assert np.allclose(auto_origin(samples, np.full(3, np.inf)), [0, -1, 0])


# ---------- relaxation ----------

def test_repulsion_pushes_close_pair_apart_tangentially():
    p = np.array([[0.0, 0, 0], [0.5, 0, 0]])
    nrm = np.array([[0.0, 0, 1], [0.0, 0, 1]])
    out, v = repulsion_step(p, np.zeros_like(p), nrm, np.array([[0, 1]]), 1.0)
    assert out[1, 0] - out[0, 0] > 0.5
    assert np.allclose(out[:, 2], 0.0)          # no normal motion
    assert np.allclose(out[0] + out[1], p[0] + p[1])  # forces cancel
    assert np.allclose(v, out - p)


def test_repulsion_leaves_perfect_lattice_interior_still():
    p = planar_hex_grid(12, 12)
    nrm = np.tile([0.0, 0, 1], (len(p), 1))
    out, _ = repulsion_step(p, np.zeros_like(p), nrm, neighbor_pairs(p, 2.0), 1.0)
    centre = np.linalg.norm(p - p.mean(axis=0), axis=1) < 2.0
    assert np.allclose(out[centre], p[centre], atol=1e-9)


def test_repulsion_move_is_clamped():
    p = np.array([[0.0, 0, 0], [1e-9, 0, 0]])
    nrm = np.array([[0.0, 0, 1], [0.0, 0, 1]])
    out, _ = repulsion_step(p, np.zeros_like(p), nrm, np.array([[0, 1]]), 1.0)
    assert np.all(np.linalg.norm(out - p, axis=1) <= 0.15 + 1e-12)


def test_repulsion_without_pairs_keeps_momentum_only():
    p = np.zeros((1, 3))
    nrm = np.array([[0.0, 0, 1]])
    out, v = repulsion_step(p, np.array([[0.01, 0, 0.5]]), nrm, np.zeros((0, 2), np.int64), 1.0)
    assert np.allclose(v, [[0.009, 0, 0]])


def test_trim_boundary():
    mask = trim_boundary(np.array([0.0, 0.49, 0.5, np.inf]), 1.0)
    assert mask.tolist() == [False, False, True, True]


def test_jitter_is_tangential_and_scaled():
    rng = np.random.default_rng(4)
    p = np.zeros((20000, 3))
    nrm = np.tile([0.0, 0, 1], (len(p), 1))
    out = jitter_tangent(p, nrm, 0.1, 2.0, rng)
    assert np.allclose(out[:, 2], 0)
    assert np.sqrt(np.mean(np.sum(out**2, axis=1))) == pytest.approx(0.2, rel=0.03)
    assert np.array_equal(jitter_tangent(p, nrm, 0.0, 2.0, rng), p)


# ---------- metrics ----------

def test_metrics_on_perfect_lattice():
    p = planar_hex_grid()
    nn = nearest_neighbor_distances(p, neighbor_pairs(p, 1.5))
    assert np.allclose(nn, 1.0)
    assert spacing_cv(nn, np.ones(len(p), bool)) == pytest.approx(0.0, abs=1e-12)


def test_metrics_isolated_point_is_inf_and_ignored():
    p = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0], [50.0, 0, 0]])
    nn = nearest_neighbor_distances(p, neighbor_pairs(p, 1.5))
    assert nn.tolist() == [1.0, 1.0, 1.0, np.inf]
    assert spacing_cv(nn, np.ones(4, bool)) == 0.0
    assert math.isnan(spacing_cv(nn, np.array([True, False, False, True])))


def test_sphere_fixture_used_by_atomics_is_consistent():
    # Guard: the sphere stub of the DDs is what the compound tests rely on.
    assert SphereSurface(1.0).area() == pytest.approx(4 * math.pi)
