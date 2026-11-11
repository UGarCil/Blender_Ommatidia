"""Compound lattice/: regularity, coverage, invariants and determinism (AGENTS §8)."""

import numpy as np
import pytest

from fixtures.analytic_surfaces import CapSurface, EllipsoidSurface, SphereSurface
from ommatidia.config import Config
from ommatidia.lattice.atomic_algorithm_scripts.neighbors import assign_nearest, neighbor_pairs
from ommatidia.lattice.lattice import build_lattice, interior_mask

SURFACES = {
    "sphere": SphereSurface(1.0),
    "ellipsoid": EllipsoidSurface(1.0, 0.7, 0.5),   # the "deformed" fixture
    "cap": CapSurface(1.0, np.pi / 3),               # open sheet
}
# Spacing CV of interior centres. Dossier §5 asked 0.05 on a sphere; growing a
# perfect patch trades a little size variation (cells shrink by up to ~10%
# toward the seams) for far fewer defects (logbook 2026-11-03), so 0.07.
MAX_CV = {"sphere": 0.07, "ellipsoid": 0.08, "cap": 0.08}


@pytest.fixture(scope="module", params=list(SURFACES))
def built(request):
    surface = SURFACES[request.param]
    lattice, metrics = build_lattice(surface, Config())
    return request.param, surface, lattice, metrics


def test_regularity(built):
    name, _, _, metrics = built
    assert metrics["SPACING_CV"] <= MAX_CV[name]
    assert 0.9 <= metrics["MEAN_NN"] / metrics["SPACING"] <= 1.05


def test_centres_on_surface_with_unit_normals(built):
    _, surface, lattice, _ = built
    x, n = surface.closest_point(lattice.centers)
    assert np.allclose(x, lattice.centers, atol=1e-9 * lattice.spacing)
    assert np.allclose(lattice.normals, n)


def test_no_two_centres_closer_than_half_spacing(built):
    _, _, lattice, metrics = built
    assert len(neighbor_pairs(lattice.centers, 0.5 * lattice.spacing)) == 0
    assert metrics["MIN_NN"] >= 0.5 * lattice.spacing


def test_open_boundary_trimmed(built):
    _, surface, lattice, _ = built
    assert np.all(surface.boundary_distance(lattice.centers) >= 0.5 * lattice.spacing)


def test_coverage_no_hole_larger_than_one_cell(built):
    """Every interior surface point is within one SPACING of a centre (a perfect
    lattice gives SPACING/sqrt(3))."""
    _, surface, lattice, _ = built
    samples = surface.sample_area_weighted(20_000, np.random.default_rng(99))
    samples = samples[surface.boundary_distance(samples) > lattice.spacing]
    nearest = assign_nearest(samples, lattice.centers, lattice.spacing)
    gap = np.linalg.norm(samples - lattice.centers[nearest], axis=1)
    assert gap.max() < lattice.spacing


def test_count_close_to_target(built):
    """Curvature compresses the patch slightly, so a closed surface gets a few % more."""
    name, _, lattice, metrics = built
    assert metrics["N_TARGET"] == 750
    assert 0.93 <= len(lattice) / 750 <= 1.07


def test_one_patch_covers_a_single_surface(built):
    _, _, lattice, metrics = built
    assert metrics["N_PATCHES"] == 1
    assert np.mean(lattice.patch >= 0) > 0.95


def test_interior_mask_closed_vs_open(built):
    name, surface, lattice, _ = built
    mask = interior_mask(lattice, surface)
    assert mask.all() if name != "cap" else (0 < mask.sum() < len(mask))


def test_determinism_same_seed_identical_other_seed_differs():
    s = SphereSurface(1.0)
    cfg = Config(target_count=200, relax_iterations=20, irregularity=0.05)
    a, _ = build_lattice(s, cfg)
    b, _ = build_lattice(s, cfg)
    c, _ = build_lattice(s, Config(target_count=200, relax_iterations=20, irregularity=0.05, seed=1))
    assert np.array_equal(a.centers, b.centers)
    assert not np.array_equal(a.centers, c.centers)


def test_spacing_mode_scales_with_surface():
    """Real Drosophila units: 16 µm facets on a sphere of radius 170 µm."""
    s = SphereSurface(170e-6)
    lattice, metrics = build_lattice(s, Config(density_mode="SPACING", spacing=16e-6,
                                               relax_iterations=50))
    expected = round(4 * np.pi * 170e-6**2 / (np.sqrt(3) / 2 * 16e-6**2))  # = 1638
    assert metrics["SPACING"] == 16e-6
    assert metrics["N_TARGET"] == expected
    assert 0.95 <= len(lattice) / expected <= 1.07
    assert metrics["SPACING_CV"] < 0.07


def test_irregularity_adds_disorder_but_keeps_centres_on_surface():
    s = SphereSurface(1.0)
    _, m0 = build_lattice(s, Config(target_count=300, relax_iterations=100))
    lat, m1 = build_lattice(s, Config(target_count=300, relax_iterations=100, irregularity=0.1))
    assert m1["SPACING_CV"] > m0["SPACING_CV"]
    assert np.allclose(np.linalg.norm(lat.centers, axis=1), 1.0)


def test_zero_iterations_and_single_centre():
    s = SphereSurface(1.0)
    lat, m = build_lattice(s, Config(target_count=1, relax_iterations=0))
    assert len(lat) == 1 and np.isnan(m["SPACING_CV"])


def test_invalid_config_rejected_before_work():
    with pytest.raises(ValueError):
        build_lattice(SphereSurface(1.0), Config(pokeyness=2.0))


# ---------- regularity first: perfect patch, defects in the seams ----------

def _perfect_radius(surface, config, origin):
    """Distance (in spacings, straight line) from origin to the nearest non-hexagonal
    interior cell."""
    from ommatidia.topology.topology import build_topology, rim_cells

    lattice, _ = build_lattice(surface, config)
    topology, _ = build_topology(lattice, surface, config)
    bad = (topology.valence != 6) & ~rim_cells(topology.triangles, len(lattice))
    return np.linalg.norm(lattice.centers[bad] - origin, axis=1).min() / lattice.spacing


def test_sphere_is_perfect_around_the_front():
    """Default centre on a closed mesh faces -Y. 750 on a unit sphere: perfect to
    ~62° (measured), i.e. ~8 spacings."""
    r = _perfect_radius(SphereSurface(1.0), Config(), np.array([0.0, -1.0, 0.0]))
    assert r > 6.5


def test_cap_is_perfect_up_to_the_rim():
    """An eye-like open cap: defects only within ~2.5 rows of the rim (at 60°)."""
    cap = CapSurface(1.0, np.pi / 3)
    lattice, _ = build_lattice(cap, Config())
    r = _perfect_radius(cap, Config(), np.array([0.0, 0.0, 1.0]))
    rim_radius = 2 * np.sin(np.pi / 6) / lattice.spacing   # chord to the rim, in spacings
    assert r > rim_radius - 3.0


def test_lattice_origin_moves_the_perfect_patch():
    origin = np.array([1.0, 0.0, 0.0])
    r = _perfect_radius(SphereSurface(1.0), Config(lattice_origin=tuple(origin)), origin)
    assert r > 6.5


def test_rotation_turns_the_rows():
    s = SphereSurface(1.0)
    a, _ = build_lattice(s, Config(target_count=300))
    b, _ = build_lattice(s, Config(target_count=300, lattice_rotation_deg=30.0))
    front = np.array([0.0, -1.0, 0.0])

    def first_ring(lat):
        return lat.centers[np.argsort(np.linalg.norm(lat.centers - front, axis=1))[1:7]]

    # the first ring turned by 30°: no first-ring centre stays where it was
    d = np.linalg.norm(first_ring(a)[:, None] - first_ring(b)[None], axis=2).min(axis=1)
    assert d.min() > 0.3 * a.spacing


class _TwoSpheres:
    """Two disjoint unit spheres: two eyes in one object."""

    def __init__(self):
        self.a, self.offset = SphereSurface(1.0), np.array([3.0, 0.0, 0.0])

    def closest_point(self, p):
        p = np.asarray(p, float).reshape(-1, 3)
        x1, n1 = self.a.closest_point(p)
        x2, n2 = self.a.closest_point(p - self.offset)
        x2 = x2 + self.offset
        second = np.linalg.norm(p - x2, axis=1) < np.linalg.norm(p - x1, axis=1)
        return np.where(second[:, None], x2, x1), np.where(second[:, None], n2, n1)

    def sample_area_weighted(self, n, rng):
        x = self.a.sample_area_weighted(n, rng)
        x[rng.random(n) < 0.5] += self.offset
        return x

    def boundary_distance(self, p):
        return np.full(len(p), np.inf)

    def area(self):
        return 2 * self.a.area()


def test_second_component_gets_its_own_patch():
    lattice, metrics = build_lattice(_TwoSpheres(), Config(target_count=600))
    assert metrics["N_PATCHES"] == 2
    on_second = lattice.centers[:, 0] > 1.5
    assert 0.4 < on_second.mean() < 0.6
    # each sphere is one patch (which one comes first depends on the auto centre)
    for side in (on_second, ~on_second):
        ids = set(np.unique(lattice.patch[side]).tolist()) - {-1}
        assert len(ids) == 1
