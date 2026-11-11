"""The analytic fixtures must obey the SURFACE invariants before anything is tested on them."""

import numpy as np
import pytest

from fixtures.analytic_surfaces import CapSurface, EllipsoidSurface, SphereSurface

SURFACES = [SphereSurface(1.0), EllipsoidSurface(1.0, 0.8, 0.6), CapSurface(1.0, np.pi / 3)]
IDS = ["sphere", "ellipsoid", "cap"]


@pytest.mark.parametrize("surface", SURFACES, ids=IDS)
def test_projection_is_idempotent_with_unit_normals(surface):
    rng = np.random.default_rng(1)
    p = rng.normal(size=(500, 3))
    x, n = surface.closest_point(p)
    x2, n2 = surface.closest_point(x)
    assert np.allclose(x, x2, atol=1e-9)
    assert np.allclose(np.linalg.norm(n, axis=1), 1.0)


@pytest.mark.parametrize("surface", SURFACES, ids=IDS)
def test_projection_is_closest(surface):
    """No sampled surface point is closer to a query than its projection."""
    rng = np.random.default_rng(2)
    q = rng.normal(size=(50, 3)) * 0.3 + rng.normal(size=(50, 3))
    x, _ = surface.closest_point(q)
    s = surface.sample_area_weighted(20_000, rng)
    d_proj = np.linalg.norm(q - x, axis=1)
    d_brute = np.linalg.norm(q[:, None, :] - s[None, :, :], axis=2).min(axis=1)
    assert np.all(d_proj <= d_brute + 1e-9)


@pytest.mark.parametrize("surface", SURFACES, ids=IDS)
def test_samples_lie_on_surface(surface):
    s = surface.sample_area_weighted(1000, np.random.default_rng(3))
    assert s.shape == (1000, 3)
    x, _ = surface.closest_point(s)
    assert np.allclose(s, x, atol=1e-9)


def test_sample_zero_points():
    assert SphereSurface().sample_area_weighted(0, np.random.default_rng(0)).shape == (0, 3)


def test_ellipsoid_area_matches_sphere_limit_and_thomsen():
    assert EllipsoidSurface(1, 1, 1).area() == pytest.approx(4 * np.pi, rel=1e-3)
    a, b, c, p = 1.0, 0.8, 0.6, 1.6075  # Knud Thomsen approximation, error < 1.1%
    thomsen = 4 * np.pi * (((a * b) ** p + (a * c) ** p + (b * c) ** p) / 3) ** (1 / p)
    assert EllipsoidSurface(a, b, c).area() == pytest.approx(thomsen, rel=0.011)


def test_ellipsoid_sampling_is_area_weighted():
    """Fraction of samples with |z| > 0.5c must equal that band's area fraction
    (estimated independently by a fine latitude-longitude quadrature)."""
    e = EllipsoidSurface(1.0, 0.8, 0.6)
    s = e.sample_area_weighted(200_000, np.random.default_rng(4))
    frac = np.mean(np.abs(s[:, 2]) > 0.3)
    th, ph = np.meshgrid(np.linspace(0, np.pi, 1201), np.linspace(0, 2 * np.pi, 1201), indexing="ij")
    def pos(t, f):
        return np.stack([np.sin(t) * np.cos(f), np.sin(t) * np.sin(f), np.cos(t)], -1) * e.axes
    dt, dp = 1e-6, 1e-6
    cross = np.cross(pos(th + dt, ph) - pos(th, ph), pos(th, ph + dp) - pos(th, ph)) / (dt * dp)
    dA = np.linalg.norm(cross, axis=-1)
    z = pos(th, ph)[..., 2]
    expected = (dA * (np.abs(z) > 0.3)).sum() / dA.sum()
    assert frac == pytest.approx(expected, abs=0.005)


def test_cap_boundary_distance():
    cap = CapSurface(2.0, np.pi / 3)
    pole = np.array([[0.0, 0.0, 5.0]])
    assert cap.boundary_distance(pole)[0] == pytest.approx(2.0 * np.pi / 3)
    rim = np.array([[np.sin(np.pi / 3), 0.0, np.cos(np.pi / 3)]]) * 2.0
    assert cap.boundary_distance(rim)[0] == pytest.approx(0.0, abs=1e-9)
    assert np.isinf(SphereSurface().boundary_distance(pole)[0])
