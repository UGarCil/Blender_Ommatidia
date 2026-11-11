"""Analytic SURFACE implementations used as test fixtures.

Each class follows the CDD. SURFACE protocol (docs/agents/data_definitions.md §3):

    surface.closest_point(points) -> (points, normals)
    surface.sample_area_weighted(n, rng) -> points
    surface.boundary_distance(points) -> distances   (inf on closed surfaces)
    surface.area() -> float

They are pure numpy, so the lattice can be tested without Blender.
"""

from __future__ import annotations

import numpy as np


def _unit_rows(v: np.ndarray) -> np.ndarray:
    """Normalise rows; a zero row becomes +Z (the projection of the centre is arbitrary)."""
    norm = np.linalg.norm(v, axis=1, keepdims=True)
    out = np.where(norm > 0, v / np.where(norm > 0, norm, 1.0), 0.0)
    out[norm[:, 0] == 0] = (0.0, 0.0, 1.0)
    return out


def _random_directions(n: int, rng: np.random.Generator) -> np.ndarray:
    return _unit_rows(rng.normal(size=(n, 3)))


# CDD. SURFACE (analytic): sphere of radius r centred at the origin. Closed.
class SphereSurface:
    def __init__(self, radius: float = 1.0):
        self.radius = float(radius)

    def closest_point(self, points):
        normals = _unit_rows(np.asarray(points, dtype=float).reshape(-1, 3))
        return normals * self.radius, normals

    def sample_area_weighted(self, n, rng):
        return _random_directions(n, rng) * self.radius

    def boundary_distance(self, points):
        return np.full(len(points), np.inf)

    def area(self):
        return 4.0 * np.pi * self.radius**2


# CDD. SURFACE (analytic): ellipsoid x²/a² + y²/b² + z²/c² = 1. Closed.
class EllipsoidSurface:
    def __init__(self, a: float = 1.0, b: float = 0.8, c: float = 0.6):
        self.axes = np.array([a, b, c], dtype=float)
        self._area = None

    def normals_at(self, surface_points):
        return _unit_rows(surface_points / self.axes**2)

    def closest_point(self, points):
        """Project with the Lagrange parameter t: x_i = p_i a_i² / (a_i² + t),
        where t solves sum((a_i p_i / (a_i² + t))²) = 1. f(t) decreases on
        (-min a², inf), so a vectorised bisection finds it."""
        p = np.asarray(points, dtype=float).reshape(-1, 3)
        a2 = self.axes**2

        def f(t):
            return np.sum((self.axes * p / (a2 + t[:, None])) ** 2, axis=1) - 1.0

        lo = np.full(len(p), -a2.min() * (1 - 1e-12))
        hi = np.full(len(p), max(np.linalg.norm(p, axis=1).max(initial=0.0), 1.0) * self.axes.max())
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            go_right = f(mid) > 0
            lo = np.where(go_right, mid, lo)
            hi = np.where(go_right, hi, mid)
        t = 0.5 * (lo + hi)
        x = p * a2 / (a2 + t[:, None])
        # Points on the minor axis inside the body have no unique projection;
        # rescale onto the surface so the invariant holds anyway.
        x /= np.sqrt(np.sum((x / self.axes) ** 2, axis=1, keepdims=True))
        return x, self.normals_at(x)

    def _area_density(self, u):
        # Area element of the map u (unit sphere) -> axes * u, divided by abc.
        return np.sqrt(np.sum((u / self.axes) ** 2, axis=1))

    def sample_area_weighted(self, n, rng):
        g_max = 1.0 / self.axes.min()
        out = []
        count = 0
        while count < n:
            u = _random_directions(2 * (n - count) + 16, rng)
            keep = rng.random(len(u)) * g_max < self._area_density(u)
            out.append(u[keep])
            count += int(keep.sum())
        return np.concatenate(out)[:n] * self.axes

    def boundary_distance(self, points):
        return np.full(len(points), np.inf)

    def area(self):
        if self._area is None:
            u = _random_directions(2_000_000, np.random.default_rng(0))
            self._area = float(4 * np.pi * np.prod(self.axes) * self._area_density(u).mean())
        return self._area


# CDD. SURFACE (analytic): spherical cap, polar angle theta <= theta_max. Open.
class CapSurface:
    def __init__(self, radius: float = 1.0, theta_max: float = np.pi / 3):
        self.radius = float(radius)
        self.theta_max = float(theta_max)

    def closest_point(self, points):
        p = np.asarray(points, dtype=float).reshape(-1, 3)
        n = _unit_rows(p)
        theta = np.arccos(np.clip(n[:, 2], -1, 1))
        outside = theta > self.theta_max
        if outside.any():
            phi = np.arctan2(n[outside, 1], n[outside, 0])
            st, ct = np.sin(self.theta_max), np.cos(self.theta_max)
            n[outside] = np.column_stack([st * np.cos(phi), st * np.sin(phi), np.full(len(phi), ct)])
        return n * self.radius, n

    def sample_area_weighted(self, n, rng):
        z = rng.uniform(np.cos(self.theta_max), 1.0, n)  # Archimedes: uniform in z
        phi = rng.uniform(0, 2 * np.pi, n)
        s = np.sqrt(1 - z**2)
        return np.column_stack([s * np.cos(phi), s * np.sin(phi), z]) * self.radius

    def boundary_distance(self, points):
        _, n = self.closest_point(points)
        theta = np.arccos(np.clip(n[:, 2], -1, 1))
        return np.maximum(self.theta_max - theta, 0.0) * self.radius

    def area(self):
        return 2 * np.pi * self.radius**2 * (1 - np.cos(self.theta_max))


surface0 = SphereSurface(1.0)  # stub
surface1 = EllipsoidSurface(1.0, 0.8, 0.6)
surface2 = CapSurface(1.0, np.pi / 3)
