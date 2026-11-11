"""Atomic: turn a Blender mesh object into a SURFACE (BVH-backed).

This is the unit boundary (data_definitions.md §1): mesh coordinates are
divided by MODEL_SCALE here, so everything downstream works in specimen metres.
"""

from __future__ import annotations

import bpy
import numpy as np
from mathutils import Vector
from mathutils.bvhtree import BVHTree
from mathutils.kdtree import KDTree

from ...topology.atomic_algorithm_scripts.triangulate import open_edges

# Boundary edges are sampled into points no farther apart than this fraction of
# sqrt(surface area); the nearest sample's segment is then measured exactly.
BOUNDARY_SAMPLE_STEP = 1e-3


# DD. MESH_ARRAYS
# mesh_arrays = (vertices (V,3) float, triangles (T,3) int, corner_normals (T,3,3) float)
# interp. a triangulated copy of an evaluated mesh, world space, in specimen
#         metres (world coordinates / MODEL_SCALE). corner_normals[t, k] is the
#         shading normal of corner k of triangle t: it follows the user's
#         smooth/flat shading (flat faces get their face normal)
# Invariant: T >= 1; normals unit length


# FD. mesh_arrays
# Signature: bpy Object, depsgraph, MODEL_SCALE, str (modifier to skip) -> MESH_ARRAYS
def mesh_arrays(obj, depsgraph, model_scale: float, skip_modifier: str = "") -> tuple:
    """Read the evaluated, triangulated mesh of obj in world space / model_scale.

    Algorithm:
        temporarily hide the named modifier (our own Ommatidia GN modifier must not
            feed the domes back into the surface they are placed on)
        evaluate the object; read vertices, loop triangles and their corner normals
        apply the world matrix (normals with its inverse transpose); divide by model_scale
        restore the modifier
    """
    if obj.type != "MESH":
        raise TypeError(f"{obj.name!r} is not a mesh object")
    modifier = obj.modifiers.get(skip_modifier) if skip_modifier else None
    was_visible = modifier.show_viewport if modifier else None
    try:
        if modifier:
            modifier.show_viewport = False
            depsgraph = bpy.context.evaluated_depsgraph_get()
            depsgraph.update()
        evaluated = obj.evaluated_get(depsgraph)
        mesh = evaluated.to_mesh()
        try:
            mesh.calc_loop_triangles()
            n_v, n_t = len(mesh.vertices), len(mesh.loop_triangles)
            if n_t == 0:
                raise ValueError(f"{obj.name!r} has no faces to cover")
            verts = np.empty(n_v * 3)
            mesh.vertices.foreach_get("co", verts)
            corner = np.empty(len(mesh.loops) * 3)
            mesh.corner_normals.foreach_get("vector", corner)
            tris = np.empty(n_t * 3, dtype=np.int64)
            mesh.loop_triangles.foreach_get("vertices", tris)
            tri_loops = np.empty(n_t * 3, dtype=np.int64)
            mesh.loop_triangles.foreach_get("loops", tri_loops)
            normals = corner.reshape(-1, 3)[tri_loops]          # (3T, 3)
        finally:
            evaluated.to_mesh_clear()
    finally:
        if modifier:
            modifier.show_viewport = was_visible
    m = np.array(obj.matrix_world)
    verts = verts.reshape(-1, 3) @ m[:3, :3].T + m[:3, 3]
    normals = normals.reshape(-1, 3) @ np.linalg.inv(m[:3, :3])  # = (M^-T n) as rows
    normals /= np.maximum(np.linalg.norm(normals, axis=1, keepdims=True), 1e-300)
    return verts / model_scale, tris.reshape(-1, 3), normals.reshape(-1, 3, 3)


def _barycentric(p: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    v0, v1, v2 = b - a, c - a, p - a
    d00 = np.sum(v0 * v0, 1); d01 = np.sum(v0 * v1, 1); d11 = np.sum(v1 * v1, 1)
    d20 = np.sum(v2 * v0, 1); d21 = np.sum(v2 * v1, 1)
    den = d00 * d11 - d01 * d01
    den = np.where(np.abs(den) > 0, den, 1.0)
    v = (d11 * d20 - d01 * d21) / den
    w = (d00 * d21 - d01 * d20) / den
    return np.column_stack([1.0 - v - w, v, w])


# A mesh counts as a thin solid (e.g. Sculpt ▸ Mask Extract with "Extract as Solid")
# when at least this share of its area has an opposite-facing layer right behind it.
THIN_SHELL_SHARE = 0.25
OPPOSITE = -0.5            # cosine: the layer behind faces the other way


# FD. inner_layer_faces
# Signature: POINTS (vertices), TRIANGLES, float (max gap) -> (T,) bool
def inner_layer_faces(vertices: np.ndarray, triangles: np.ndarray, max_gap: float) -> np.ndarray:
    """Faces of the inner skin of a thin solid; none if the mesh is not one.

    Two skins closer than about one facet cannot both carry ommatidia: the
    lattice would hop between them. The outer skin is the one facing away from
    the shape's centre (an eye bulges outward).

    Algorithm:
        for each face: cast a ray from its centre backwards (−normal) up to max_gap
        doubled = the ray hits a face that faces the other way
        if doubled faces cover < THIN_SHELL_SHARE of the area: not a thin solid
        inner = doubled faces whose normal points toward the area-weighted centroid
    """
    a, b, c = (vertices[triangles[:, k]] for k in range(3))
    cross = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(cross, axis=1)
    normal = cross / np.maximum(2 * area[:, None], 1e-300)
    centre = (a + b + c) / 3.0
    bvh = BVHTree.FromPolygons([Vector(v) for v in vertices], triangles.tolist(), all_triangles=True)
    doubled = np.zeros(len(triangles), dtype=bool)
    eps = 1e-6 * max_gap
    for f in range(len(triangles)):
        n = Vector(normal[f])
        hit, hit_normal, index, _ = bvh.ray_cast(Vector(centre[f]) - eps * n, -n, max_gap)
        if hit is not None and index != f and hit_normal.dot(n) < OPPOSITE:
            doubled[f] = True
    total = area.sum()
    if total == 0 or area[doubled].sum() < THIN_SHELL_SHARE * total:
        return np.zeros(len(triangles), dtype=bool)
    centroid = (area[:, None] * centre).sum(axis=0) / total
    return doubled & (np.sum((centre - centroid) * normal, axis=1) < 0)


# CDD. SURFACE (Blender implementation)
# surface = BvhSurface(vertices, triangles, corner_normals)
# interp. a triangle mesh in specimen metres; see data_definitions.md §3
class BvhSurface:
    def __init__(self, vertices: np.ndarray, triangles: np.ndarray, corner_normals: np.ndarray):
        self.vertices = np.asarray(vertices, dtype=float)
        self.triangles = np.asarray(triangles, dtype=np.int64)
        self.corner_normals = np.asarray(corner_normals, dtype=float).reshape(-1, 3, 3)
        self._bvh = BVHTree.FromPolygons([Vector(v) for v in self.vertices],
                                         self.triangles.tolist(), all_triangles=True)
        a, b, c = (self.vertices[self.triangles[:, k]] for k in range(3))
        cross = np.cross(b - a, c - a)
        self._areas = 0.5 * np.linalg.norm(cross, axis=1)
        self._face_normals = cross / np.maximum(2 * self._areas[:, None], 1e-300)
        self._area = float(self._areas.sum())
        self._boundary = open_edges(self.triangles)
        self._build_boundary_index()

    def _build_boundary_index(self):
        """KD-tree over points sampled along the open edges; each sample remembers its segment."""
        self._kd = None
        if len(self._boundary) == 0:
            return
        step = BOUNDARY_SAMPLE_STEP * np.sqrt(self._area)
        p0 = self.vertices[self._boundary[:, 0]]
        p1 = self.vertices[self._boundary[:, 1]]
        pieces = np.maximum(1, np.ceil(np.linalg.norm(p1 - p0, axis=1) / step)).astype(int)
        seg = np.repeat(np.arange(len(p0)), pieces + 1)
        t = np.concatenate([np.linspace(0, 1, k + 1) for k in pieces])
        samples = p0[seg] + t[:, None] * (p1[seg] - p0[seg])
        self._kd = KDTree(len(samples))
        for i, s in enumerate(samples):
            self._kd.insert(Vector(s), i)
        self._kd.balance()
        self._sample_segment = seg

    def area(self) -> float:
        return self._area

    def closest_point(self, points):
        """Projection onto the mesh, with the triangle's corner normals interpolated
        barycentrically (smooth where the mesh is smooth-shaded, flat where flat)."""
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        loc = np.empty_like(points)
        tri = np.empty(len(points), dtype=np.int64)
        for i, p in enumerate(points):
            hit, _, index, _ = self._bvh.find_nearest(Vector(p))
            loc[i] = hit
            tri[i] = index
        t = self.triangles[tri]
        w = _barycentric(loc, self.vertices[t[:, 0]], self.vertices[t[:, 1]], self.vertices[t[:, 2]])
        w = np.clip(w, 0.0, 1.0)
        n = np.einsum("ik,ikj->ij", w, self.corner_normals[tri])
        length = np.linalg.norm(n, axis=1, keepdims=True)
        n = np.where(length > 1e-12, n / np.maximum(length, 1e-300), self._face_normals[tri])
        return loc, n

    def sample_area_weighted(self, n, rng):
        if n == 0:
            return np.zeros((0, 3))
        tri = rng.choice(len(self.triangles), size=n, p=self._areas / self._area)
        r1, r2 = rng.random(n), rng.random(n)
        s = np.sqrt(r1)
        t = self.triangles[tri]
        a, b, c = (self.vertices[t[:, k]] for k in range(3))
        return (1 - s)[:, None] * a + (s * (1 - r2))[:, None] * b + (s * r2)[:, None] * c

    def boundary_distance(self, points):
        """Euclidean distance to the nearest open edge; inf on a closed mesh."""
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        if self._kd is None:
            return np.full(len(points), np.inf)
        seg = np.array([self._sample_segment[self._kd.find(Vector(p))[1]] for p in points],
                       dtype=np.int64)
        p0 = self.vertices[self._boundary[seg, 0]]
        d = self.vertices[self._boundary[seg, 1]] - p0
        t = np.clip(np.sum((points - p0) * d, 1) / np.maximum(np.sum(d * d, 1), 1e-300), 0, 1)
        return np.linalg.norm(points - (p0 + t[:, None] * d), axis=1)
