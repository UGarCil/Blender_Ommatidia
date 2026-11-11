"""Compound: the visible ommatidia, a composite mesh of domed facets plus bristles.

Pure numpy: never imports bpy or mathutils. See facets.md.

Two stages, so the sliders can act instantly:
- build_facet_mesh (once per Generate, needs the SURFACE): lays out every
  vertex on the surface and stores how far it is from its facet's centre.
- shape_positions (every slider change, no SURFACE): heights for the domes,
  positions for the bristles.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import Config
from .atomic_algorithm_scripts.bristles import bristle_positions, bristle_template
from .atomic_algorithm_scripts.cell_loops import cell_loops
from .atomic_algorithm_scripts.dome_mesh import facet_layout, facet_resolution
from .atomic_algorithm_scripts.dome_profile import dome_heights, rim_fade

KIND_FACET, KIND_BRISTLE = 0, 1
# Facets float this far above the input surface (x spacing), so the valleys do
# not flicker against the input mesh when both are visible.
LIFT = 0.01
NO_RIM = 1e6   # rim distance (in spacings) stored for centres on closed surfaces


# CDD. FACET_MESH
# facet_mesh = FacetMesh(base, normal, kind, rho_round, rho_hex, radius, rim,
#                        axial, angle, loop_start, loop_vertices, face_kind, spacing)
# interp. the composite before shaping, in specimen metres. Per vertex:
#   base, normal   point on the surface and its normal (for bristles: the corner site)
#   kind           KIND_FACET or KIND_BRISTLE
#   rho_round, rho_hex, radius   see FACET_LAYOUT (radius = facet inner radius, 0 on walls)
#   rim            boundary distance of the vertex's facet centre / SPACING (NO_RIM if closed)
#   axial, angle   bristle template coordinates (0 for facet vertices)
# Faces: loop_start / loop_vertices as in FACET_LAYOUT; face_kind per face.
# Invariant: shape_positions(facet_mesh, config) is defined for every valid CONFIG
@dataclass
class FacetMesh:
    base: np.ndarray
    normal: np.ndarray
    kind: np.ndarray
    rho_round: np.ndarray
    rho_hex: np.ndarray
    radius: np.ndarray
    rim: np.ndarray
    axial: np.ndarray
    angle: np.ndarray
    loop_start: np.ndarray
    loop_vertices: np.ndarray
    face_kind: np.ndarray
    spacing: float

    def __len__(self) -> int:
        return len(self.base)


# FD. build_facet_mesh
# Signature: LATTICE, TOPOLOGY, SURFACE, CONFIG -> FACET_MESH
def build_facet_mesh(lattice, topology, surface, config: Config) -> FacetMesh:
    """Lay out the facets and bristles on the surface.

    Algorithm:
        corner loop around each centre (no facet where it does not close)
        facet layout at the CONFIG resolution; project its vertices onto SURFACE
        rim distance of each facet's centre
        one bristle template per bristle site (corner with a bristle)
        concatenate facets and bristles
    """
    sides, rings = facet_resolution(config.subdivision_level)
    loops = cell_loops(topology.triangles, len(lattice))
    layout, radius = facet_layout(lattice.centers, topology.cell_points, loops,
                                  topology.cell_edges, sides, rings)
    if len(layout.points):
        base, normal = surface.closest_point(layout.points)
    else:
        base, normal = np.zeros((0, 3)), np.zeros((0, 3))
    centre_rim = surface.boundary_distance(lattice.centers) / lattice.spacing
    centre_rim = np.where(np.isfinite(centre_rim), centre_rim, NO_RIM)
    owner = layout.owner
    n_facet = len(base)
    rim = np.where(owner >= 0, centre_rim[np.maximum(owner, 0)], NO_RIM)

    sites = topology.cell_points[topology.bristle_mask]
    site_normals = topology.cell_normals[topology.bristle_mask]
    axial_t, angle_t, faces_t = bristle_template()
    k = len(axial_t)
    nb = len(sites)
    b_base = np.repeat(sites, k, axis=0).reshape(-1, 3)
    b_normal = np.repeat(site_normals, k, axis=0).reshape(-1, 3)
    offsets = n_facet + k * np.arange(nb)
    b_faces = faces_t[None] + offsets[:, None, None]
    b_faces = np.where(faces_t[None] < 0, -1, b_faces).reshape(-1, 4)
    sizes = np.where(b_faces[:, 3] < 0, 3, 4)
    b_loops = b_faces[b_faces >= 0]
    b_start = len(layout.loop_vertices) + np.concatenate([[0], np.cumsum(sizes)[:-1]]) if nb else np.zeros(0, int)

    zeros_b = np.zeros(nb * k)
    return FacetMesh(
        base=np.concatenate([base, b_base]),
        normal=np.concatenate([normal, b_normal]),
        kind=np.concatenate([np.full(n_facet, KIND_FACET), np.full(nb * k, KIND_BRISTLE)]).astype(np.int64),
        rho_round=np.concatenate([layout.rho_round, zeros_b]),
        rho_hex=np.concatenate([layout.rho_hex, zeros_b]),
        radius=np.concatenate([radius[np.maximum(owner, 0)] * (owner >= 0), zeros_b]),
        rim=np.concatenate([rim, np.full(nb * k, NO_RIM)]),
        axial=np.concatenate([np.zeros(n_facet), np.tile(axial_t, nb)]),
        angle=np.concatenate([np.zeros(n_facet), np.tile(angle_t, nb)]),
        loop_start=np.concatenate([layout.loop_start, b_start]).astype(np.int64),
        loop_vertices=np.concatenate([layout.loop_vertices, b_loops]).astype(np.int64),
        face_kind=np.concatenate([np.full(len(layout.loop_start), KIND_FACET),
                                  np.full(len(sizes) if nb else 0, KIND_BRISTLE)]).astype(np.int64),
        spacing=lattice.spacing,
    )


# FD. real_length_scale
# Signature: CONFIG, SPACING (lattice) -> float
def real_length_scale(config: Config, spacing: float) -> float:
    """Specimen metres per real metre, for lengths given in real units (bristles).

    SPACING mode: the specimen is real-size (model_scale did the conversion): 1.
    COUNT mode: the size of the specimen is unknown, so one lattice spacing is
    taken to be CONFIG.spacing real metres (16 µm by default). A unit sphere with
    750 facets then gets bristles in proportion to its facets, at any model scale.
    """
    return spacing / config.spacing if config.density_mode == "COUNT" else 1.0


# FD. shape_positions
# Signature: FACET_MESH, CONFIG -> POINTS (specimen metres)
def shape_positions(mesh: FacetMesh, config: Config) -> np.ndarray:
    """Where every vertex goes for the current shape settings (fast: no SURFACE).

    Algorithm:
        facets: base + normal · (LIFT·spacing + dome height · rim fade)
        bristles: template placed at its site, tilted toward the comb direction;
            real lengths × real_length_scale
    """
    out = np.empty_like(mesh.base)
    f = mesh.kind == KIND_FACET
    h = dome_heights(mesh.rho_round[f], mesh.rho_hex[f], mesh.radius[f], config.inflation,
                     config.pokeyness, config.facet_shape, config.crease_sharpness)
    h = LIFT * mesh.spacing + h * rim_fade(mesh.rim[f], config.rim_width)
    out[f] = mesh.base[f] + h[:, None] * mesh.normal[f]
    b = ~f
    if b.any():
        k = real_length_scale(config, mesh.spacing)
        out[b] = bristle_positions(mesh.base[b], mesh.normal[b], mesh.axial[b], mesh.angle[b],
                                   k * config.bristle_length, k * config.bristle_radius,
                                   config.bristle_tilt_deg, np.asarray(config.comb_direction))
        out[b] += LIFT * mesh.spacing * mesh.normal[b]
    return out
