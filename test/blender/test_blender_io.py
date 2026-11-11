"""blender_io/, pipeline and ui/: the BVH SURFACE, units, the lattice and ommatidia
objects, settings, live updates, operator.

Run with test/run_blender_tests.py inside Blender.
"""

import json
import math
import time
from dataclasses import fields

import bmesh
import bpy
import numpy as np

import ommatidia
from ommatidia.blender_io.atomic_algorithm_scripts.config_props import (
    DERIVED_FIELDS,
    config_to_props,
    props_to_config,
)
from ommatidia.blender_io.atomic_algorithm_scripts.facets_to_object import facets_object_name
from ommatidia.blender_io.blender_io import MODIFIER_NAME, surface_from_object, write_lattice_object
from ommatidia.blender_io.atomic_algorithm_scripts.lattice_to_object import lattice_object_name
from ommatidia.config import Config
from ommatidia.lattice.lattice import build_lattice
from ommatidia.pipeline import generate as run_pipeline
from ommatidia.pipeline import reshape
from ommatidia.topology.topology import build_topology
from ommatidia.ui import live


def approx(a, b, rel=1e-6, abs_=0.0):
    return abs(a - b) <= max(rel * abs(b), abs_)


def ico(radius=1.0, subdivisions=5, location=(0, 0, 0)):
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=subdivisions, radius=radius, location=location)
    return bpy.context.active_object


def open_cap(radius=1.0, z_cut=0.3):
    """An open sheet: an icosphere with everything below z_cut deleted."""
    obj = ico(radius, 6)
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bmesh.ops.delete(bm, geom=[v for v in bm.verts if v.co.z < z_cut * radius], context="VERTS")
    bm.to_mesh(obj.data)
    bm.free()
    return obj


def generate(obj, config):
    surface = surface_from_object(obj, config)
    lattice, metrics = build_lattice(surface, config)
    topology, tm = build_topology(lattice, surface, config)
    metrics.update(tm)
    return surface, lattice, topology, metrics, write_lattice_object(obj, lattice, topology, config, metrics)


# ---------- SURFACE on a real mesh ----------

def test_bvh_surface_sphere_invariants():
    obj = ico(2.0, 5)
    bpy.ops.object.shade_smooth()      # primitives are flat-shaded; smooth gives radial normals
    surface = surface_from_object(obj, Config())
    assert approx(surface.area(), 4 * math.pi * 4, rel=0.01), surface.area()
    rng = np.random.default_rng(0)
    q = rng.normal(size=(300, 3)) * 3
    x, n = surface.closest_point(q)
    x2, _ = surface.closest_point(x)
    assert np.allclose(x, x2, atol=1e-6)
    assert np.allclose(np.linalg.norm(n, axis=1), 1.0)
    radial = x / np.linalg.norm(x, axis=1, keepdims=True)
    assert np.min(np.sum(n * radial, axis=1)) > 0.999   # smooth, outward
    s = surface.sample_area_weighted(500, rng)
    assert np.allclose(surface.closest_point(s)[0], s, atol=1e-6)
    assert np.all(np.isinf(surface.boundary_distance(q)))


def test_model_scale_converts_to_specimen_metres():
    obj = ico(2.0, 4)
    a1 = surface_from_object(obj, Config(model_scale=1.0)).area()
    a2 = surface_from_object(obj, Config(model_scale=1000.0)).area()
    assert approx(a2, a1 / 1e6, rel=1e-9)


def test_world_transform_and_normals_under_nonuniform_scale():
    obj = ico(1.0, 5, location=(5, 0, 0))
    obj.scale = (1.0, 1.0, 3.0)
    bpy.context.view_layer.update()
    surface = surface_from_object(obj, Config())
    x, n = surface.closest_point(surface.sample_area_weighted(300, np.random.default_rng(1)))
    local = x - np.array([5.0, 0, 0])
    assert np.allclose(np.sum((local / [1, 1, 3]) ** 2, axis=1), 1.0, atol=0.01)
    gradient = local / np.array([1.0, 1.0, 9.0])     # ellipsoid normal direction
    gradient /= np.linalg.norm(gradient, axis=1, keepdims=True)
    assert np.min(np.sum(n * gradient, axis=1)) > 0.995


def test_open_sheet_boundary_distance():
    bpy.ops.mesh.primitive_grid_add(x_subdivisions=20, y_subdivisions=20, size=2.0)
    surface = surface_from_object(bpy.context.active_object, Config())
    d = surface.boundary_distance(np.array([[0, 0, 0], [0.9, 0, 0], [0.5, 0.5, 0], [-1, -1, 0]]))
    assert np.allclose(d, [1.0, 0.1, 0.5, 0.0], atol=1e-9), d


def test_own_modifier_is_ignored_and_restored():
    obj = ico(1.0, 4)
    area = surface_from_object(obj, Config()).area()
    mod = obj.modifiers.new(MODIFIER_NAME, "DISPLACE")
    mod.strength = 0.5
    assert approx(surface_from_object(obj, Config()).area(), area, rel=1e-9)
    assert mod.show_viewport
    other = obj.modifiers.new("UserDisplace", "DISPLACE")   # user modifiers do count
    other.strength = 0.5
    assert surface_from_object(obj, Config()).area() > 1.5 * area


def test_mesh_without_faces_is_rejected():
    mesh = bpy.data.meshes.new("pts")
    mesh.from_pydata([(0, 0, 0), (1, 0, 0)], [], [])
    obj = bpy.data.objects.new("pts", mesh)
    bpy.context.scene.collection.objects.link(obj)
    try:
        surface_from_object(obj, Config())
    except ValueError as error:
        assert "no faces" in str(error)
    else:
        raise AssertionError("expected ValueError")


# ---------- end to end ----------

def test_closed_mesh_end_to_end():
    obj = ico(1.0, 6)
    config = Config(model_scale=1000.0)  # COUNT mode: model_scale must cancel out
    _, lattice, topology, m, lat_obj = generate(obj, config)
    assert 0.93 <= m["N"] / 750 <= 1.07 and m["HOLES_ON_CLOSED_SURFACE"] == 0
    assert m["SPACING_CV"] < 0.07 and m["HEX_FRACTION"] >= 0.92, m
    assert 0.9 < m["N_BRISTLES"] / m["N"] < 1.1
    mesh = lat_obj.data
    assert len(mesh.vertices) == m["N"] + m["N_TRIANGLES"]
    assert len(mesh.edges) == len(topology.cell_edges)
    kinds = np.empty(len(mesh.vertices), dtype=np.int32)
    mesh.attributes["omm_kind"].data.foreach_get("value", kinds)
    assert (kinds == 0).sum() == m["N"]
    bristles = np.empty(len(mesh.vertices), dtype=bool)
    mesh.attributes["omm_bristle"].data.foreach_get("value", bristles)
    assert bristles.sum() == m["N_BRISTLES"] and not bristles[kinds == 0].any()
    for name in ("omm_normal", "omm_valence"):
        assert name in mesh.attributes
    # Positions are back in model units: on the unit icosphere (model_scale cancels).
    co = np.empty(len(mesh.vertices) * 3)
    mesh.vertices.foreach_get("co", co)
    assert np.allclose(np.linalg.norm(co.reshape(-1, 3)[kinds == 0], axis=1), 1.0, atol=0.01)
    assert Config.from_json(lat_obj["omm_config"]) == config
    assert json.loads(lat_obj["omm_metrics"])["N"] == m["N"]
    assert lat_obj["omm_target"] == obj.name


def test_regenerating_replaces_the_lattice_object():
    obj = ico(1.0, 5)
    generate(obj, Config(target_count=200, relax_iterations=20))
    meshes_before = len(bpy.data.meshes)
    _, _, _, _, lat_obj = generate(obj, Config(target_count=300, relax_iterations=20))
    assert [o.name for o in bpy.data.objects].count(lattice_object_name(obj.name)) == 1
    assert len(bpy.data.meshes) == meshes_before        # old mesh removed, not leaked
    assert json.loads(lat_obj["omm_metrics"])["N"] == m_n(lat_obj)


def m_n(lat_obj):
    return json.loads(lat_obj["omm_metrics"])["N"]


def test_open_sheet_end_to_end():
    obj = open_cap(1.0, 0.3)
    surface, lattice, _, m, _ = generate(obj, Config(target_count=500))
    assert m["N"] < m["N_TARGET"] == 500                 # rim half-cells trimmed
    assert np.all(surface.boundary_distance(lattice.centers) >= 0.5 * lattice.spacing)
    assert m["N_OPEN_EDGES"] > 0 and m["HEX_FRACTION"] >= 0.92, m


def test_determinism_on_mesh():
    obj = ico(1.0, 5)
    cfg = Config(target_count=200, relax_iterations=30)
    a = build_lattice(surface_from_object(obj, cfg), cfg)[0].centers
    b = build_lattice(surface_from_object(obj, cfg), cfg)[0].centers
    assert np.array_equal(a, b)


# ---------- UI: settings and operator ----------

def _registered():
    try:
        ommatidia.register()
    except ValueError:  # already registered in this Blender session
        pass


def test_settings_cover_every_config_field_and_round_trip():
    _registered()
    obj = ico(1.0, 3)
    props = obj.ommatidia
    names = {f.name for f in fields(Config)} - set(DERIVED_FIELDS)
    assert names <= set(props.bl_rna.properties.keys())
    assert props_to_config(props) == Config() or _close(props_to_config(props), Config())
    wanted = Config(density_mode="SPACING", spacing=20e-6, seed=5, bristle_pattern="RANDOM",
                    comb_direction=(1.0, 0.0, 0.0), bristle_radius=0.7e-6, model_scale=250.0,
                    lattice_origin=(0.001, 0.0, -0.002), lattice_rotation_deg=15.0, inflation=0.8)
    config_to_props(wanted, props)
    assert _close(props_to_config(props), wanted)
    assert abs(props.spacing - 20.0) < 1e-5            # the UI shows µm


def _close(a, b):
    """Equal up to float32 storage of Blender float properties."""
    for f in fields(Config):
        x, y = getattr(a, f.name), getattr(b, f.name)
        if isinstance(y, float):
            if not approx(x, y, rel=1e-6, abs_=1e-12):
                return False
        elif isinstance(y, tuple):
            if x is None or not np.allclose(x, y, rtol=1e-6):
                return False
        elif x != y:
            return False
    return True


def test_operator_generates_and_reports():
    _registered()
    obj = ico(1.0, 5)
    obj.ommatidia.target_count = 300
    obj.ommatidia.relax_iterations = 50
    assert bpy.ops.ommatidia.generate() == {"FINISHED"}
    lat_obj = bpy.data.objects[lattice_object_name(obj.name)]
    assert 0.9 < json.loads(lat_obj["omm_metrics"])["N"] / 300 < 1.1
    assert bpy.data.objects.get(facets_object_name(obj.name)) is not None
    bpy.context.view_layer.objects.active = lat_obj
    assert not bpy.ops.ommatidia.generate.poll()         # a lattice object is not a target


def test_operator_reports_invalid_density_instead_of_freezing():
    _registered()
    obj = ico(1.0, 3)    # 1 m sphere at model_scale 1 = a 1 m "eye": far too many 16 µm facets
    obj.ommatidia.density_mode = "SPACING"
    try:
        result = bpy.ops.ommatidia.generate()
    except RuntimeError as error:   # operator ERROR reports raise in background mode
        assert "max_count" in str(error)
    else:
        assert result == {"CANCELLED"}


# ---------- performance (logged, never asserted; AGENTS §8) ----------

def test_performance_800_on_200k_vertex_mesh():
    obj = ico(1.0, 8)                       # 163,842 vertices (Blender level 1 = icosahedron)
    obj.scale = (1.0, 0.8, 0.6)
    bpy.context.view_layer.update()
    config = Config(target_count=800)
    t0 = time.perf_counter()
    surface = surface_from_object(obj, config)
    t1 = time.perf_counter()
    lattice, m = build_lattice(surface, config)
    t2 = time.perf_counter()
    _, tm = build_topology(lattice, surface, config)
    t3 = time.perf_counter()
    s = m["SECONDS"]
    print(f"  PERF {len(obj.data.vertices)} verts: surface {t1 - t0:.2f}s, lattice {t2 - t1:.2f}s "
          f"(sampling {s['sampling']:.2f}, growth {s['growth']:.2f}, seams {s['seams']:.2f}), "
          f"topology {t3 - t2:.2f}s; hex {tm['HEX_FRACTION']:.3f}, cv {m['SPACING_CV']:.3f}")


# ---------- shading and sharp creases ----------

def test_flat_faces_get_face_normals():
    """Corner normals follow the user's shading: the default cube is flat, so the
    normal anywhere on its top face is +Z, not an interpolated corner diagonal."""
    bpy.ops.mesh.primitive_cube_add(size=2.0)
    surface = surface_from_object(bpy.context.active_object, Config())
    _, n = surface.closest_point(np.array([[0.3, 0.2, 3.0], [-0.9, 0.9, 2.0]]))
    assert np.allclose(n, [[0, 0, 1], [0, 0, 1]], atol=1e-6)


def test_smoothed_cube_is_watertight_and_hard_cube_reports_creases():
    bpy.ops.mesh.primitive_cube_add(size=2.0)
    cube = bpy.context.active_object
    cube.modifiers.new("Smooth", "SUBSURF").levels = 3
    _, _, _, m, _ = generate(cube, Config(target_count=600))
    assert m["HOLES_ON_CLOSED_SURFACE"] == 0 and m["N_OVERLAP_EDGES"] == 0, m
    assert m["HEX_FRACTION"] >= 0.88
    cube.modifiers.clear()
    _, _, _, m, _ = generate(cube, Config(target_count=600))
    print(f"  HARD CUBE: hex {m['HEX_FRACTION']:.3f}, holes {m['HOLES_ON_CLOSED_SURFACE']}, "
          f"overlaps {m['N_OVERLAP_EDGES']} (known limitation at 90° creases)")
    assert m["HEX_FRACTION"] >= 0.7


# ---------- the visible ommatidia object ----------

def _facet_faces(obj):
    mesh = obj.data
    kind = np.empty(len(mesh.polygons), dtype=np.int32)
    mesh.polygons.foreach_get("material_index", kind)
    return [tuple(p.vertices) for p, k in zip(mesh.polygons, kind) if k == 0]


def test_pipeline_builds_a_watertight_parented_ommatidia_object():
    obj = ico(1.0, 6)
    bpy.ops.object.shade_smooth()
    m = run_pipeline(obj, Config(target_count=400))
    facets = bpy.data.objects[facets_object_name(obj.name)]
    assert facets.parent == obj and "omm_target" in facets
    assert len(facets.material_slots) == 2 and facets.material_slots[0].material is not None
    assert facets.data.validate() is False                # nothing to fix
    edges = {}
    for face in _facet_faces(facets):
        for a, b in zip(face, face[1:] + face[:1]):
            key = (min(a, b), max(a, b))
            edges[key] = edges.get(key, 0) + 1
    assert set(edges.values()) == {2}                     # closed: every wall shared by 2 faces
    co = np.empty(len(facets.data.vertices) * 3)
    facets.data.vertices.foreach_get("co", co)
    r = np.linalg.norm(co.reshape(-1, 3), axis=1)
    assert r.min() > 0.99 and r.max() < 1.2               # domes stand on the sphere
    assert m["N_FACET_VERTICES"] == len(facets.data.vertices)
    lattice = bpy.data.objects[lattice_object_name(obj.name)]
    assert lattice.parent == obj and lattice.hide_get()


def test_ommatidia_follow_the_target_when_it_moves():
    obj = ico(1.0, 5)
    run_pipeline(obj, Config(target_count=200))
    facets = bpy.data.objects[facets_object_name(obj.name)]
    before = facets.matrix_world.copy()
    obj.location = (3.0, 0.0, 0.0)
    bpy.context.view_layer.update()
    assert np.allclose(np.array(facets.matrix_world.translation) - np.array(before.translation), [3, 0, 0])


def test_reshape_is_fast_and_moves_only_vertices():
    obj = ico(1.0, 6)
    run_pipeline(obj, Config(target_count=750))
    facets = bpy.data.objects[facets_object_name(obj.name)]
    n_faces = len(facets.data.polygons)

    def heights():
        co = np.empty(len(facets.data.vertices) * 3)
        facets.data.vertices.foreach_get("co", co)
        return co.reshape(-1, 3)

    low = heights()
    t0 = time.perf_counter()
    assert reshape(obj, Config(target_count=750, inflation=1.0, bristle_length=20e-6))
    seconds = time.perf_counter() - t0
    high = heights()
    print(f"  RESHAPE {len(low)} vertices in {seconds * 1000:.0f} ms")
    assert seconds < 1.0 and len(facets.data.polygons) == n_faces
    assert np.abs(high - low).max() > 0
    assert json.loads(bpy.data.objects[lattice_object_name(obj.name)]["omm_config"])["inflation"] == 1.0


def test_materials_survive_regeneration():
    obj = ico(1.0, 5)
    run_pipeline(obj, Config(target_count=200))
    facets = bpy.data.objects[facets_object_name(obj.name)]
    mine = bpy.data.materials.new("MyCornea")
    facets.material_slots[0].material = mine
    run_pipeline(obj, Config(target_count=250))
    assert facets.material_slots[0].material == mine


def test_bristles_none_has_no_bristle_faces():
    obj = ico(1.0, 5)
    run_pipeline(obj, Config(target_count=200, bristle_pattern="NONE"))
    facets = bpy.data.objects[facets_object_name(obj.name)]
    kind = np.empty(len(facets.data.polygons), dtype=np.int32)
    facets.data.polygons.foreach_get("material_index", kind)
    assert not (kind == 1).any()


# ---------- live updates ----------

def test_shape_setting_updates_at_once():
    _registered()
    obj = ico(1.0, 5)
    obj.ommatidia.target_count = 200
    assert bpy.ops.ommatidia.generate() == {"FINISHED"}
    facets = bpy.data.objects[facets_object_name(obj.name)]
    co0 = np.empty(len(facets.data.vertices) * 3)
    facets.data.vertices.foreach_get("co", co0)
    obj.ommatidia.inflation = 0.9                     # update callback reshapes now
    co1 = np.empty_like(co0)
    facets.data.vertices.foreach_get("co", co1)
    assert np.abs(co1 - co0).max() > 0


def test_lattice_setting_regenerates_after_a_pause():
    _registered()
    obj = ico(1.0, 5)
    obj.ommatidia.target_count = 200
    obj.ommatidia.seed = 0
    obj.ommatidia.seed = 0                      # before the first Generate: nothing scheduled
    assert not bpy.app.timers.is_registered(live._regenerate_pending)
    assert bpy.ops.ommatidia.generate() == {"FINISHED"}
    obj.ommatidia.target_count = 300
    obj.ommatidia.target_count = 320            # a second change restarts the pause
    assert bpy.app.timers.is_registered(live._regenerate_pending)
    live.cancel_pending()                       # timers do not fire in background mode
    live.lattice_changed(obj)
    live._regenerate_pending()                  # ... so run the scheduled work directly
    lat = bpy.data.objects[lattice_object_name(obj.name)]
    assert 0.9 < json.loads(lat["omm_metrics"])["N"] / 320 < 1.1
    live.cancel_pending()


def test_live_update_off_does_not_schedule():
    _registered()
    obj = ico(1.0, 4)
    obj.ommatidia.target_count = 100
    bpy.ops.ommatidia.generate()
    obj.ommatidia.live_update = False
    obj.ommatidia.seed = 3
    assert not bpy.app.timers.is_registered(live._regenerate_pending)


def test_cursor_centre_moves_the_regular_patch():
    _registered()
    obj = ico(1.0, 6)
    bpy.ops.object.shade_smooth()
    bpy.context.scene.cursor.location = (0.0, 0.0, 2.0)
    obj.ommatidia.origin_source = "CURSOR"
    config = props_to_config(obj.ommatidia, bpy.context.scene)
    assert np.allclose(config.lattice_origin, (0, 0, 2))
    surface, lattice, topology, _, _ = generate(obj, config)
    hexagonal = topology.valence == 6
    top = lattice.centers[:, 2] > 0.8                  # within ~37° of the cursor
    assert hexagonal[top].all()


def test_invalid_setting_shows_status_instead_of_raising():
    _registered()
    obj = ico(1.0, 4)
    obj.ommatidia.target_count = 100
    bpy.ops.ommatidia.generate()
    obj.ommatidia.comb_direction = (0.0, 0.0, 0.0)     # invalid; reshape reports it
    assert "comb_direction" in obj.ommatidia.status


# ---------- thin solids (Sculpt > Mask Extract with "Extract as Solid") ----------

def solid_cap(thickness):
    obj = open_cap(1.0, 0.3)
    bpy.ops.object.shade_smooth()
    mod = obj.modifiers.new("Solid", "SOLIDIFY")
    mod.thickness, mod.offset = thickness, -1.0
    bpy.ops.object.modifier_apply(modifier="Solid")
    return obj


def test_thin_solid_is_covered_on_its_outer_skin_only():
    for thickness in (0.002, 0.02):
        obj = solid_cap(thickness)
        m = run_pipeline(obj, Config(target_count=600))
        assert m["DROPPED_INNER_FACES"] > 0
        lattice_obj = bpy.data.objects[lattice_object_name(obj.name)]
        mesh = lattice_obj.data
        kinds = np.empty(len(mesh.vertices), dtype=np.int32)
        mesh.attributes["omm_kind"].data.foreach_get("value", kinds)
        normals = np.empty(len(mesh.vertices) * 3, dtype=np.float32)
        mesh.attributes["omm_normal"].data.foreach_get("vector", normals)
        co = np.empty(len(mesh.vertices) * 3)
        mesh.vertices.foreach_get("co", co)
        c, n = co.reshape(-1, 3)[kinds == 0], normals.reshape(-1, 3)[kinds == 0]
        outward = np.sum(c / np.linalg.norm(c, axis=1, keepdims=True) * n, axis=1)
        assert (outward > 0.5).mean() > 0.97, thickness          # on the outer skin
        assert m["N_PATCHES"] == 1 and m["HEX_FRACTION"] >= 0.93, m
        assert 0.85 < m["N"] / 600 < 1.0                          # count is not split over two skins


def test_closed_and_open_meshes_are_not_thin_solids():
    assert surface_from_object(ico(1.0, 5), Config()).dropped_inner_faces == 0
    assert surface_from_object(open_cap(1.0, 0.3), Config()).dropped_inner_faces == 0
