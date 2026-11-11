"""Ommatidia, headless: cover a mesh in a .blend with ommatidia and save the result.

    blender -b input.blend --python main.py -- path/to/config.json

The config path is optional (defaults: Drosophila, see src/ommatidia/config.py).
Each run writes 00_results/<time>_seed<seed>/ with scene.blend, config.json
and metrics.json. See docs/user/README.md.
"""

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

import bpy  # noqa: E402  (only available inside Blender)

from ommatidia.blender_io.atomic_algorithm_scripts.facets_to_object import write_facets_object  # noqa: E402
from ommatidia.blender_io.blender_io import surface_from_object, write_lattice_object  # noqa: E402
from ommatidia.config import Config  # noqa: E402
from ommatidia.facets.facets import build_facet_mesh  # noqa: E402
from ommatidia.lattice.lattice import build_lattice  # noqa: E402
from ommatidia.material.material import assign_materials  # noqa: E402
from ommatidia.topology.topology import build_topology  # noqa: E402


# FD. config_path_from_argv
# Signature: None -> path | None
def config_path_from_argv():
    """The single positional argument after Blender's '--', if any."""
    args = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    return Path(args[0]) if args else None


# FD. pick_target
# Signature: str | None -> bpy Object
def pick_target(name):
    """The named object, else the active mesh, else the only mesh in the file."""
    if name:
        return bpy.data.objects[name]
    active = bpy.context.view_layer.objects.active
    if active is not None and active.type == "MESH":
        return active
    meshes = [o for o in bpy.context.scene.objects if o.type == "MESH" and "omm_target" not in o]
    if len(meshes) != 1:
        raise SystemExit(f"Set TARGET_NAME in main.py: found {len(meshes)} mesh objects")
    return meshes[0]


if __name__ == "__main__":
    TARGET_NAME = None                    # None = active mesh, or the only mesh
    RESULTS_ROOT = ROOT / "00_results"

    # 1. configuration
    config_path = config_path_from_argv()
    config = Config.load(config_path) if config_path else Config()
    config.validate()

    # 2. the surface to cover (specimen metres)
    target = pick_target(TARGET_NAME)
    surface = surface_from_object(target, config)

    # 3. ommatidium centres (a regular patch, seams filled), then facet outlines and bristle sites
    lattice, metrics = build_lattice(surface, config)
    topology, topology_metrics = build_topology(lattice, surface, config)
    metrics.update(topology_metrics)
    metrics["DROPPED_INNER_FACES"] = surface.dropped_inner_faces   # thin solid: inner skin

    # 4. the visible ommatidia: domed facets and bristles, with the SEM material
    facet_mesh = build_facet_mesh(lattice, topology, surface, config)
    ommatidia_object = write_facets_object(target, facet_mesh, config)
    assign_materials(ommatidia_object)

    # 5. the lattice itself (hidden; holds the config and metrics)
    lattice_object = write_lattice_object(target, lattice, topology, config, metrics)

    # 6. results: scene, config and metrics side by side
    run_dir = RESULTS_ROOT / f"{time.strftime('%Y%m%d-%H%M%S')}_seed{config.seed}"
    run_dir.mkdir(parents=True, exist_ok=True)
    config.save(run_dir / "config.json")
    metrics.update({"TARGET": target.name, "LATTICE_OBJECT": lattice_object.name,
                    "OMMATIDIA_OBJECT": ommatidia_object.name, "BLENDER": bpy.app.version_string})
    (run_dir / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(run_dir / "scene.blend"), copy=True)
    print(f"[ommatidia] {metrics['N']} ommatidia, {metrics['HEX_FRACTION']:.1%} hexagonal "
          f"-> {run_dir}")
