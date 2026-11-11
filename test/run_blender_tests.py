"""Run the Blender-side tests headless (Blender's Python has no pytest):

    blender -b --factory-startup --python test/run_blender_tests.py [-- name_filter]

Every function named test_* in test/blender/test_*.py runs in an empty scene.
Exit code 1 if any test fails.
"""

import importlib.util
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT / "src"), str(HERE), str(HERE / "blender")]

import bpy  # noqa: E402

name_filter = sys.argv[sys.argv.index("--") + 1] if "--" in sys.argv[:-1] else ""
failures, passed = [], 0
for path in sorted((HERE / "blender").glob("test_*.py")):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for name in [n for n in dir(module) if n.startswith("test_") and name_filter in n]:
        bpy.ops.wm.read_factory_settings(use_empty=True)
        t0 = time.perf_counter()
        try:
            getattr(module, name)()
            passed += 1
            print(f"PASS {path.stem}::{name} ({time.perf_counter() - t0:.1f}s)")
        except Exception:
            failures.append(f"{path.stem}::{name}")
            print(f"FAIL {path.stem}::{name}\n{traceback.format_exc()}")

print(f"\n{passed} passed, {len(failures)} failed  (Blender {bpy.app.version_string})")
for f in failures:
    print("  FAILED", f)
sys.exit(1 if failures else 0)
