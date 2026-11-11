"""CONFIG: every user-settable parameter of a run.

See docs/agents/data_definitions.md §7. Imports nothing from the package.
Units: specimen SI metres for lengths of independent structures; facet-tied
lengths (rim_width) are multiples of SPACING; dome height is a fraction of the
facet radius (inflation); see §1 of the DDs.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields
from pathlib import Path

# Allowed values of the enum-like fields.
DENSITY_MODES = ("SPACING", "COUNT")
BRISTLE_PATTERNS = ("ALTERNATING", "RANDOM", "NONE")
NIPPLE_MODES = ("OFF", "SHADER", "GEOMETRY")


# CDD. CONFIG
# config = Config(...)
# interp. every user-settable parameter of a run. Serialised to config.json next
#         to every result and stored as JSON on the generated lattice object.
#         The UI PropertyGroup is a view of it (converted only in config_props.py).
# Invariant: validate() returns without raising.
@dataclass
class Config:
    # density
    density_mode: str = "COUNT"            # DENSITY_MODE
    spacing: float = 16e-6                 # SPACING (m); dossier §3, estimated from SEM
    # ~700 (male) to ~750 (female) ommatidia per eye: "An anatomical atlas of
    # Drosophila melanogaster - the wild-type" (compound-eye overview).
    target_count: int = 750                # TARGET_COUNT
    max_count: int = 50_000                # MAX_COUNT
    model_scale: float = 1.0               # MODEL_SCALE (Blender units per specimen metre)
    # lattice
    seed: int = 0                          # SEED (only the seams are random)
    # Where the perfect hexagonal patch starts (specimen metres, world axes).
    # None = automatic: the point farthest from the rim on an open sheet, the
    # surface point facing -Y (Blender's front view) on a closed mesh.
    lattice_origin: tuple | None = None    # LATTICE_ORIGIN
    lattice_rotation_deg: float = 0.0      # LATTICE_ROTATION (row direction about the normal)
    relax_iterations: int = 200            # RELAX_ITERATIONS (seam fill only)
    irregularity: float = 0.0              # IRREGULARITY (x spacing)
    # facet shape
    facet_shape: float = 0.0               # FACET_SHAPE (0 dome, 1 flat facet)
    inflation: float = 0.4                 # INFLATION (dome height / facet radius; 1 = semi-dome); guess
    pokeyness: float = 0.0                 # POKEYNESS (0 round cap, 1 pointed cone)
    crease_sharpness: float = 0.5          # CREASE_SHARPNESS
    rim_width: float = 1.5                 # RIM_WIDTH (x spacing)
    subdivision_level: int = -1            # facet mesh resolution; -1 = auto (2)
    # bristles: one per ommatidium in Drosophila (same atlas) -> ALTERNATING
    bristle_pattern: str = "ALTERNATING"   # BRISTLE_PATTERN
    bristle_frequency: float = 0.5         # BRISTLE_FREQUENCY (RANDOM only)
    bristle_length: float = 10e-6          # guess from SEM, dossier §8
    bristle_radius: float = 0.4e-6         # guess
    bristle_tilt_deg: float = 20.0         # guess (SEM: bristles fan out nearly along the normal)
    comb_direction: tuple = (0.0, 0.0, -1.0)
    socket_radius: float = 1.0e-6          # guess
    socket_depth: float = 0.3e-6           # guess
    # corneal nipples
    nipple_mode: str = "SHADER"            # NIPPLE_MODE
    nipple_spacing: float = 0.2e-6         # needs citation, dossier §8
    nipple_height: float = 0.05e-6         # needs citation

    # FD. validate
    # Signature: CONFIG -> None (raises ValueError)
    def validate(self) -> None:
        """Check every invariant of CONFIG and report all violations at once.

        Algorithm:
            for each field, evaluate its invariant; collect a message if broken
            if any messages: raise ValueError listing them
        """
        problems = [msg for ok, msg in self._checks() if not ok]
        if problems:
            raise ValueError("Invalid CONFIG:\n  - " + "\n  - ".join(problems))

    def _checks(self) -> list[tuple[bool, str]]:
        c = self
        return [
            (c.density_mode in DENSITY_MODES, f"density_mode must be one of {DENSITY_MODES}"),
            (_is_real(c.spacing) and c.spacing > 0, "spacing must be > 0"),
            (_is_int(c.target_count) and c.target_count >= 1, "target_count must be an int >= 1"),
            (_is_int(c.max_count) and c.max_count >= 1, "max_count must be an int >= 1"),
            (_is_real(c.model_scale) and c.model_scale > 0, "model_scale must be > 0"),
            (_is_int(c.seed) and c.seed >= 0, "seed must be an int >= 0"),
            (c.lattice_origin is None or _is_vec3(c.lattice_origin),
             "lattice_origin must be None or 3 finite numbers"),
            (_is_real(c.lattice_rotation_deg), "lattice_rotation_deg must be a finite number"),
            (_is_int(c.relax_iterations) and c.relax_iterations >= 0, "relax_iterations must be an int >= 0"),
            (_in_range(c.irregularity, 0.0, 0.3), "irregularity must be in [0, 0.3]"),
            (_in_range(c.facet_shape, 0.0, 1.0), "facet_shape must be in [0, 1]"),
            (_in_range(c.inflation, 0.0, 1.0), "inflation must be in [0, 1]"),
            (_in_range(c.pokeyness, 0.0, 1.0), "pokeyness must be in [0, 1]"),
            (_in_range(c.crease_sharpness, 0.0, 1.0), "crease_sharpness must be in [0, 1]"),
            (_is_real(c.rim_width) and c.rim_width >= 0, "rim_width must be >= 0"),
            (_is_int(c.subdivision_level) and -1 <= c.subdivision_level <= 6,
             "subdivision_level must be an int in [-1, 6]"),
            (c.bristle_pattern in BRISTLE_PATTERNS, f"bristle_pattern must be one of {BRISTLE_PATTERNS}"),
            (_in_range(c.bristle_frequency, 0.0, 1.0), "bristle_frequency must be in [0, 1]"),
            (_is_real(c.bristle_length) and c.bristle_length > 0, "bristle_length must be > 0"),
            (_is_real(c.bristle_radius) and c.bristle_radius > 0, "bristle_radius must be > 0"),
            (_in_range(c.bristle_tilt_deg, 0.0, 90.0), "bristle_tilt_deg must be in [0, 90]"),
            (_is_nonzero_vec3(c.comb_direction), "comb_direction must be 3 finite numbers, not all zero"),
            (_is_real(c.socket_radius) and c.socket_radius >= 0, "socket_radius must be >= 0"),
            (_is_real(c.socket_depth) and c.socket_depth >= 0, "socket_depth must be >= 0"),
            (c.nipple_mode in NIPPLE_MODES, f"nipple_mode must be one of {NIPPLE_MODES}"),
            (_is_real(c.nipple_spacing) and c.nipple_spacing > 0, "nipple_spacing must be > 0"),
            (_is_real(c.nipple_height) and c.nipple_height >= 0, "nipple_height must be >= 0"),
        ]

    # FD. to_json
    # Signature: CONFIG -> str
    def to_json(self) -> str:
        """Serialise CONFIG as indented JSON (tuples become lists)."""
        return json.dumps(asdict(self), indent=2)

    # FD. from_json
    # Signature: str -> CONFIG
    @classmethod
    def from_json(cls, text: str) -> "Config":
        """Parse JSON into a CONFIG and validate it.

        Algorithm:
            parse the JSON object
            reject unknown keys (a typo must not silently fall back to a default)
            convert vector fields back to tuples
            build Config, validate, return
        """
        data = json.loads(text)
        known = {f.name for f in fields(cls)}
        unknown = sorted(set(data) - known)
        if unknown:
            raise ValueError(f"Unknown CONFIG keys: {unknown}")
        for key in ("comb_direction", "lattice_origin"):
            if data.get(key) is not None:
                data[key] = tuple(data[key])
        config = cls(**data)
        config.validate()
        return config

    # FD. save
    # Signature: CONFIG, path -> None
    def save(self, path: str | Path) -> None:
        """Write CONFIG as config JSON to path."""
        Path(path).write_text(self.to_json(), encoding="utf-8")

    # FD. load
    # Signature: path -> CONFIG
    @classmethod
    def load(cls, path: str | Path) -> "Config":
        """Read and validate a CONFIG from a JSON file."""
        return cls.from_json(Path(path).read_text(encoding="utf-8"))


def _is_real(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and x == x and abs(x) != float("inf")


def _is_int(x) -> bool:
    return isinstance(x, int) and not isinstance(x, bool)


def _in_range(x, lo: float, hi: float) -> bool:
    return _is_real(x) and lo <= x <= hi


def _is_vec3(v) -> bool:
    return isinstance(v, (tuple, list)) and len(v) == 3 and all(_is_real(x) for x in v)


def _is_nonzero_vec3(v) -> bool:
    return _is_vec3(v) and any(x != 0 for x in v)


config0 = Config()  # stub: Drosophila defaults
config1 = Config(density_mode="SPACING", spacing=25e-6, facet_shape=1.0,
                 inflation=0.1, bristle_pattern="NONE")  # wasp-like
