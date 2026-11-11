"""CONFIG: defaults are valid, every invariant rejects a mutation, JSON round-trips."""

import json

import pytest

from ommatidia.config import Config


def test_defaults_validate():
    Config().validate()


def test_json_round_trip(tmp_path):
    c = Config(seed=7, comb_direction=(1.0, 0.0, 0.0), density_mode="SPACING")
    path = tmp_path / "config.json"
    c.save(path)
    assert Config.load(path) == c


def test_unknown_key_rejected():
    data = json.loads(Config().to_json())
    data["spacnig"] = 1.0
    with pytest.raises(ValueError, match="spacnig"):
        Config.from_json(json.dumps(data))


@pytest.mark.parametrize(
    "field, bad",
    [
        ("density_mode", "AREA"),
        ("spacing", 0.0),
        ("spacing", float("nan")),
        ("target_count", 0),
        ("target_count", 7.5),
        ("max_count", 0),
        ("model_scale", -1.0),
        ("seed", -1),
        ("seed", True),
        ("relax_iterations", -1),
        ("irregularity", 0.31),
        ("facet_shape", 1.01),
        ("pokeyness", 1.01),
        ("inflation", -0.1),
        ("subdivision_level", 7),
        ("lattice_origin", (0.0, 1.0)),
        ("crease_sharpness", -0.1),
        ("rim_width", -1.0),
        ("subdivision_level", -2),
        ("bristle_pattern", "ALL"),
        ("bristle_frequency", 2.0),
        ("bristle_length", 0.0),
        ("bristle_radius", 0.0),
        ("bristle_tilt_deg", 91.0),
        ("comb_direction", (0.0, 0.0, 0.0)),
        ("comb_direction", (1.0, 0.0)),
        ("socket_radius", -1e-6),
        ("socket_depth", -1e-6),
        ("nipple_mode", "ON"),
        ("nipple_spacing", 0.0),
        ("nipple_height", -1.0),
    ],
)
def test_each_invariant_rejects_mutation(field, bad):
    c = Config(**{field: bad})
    with pytest.raises(ValueError, match=field):
        c.validate()


def test_all_problems_reported_at_once():
    with pytest.raises(ValueError) as err:
        Config(spacing=-1.0, seed=-1).validate()
    assert "spacing" in str(err.value) and "seed" in str(err.value)
