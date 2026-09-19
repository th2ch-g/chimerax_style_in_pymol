"""Exercise representations and lifecycle using a real isolated PyMOL instance."""

from unittest.mock import patch

import numpy as np
import pymol2
import pytest
from samples import load

from chimerax_style_in_pymol import __init_plugin__, chimerax_style
from chimerax_style_in_pymol.controller import storage
from chimerax_style_in_pymol.registry import STYLES


@pytest.fixture(scope="session")
def cmd():
    import pymol

    pm = pymol2.SingletonPyMOL()
    pm.start()
    yield pymol.cmd
    pm.stop()


@pytest.mark.parametrize("style", STYLES)
def test_every_style(cmd, style):
    data, params = load(cmd, style)
    result = chimerax_style(
        style, data=data, params=params, quality="low", quiet=1, _self=cmd
    )
    assert result in cmd.get_names("all")
    assert storage(cmd)["views"]["chimerax"]["objects"]
    chimerax_style("reset", _self=cmd)
    assert not [o for o in cmd.get_names("all") if o.startswith("_cxs_")]


def snapshot(cmd, selection="sample"):
    rows = []
    cmd.iterate(
        selection,
        "rows.append((index,reps,color,ss,ID,name,resi))",
        space={"rows": rows},
    )
    return rows, cmd.get_coords(selection).copy()


def test_preservation_overlap_and_reset(cmd):
    load(cmd, "cartoon")
    cmd.show("lines", "sample")
    cmd.show("spheres", "sample and resi 1")
    cmd.color("red", "sample and resi 3")
    before, xyz = snapshot(cmd)
    ambient = cmd.get("ambient")
    chimerax_style("cartoon", selection="sample and resi 1-12", name="one", _self=cmd)
    chimerax_style("ball", selection="sample and resi 6-18", name="two", _self=cmd)
    chimerax_style("surface", selection="sample and resi 1-12", name="one", _self=cmd)
    chimerax_style("reset", name="one", _self=cmd)
    hidden = []
    cmd.iterate("sample and resi 6-18", "hidden.append(reps)", space={"hidden": hidden})
    assert not any(hidden)
    chimerax_style("reset", name="two", _self=cmd)
    after, coords = snapshot(cmd)
    assert after == before
    np.testing.assert_array_equal(xyz, coords)
    assert cmd.get("ambient") == ambient


def test_failure_keeps_previous_view(cmd):
    load(cmd, "cartoon")
    chimerax_style("cartoon", _self=cmd)
    objects = cmd.get_names("all")
    with pytest.raises(ValueError, match="ANISOU"):
        chimerax_style("aniso", _self=cmd)
    assert cmd.get_names("all") == objects
    with pytest.raises(ValueError):
        chimerax_style("ball", params={"stick_radius": -1}, _self=cmd)
    assert cmd.get_names("all") == objects


def test_states_refresh_session_and_pick(cmd):
    load(cmd, "cartoon")
    cmd.create("sample", "sample", 1, 2)
    cmd.translate([3, 1, 0], "sample", state=2)
    before, _ = snapshot(cmd)
    chimerax_style("cartoon", _self=cmd)
    entry = storage(cmd)["views"]["chimerax"]
    assert all(cmd.count_states(o) == 2 for o in entry["objects"])
    saved = cmd.get_session()
    cmd.reinitialize()
    cmd.set_session(saved)
    chimerax_style("refresh", _self=cmd)
    chimerax_style("ball", _self=cmd)
    obj = next(iter(storage(cmd)["views"]["chimerax"]["copies"]))
    assert chimerax_style("select", selection=f"{obj} and index 1", _self=cmd) == 1
    assert cmd.index("sele") == [("sample", 1)]
    chimerax_style("reset", _self=cmd)
    assert snapshot(cmd)[0] == before


def test_registration_and_parser(cmd):
    load(cmd, "cartoon")
    __init_plugin__(_self=cmd)
    cmd.do("chimerax_style ball, selection=sample, name=parser, lighting=keep")
    assert "parser" in storage(cmd)["views"]


@pytest.mark.parametrize("color", ["auto", "element", "model", "keep", "red"])
def test_atomic_palettes_do_not_copy_coordinates_or_topology(cmd, color):
    load(cmd, "cartoon")
    before = snapshot(cmd)
    with patch.object(cmd, "get_model", wraps=cmd.get_model) as read_model:
        chimerax_style("ball", color=color, quiet=1, _self=cmd)
    assert read_model.call_count == 0
    chimerax_style("reset", _self=cmd)
    assert snapshot(cmd)[0] == before[0]
    np.testing.assert_array_equal(snapshot(cmd)[1], before[1])


def test_cartoon_keep_colors_reads_native_color_indices(cmd):
    load(cmd, "cartoon")
    cmd.color("red", "sample and resi 1-4")
    before = snapshot(cmd)
    chimerax_style("default", color="keep", quiet=1, _self=cmd)
    chimerax_style("reset", _self=cmd)
    assert snapshot(cmd)[0] == before[0]
    np.testing.assert_array_equal(snapshot(cmd)[1], before[1])


def test_existing_group_and_unrelated_object(cmd):
    load(cmd, "cartoon")
    cmd.create("other", "sample")
    cmd.color("yellow", "other")
    before = snapshot(cmd, "other")
    cmd.group("chimerax_chimerax", "other")
    with pytest.raises(ValueError, match="not owned"):
        chimerax_style("cartoon", selection="sample", _self=cmd)
    assert snapshot(cmd, "other")[0] == before[0]


def test_ray_export(cmd, tmp_path):
    from PIL import Image

    load(cmd, "cartoon")
    chimerax_style("cartoon", _self=cmd)
    cmd.orient("chimerax_chimerax")
    path = tmp_path / "render.png"
    chimerax_style("ray", filename=path, width=320, height=240, _self=cmd)
    pixels = np.array(Image.open(path).convert("RGB"))
    assert pixels.shape == (240, 320, 3)
    assert np.count_nonzero(pixels.min(axis=2) < 220) > 500


def test_map_object_transform(cmd):
    from chimerax_style_in_pymol.mapio import from_pymol

    load(cmd, "cartoon")
    cmd.map_new("density", "gaussian", 1.2, "sample")
    grid = from_pymol(cmd, "density", 1, 8000000)
    corners = grid.world([[0, 0, 0], np.array(grid.values.shape) - 1])
    np.testing.assert_allclose(corners, cmd.get_extent("density"), atol=1e-4)
    chimerax_style("volume-surface", selection="density", _self=cmd)


def test_user_lighting_change_is_preserved(cmd):
    load(cmd, "cartoon")
    chimerax_style("full", _self=cmd)
    cmd.set("ambient", 0.123)
    chimerax_style("reset", _self=cmd)
    assert float(cmd.get("ambient")) == pytest.approx(0.123, abs=1e-4)


def test_projection_and_transformed_sources(cmd):
    load(cmd, "cartoon")
    cmd.set("orthoscopic", 0)
    transform = np.eye(4)
    transform[:3, 3] = [12, -4, 7]
    cmd.transform_object("sample", transform.ravel().tolist())
    before = cmd.get_coords("sample").copy()
    camera = cmd.get_view()
    chimerax_style("cartoon", _self=cmd)
    obj = next(iter(storage(cmd)["views"]["chimerax"]["copies"]))
    np.testing.assert_allclose(cmd.get_coords(obj), before, atol=1e-5)
    np.testing.assert_allclose(cmd.get_view()[:17], camera[:17], atol=1e-4)
    assert cmd.get_setting_int("orthoscopic") == 1
    chimerax_style("reset", _self=cmd)
    assert cmd.get_setting_int("orthoscopic") == 0


@pytest.mark.parametrize("payload", [None, -1, b"incomplete PNG"])
def test_failed_export_preserves_existing_file(cmd, tmp_path, monkeypatch, payload):
    path = tmp_path / "existing.png"
    path.write_bytes(b"previous image")
    monkeypatch.setattr(cmd, "png", lambda *args, **kwargs: payload)
    with pytest.raises(RuntimeError, match="did not write a PNG"):
        chimerax_style("png", filename=path, width=100, height=100, _self=cmd)
    assert path.read_bytes() == b"previous image"
    assert list(tmp_path.iterdir()) == [path]


def test_publication_background_uses_python_scalars_and_restores(cmd):
    from chimerax_style_in_pymol import chimerax_style

    cmd.pseudoatom("background_probe", pos=[0, 0, 0])
    cmd.bg_color("black")
    before = cmd.get_setting_tuple("bg_rgb")
    chimerax_style("publication", selection="background_probe", _self=cmd)
    np.testing.assert_allclose(cmd.get_color_tuple(cmd.get("bg_rgb")), [1, 1, 1])
    chimerax_style("reset", _self=cmd)
    assert cmd.get_setting_tuple("bg_rgb") == before


@pytest.mark.parametrize(
    "quality, spacing", [("high", 0.5), ("medium", 0.75), ("low", 1.0)]
)
def test_surface_sampling_is_scoped_and_restored(cmd, quality, spacing):
    load(cmd, "surface-mesh")
    before = snapshot(cmd)
    original = cmd.get_setting_float("surface_best")
    chimerax_style("surface-mesh", quality=quality, quiet=1, _self=cmd)
    obj = next(iter(storage(cmd)["views"]["chimerax"]["copies"]))
    native_quality = cmd.get_setting_int("surface_quality", obj)
    actual = cmd.get_setting_float(
        "surface_normal" if native_quality == 0 else "surface_best", obj
    )
    if native_quality == 2:
        actual /= 2
    assert actual == pytest.approx(spacing)
    assert cmd.get_setting_float("surface_best") == original
    chimerax_style("reset", _self=cmd)
    after = snapshot(cmd)
    assert before[0] == after[0]
    np.testing.assert_array_equal(before[1], after[1])
