"""Scientific geometry and input-contract regression checks."""

import subprocess
import sys

import numpy as np
import pytest
from samples import synthetic_protein

from chimerax_style_in_pymol.geometry import aniso, cartoon, ring_cycles
from chimerax_style_in_pymol.registry import parameters
from chimerax_style_in_pymol.volume import contour, read


def test_import_is_inert():
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import chimerax_style_in_pymol; assert 'pymol' not in sys.modules",
        ],
        check=True,
    )


def test_aniso_uses_tensor_axes_and_probability():
    from chempy import Atom
    from chempy.models import Indexed
    from scipy.stats import chi2

    model = Indexed()
    atom = Atom()
    atom.coord = [0, 0, 0]
    atom.u_aniso = [1, 4, 9, 0, 0, 0]
    model.atom = [atom]
    mesh = aniso(model, np.array([[1, 0, 0]]), "aniso", parameters(None), 32)
    np.testing.assert_allclose(np.max(abs(mesh.vertices), axis=0), [1, 2, 3], atol=1e-6)
    scaled = aniso(
        model, np.array([[1, 0, 0]]), "aniso", parameters({"probability": 0.5}), 32
    )
    np.testing.assert_allclose(
        np.max(abs(scaled.vertices), axis=0),
        np.array([1, 2, 3]) * np.sqrt(chi2.ppf(0.5, 3)),
        rtol=1e-6,
    )
    atom.u_aniso = [-1, 2, 3, 0, 0, 0]
    with pytest.raises(ValueError, match="positive semidefinite"):
        aniso(model, np.array([[1, 0, 0]]), "aniso", parameters(None), 16)


def test_gap_does_not_create_a_bridge():
    model = synthetic_protein()
    for atom in model.atom[36:]:
        atom.coord = (np.asarray(atom.coord) + [80, 0, 0]).tolist()
    mesh = cartoon(
        model,
        np.tile([0.2, 0.5, 1], (len(model.atom), 1)),
        "cartoon",
        parameters(None),
        12,
    )
    span = np.ptp(mesh.vertices[mesh.faces], axis=1)
    assert span.max() < 20


def test_skew_grid_contour_transform():
    x, y, z = np.mgrid[-1:1:15j, -1:1:15j, -1:1:15j]
    values = 1 - (x * x + y * y + z * z)
    affine = np.array(
        [[0.4, 0.1, 0, 7], [0, 0.5, 0.2, -4], [0, 0, 0.6, 3], [0, 0, 0, 1]]
    )
    grid = read({"values": values, "transform": affine}, 10000)
    mesh = contour(grid, 0.4, [0.4, 0.6, 1])
    assert np.isfinite(mesh.vertices).all()
    np.testing.assert_allclose(
        mesh.vertices.mean(axis=0), grid.world([7, 7, 7]), atol=0.1
    )
    np.testing.assert_allclose(np.linalg.norm(mesh.normals, axis=1), 1, atol=1e-6)
    with pytest.raises(ValueError, match="data range"):
        contour(grid, 9, [0, 0, 0])


@pytest.mark.parametrize(
    "params",
    [
        {"stick_radius": -1},
        {"ball_scale": float("nan")},
        {"probability": 1},
        {"stride": 1.3},
        {"typo": 1},
        {"show_hydrogens": "false"},
        {"helix_mode": "random"},
    ],
)
def test_bad_parameters(params):
    with pytest.raises(ValueError):
        parameters(params)


def test_rings_follow_bonds_only():
    from chempy import Atom, Bond
    from chempy.models import Indexed

    model = Indexed()
    for i in range(6):
        a = Atom()
        a.resi = "1"
        a.chain = "A"
        model.atom.append(a)
        b = Bond()
        b.index = [i, (i + 1) % 6]
        model.bond.append(b)
    assert len(ring_cycles(model)) == 1
    model.bond.pop()
    assert ring_cycles(model) == []


def test_mrc_origin_start_and_voxel_axes(tmp_path):
    import gemmi

    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = gemmi.FloatGrid(6, 8, 10)
    ccp4.grid.set_unit_cell(gemmi.UnitCell(12, 24, 40, 90, 90, 90))
    np.array(ccp4.grid, copy=False)[:] = np.arange(480).reshape(6, 8, 10)
    ccp4.update_ccp4_header()
    for k, value in zip((5, 6, 7), (2, 3, 4)):
        ccp4.set_header_i32(k, value)
    for k, value in zip((50, 51, 52), (10, 20, 30)):
        ccp4.set_header_float(k, value)
    path = tmp_path / "offset.mrc"
    ccp4.write_ccp4_map(str(path))
    result = read({"file": path}, 1000)
    np.testing.assert_allclose(result.world([0, 0, 0]), [14, 29, 46])
    np.testing.assert_allclose(result.world([1, 1, 1]), [16, 32, 50])
    np.testing.assert_array_equal(result.values, np.array(ccp4.grid))


def test_zero_dash_count_creates_a_solid_pseudobond():
    from chimerax_style_in_pymol.primitives import dashed

    p = parameters({"dashes": 0})
    mesh = dashed([0, 0, 0], [2, 0, 0], 0.1, [1, 0, 0], count=p["dashes"])
    assert len(mesh.faces) > 0


def test_batched_cgo_preserves_opacity_groups_and_triangle_attributes():
    from pymol import cgo as cg

    from chimerax_style_in_pymol.geometry import cgo
    from chimerax_style_in_pymol.mesh import Mesh

    vertices = np.arange(27, dtype=np.float32).reshape(9, 3) / 10
    normals = np.tile([0, 0, 1], (9, 1))
    colors = np.tile([0.2, 0.4, 0.8], (9, 1))
    mesh = Mesh(
        vertices,
        normals,
        colors,
        np.arange(9).reshape(3, 3),
        np.zeros(9),
        alphas=np.repeat([0.75, 0, 0.25], 3),
    )
    values = cgo(mesh)
    offset = 2
    for face, opacity in ((2, 0.25), (0, 0.75)):
        assert values[offset : offset + 4] == [
            cg.ALPHA,
            opacity,
            cg.BEGIN,
            cg.TRIANGLES,
        ]
        rows = np.array(values[offset + 4 : offset + 40]).reshape(3, 12)
        np.testing.assert_array_equal(rows[:, 0], cg.COLOR)
        np.testing.assert_array_equal(rows[:, 4], cg.NORMAL)
        np.testing.assert_array_equal(rows[:, 8], cg.VERTEX)
        ids = mesh.faces[face]
        np.testing.assert_array_equal(rows[:, 1:4], mesh.colors[ids])
        np.testing.assert_array_equal(rows[:, 5:8], mesh.normals[ids])
        np.testing.assert_array_equal(rows[:, 9:12], mesh.vertices[ids])
        assert values[offset + 40] == cg.END
        offset += 41
    assert offset == len(values)
