"""Numeric ChimeraX palettes, residue coloring, and deterministic fallbacks."""

from random import getstate
from types import SimpleNamespace

import numpy as np

from chimerax_style_in_pymol.colors import chain_color, colors, model_color, rgb
from chimerax_style_in_pymol.palettes import ELEMENTS


class Command:
    def get(self, name):
        return "black"

    def get_color_tuple(self, value):
        return (0, 0, 0) if value == "black" else (0.1, 0.2, 0.3)


def atom(**kwargs):
    return SimpleNamespace(
        segi="",
        chain=kwargs.get("chain", "A"),
        resi=kwargs.get("resi", "1"),
        resn=kwargs.get("resn", "ALA"),
        ss="",
        symbol=kwargs.get("symbol", "C"),
        b=kwargs.get("b", 0),
        color=0,
    )


def test_native_model_and_named_colors():
    expected = [
        [210, 180, 140],
        [135, 206, 235],
        [221, 160, 221],
        [144, 238, 144],
        [250, 128, 114],
        [211, 211, 211],
        [255, 20, 147],
        [255, 215, 0],
        [30, 144, 255],
        [128, 0, 128],
    ]
    for index, color in enumerate(expected):
        np.testing.assert_allclose(model_color(index, [0, 0, 0]), np.array(color) / 255)
    np.testing.assert_allclose(
        rgb("salmon", Command()), [250 / 255, 128 / 255, 114 / 255]
    )
    state = getstate()
    assert not np.array_equal(
        model_color(0, np.array(expected[0]) / 255), np.array(expected[0]) / 255
    )
    np.testing.assert_array_equal(
        model_color(12, [0, 0, 0]), model_color(12, [0, 0, 0])
    )
    assert getstate() == state


def test_chain_colors_follow_ids_without_global_rng_side_effects():
    model = SimpleNamespace(
        atom=[atom(chain="B"), atom(chain="A"), atom(chain="b"), atom(chain="AA")]
    )
    before = getstate()
    col = colors(model, "chain", Command())
    np.testing.assert_allclose(
        col[:3] * 255, [[240, 128, 128], [123, 104, 238], [240, 128, 128]]
    )
    np.testing.assert_array_equal(col[3], chain_color("aa"))
    assert getstate() == before


def test_element_palette_covers_metals_and_keeps_carbon_model_color():
    assert len(ELEMENTS) == 109
    model = SimpleNamespace(atom=[atom(symbol=s) for s in ["C", "O", "CU", "AU", "X"]])
    col = colors(model, "auto", Command(), atomic=True)
    np.testing.assert_allclose(
        col * 255,
        [
            [210, 180, 140],
            [255, 13, 13],
            [200, 128, 51],
            [255, 209, 35],
            [180, 180, 180],
        ],
    )
    np.testing.assert_allclose(
        colors(model, "element", Command())[0] * 255, [144, 144, 144]
    )


def test_rainbow_uses_residues_and_restarts_for_each_chain():
    model = SimpleNamespace(
        atom=[
            atom(chain=c, resi=r)
            for c, r in [
                ("A", "1"),
                ("A", "1"),
                ("A", "2"),
                ("A", "3"),
                ("B", "1"),
                ("B", "2"),
            ]
        ]
    )
    np.testing.assert_allclose(
        colors(model, "rainbow", Command()),
        [[0, 0, 1], [0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, 1], [1, 0, 0]],
    )


def test_bfactor_uses_blue_white_red_and_averages_cartoon_residues():
    model = SimpleNamespace(
        atom=[atom(resi=r, b=b) for r, b in [("1", 0), ("1", 2), ("2", 2), ("3", 3)]]
    )
    np.testing.assert_allclose(
        colors(model, "bfactor", Command()),
        [[0, 0, 1], [0, 0, 1], [1, 1, 1], [1, 0, 0]],
    )
    model.atom = [atom(), atom()]
    np.testing.assert_allclose(colors(model, "bfactor", Command()), 1)


def test_nucleotide_colors_preserve_non_nucleic_atoms():
    model = SimpleNamespace(
        atom=[atom(resn=r) for r in ["DA", "DT", "DG", "DC", "U", "PSU", "ALA"]]
    )
    col = colors(model, "nucleotide", Command())
    np.testing.assert_allclose(
        col * 255,
        [
            [255, 64, 64],
            [64, 64, 255],
            [64, 255, 64],
            [255, 255, 64],
            [64, 255, 255],
            [211, 211, 211],
            [210, 180, 140],
        ],
    )
