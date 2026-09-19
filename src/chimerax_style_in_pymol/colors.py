"""ChimeraX-like model, element, nucleotide, and attribute coloring."""

from functools import lru_cache
from random import Random

import numpy as np

from .palettes import CHAINS, ELEMENTS, NAMED

MODEL_COLORS = (
    "tan",
    "sky blue",
    "plum",
    "light green",
    "salmon",
    "light gray",
    "deep pink",
    "gold",
    "dodger blue",
    "purple",
)
BASE_COLORS = {
    "A": "#ff4040",
    "G": "#40ff40",
    "C": "#ffff40",
    "T": "#4040ff",
    "U": "#40ffff",
    "I": "#006400",
    "PSU": "#d3d3d3",
}


def rgb(value, cmd=None):
    if isinstance(value, (list, tuple, np.ndarray)):
        result = np.asarray(value, float)
    elif str(value).startswith("#") and len(str(value)) == 7:
        result = np.array([int(str(value)[i : i + 2], 16) / 255 for i in (1, 3, 5)])
    elif str(value).lower() in NAMED:
        result = np.asarray(NAMED[str(value).lower()], float) / 255
    elif cmd is not None:
        result = np.asarray(cmd.get_color_tuple(str(value)), float)
    else:
        known = {
            "white": "#ffffff",
            "black": "#000000",
            "red": "#ff0000",
            "blue": "#0000ff",
            "green": "#00ff00",
            "yellow": "#ffff00",
            "gray": "#808080",
        }
        if value not in known:
            raise ValueError(f"Use an RGB triple or #rrggbb color: {value}")
        return rgb(known[value])
    if (
        result.shape != (3,)
        or not np.isfinite(result).all()
        or np.any((result < 0) | (result > 1))
    ):
        raise ValueError("Colors must contain three finite values between 0 and 1")
    return result


@lru_cache(maxsize=512)
def chain_color(chain):
    key = str(chain).lower()
    color = CHAINS.get(key)
    if color is None:
        random = Random(key)
        color = tuple(random.randint(128, 255) for _ in range(3))
    return np.asarray(color, float) / 255


def model_color(index, background):
    if index < len(MODEL_COLORS):
        color = rgb(MODEL_COLORS[index])
        if not np.array_equal(
            np.rint(color * 255), np.rint(np.asarray(background) * 255)
        ):
            return color
    random = Random(index + 1)
    candidates = np.array([[random.random() for _ in range(3)] for _ in range(7)])
    avoid = np.array(
        [rgb(c) for c in MODEL_COLORS] + [[0, 0, 0], [0, 1, 0], [1, 1, 1], background]
    )
    distance = (abs(candidates[:, None] - avoid) * [1, 1, 0.5]).sum(axis=2).min(axis=1)
    return np.rint(candidates[distance.argmax()] * 255) / 255


def colors(model, mode, cmd, model_index=0, atomic=False):
    atoms = model.atom
    if not atoms:
        return np.empty((0, 3))
    background = cmd.get_color_tuple(cmd.get("bg_rgb"))
    base = model_color(model_index, background)
    result = np.tile(base, (len(atoms), 1))
    if mode == "keep":
        palette = {
            index: cmd.get_color_tuple(index) for index in {a.color for a in atoms}
        }
        return np.array([palette[a.color] for a in atoms])
    if mode == "chain":
        return np.array([chain_color(a.chain) for a in atoms])
    if mode == "secondary-structure":
        palette = {"H": rgb("#ff6961"), "S": rgb("#ffd966"), "": rgb("#a6c8ff")}
        return np.array([palette.get(a.ss, palette[""]) for a in atoms])
    if mode == "nucleotide":
        import gemmi

        palette = {key: rgb(value) for key, value in BASE_COLORS.items()}
        for i, atom in enumerate(atoms):
            info = gemmi.find_tabulated_residue(atom.resn)
            if info.is_nucleic_acid():
                result[i] = palette.get(
                    atom.resn,
                    palette.get(info.one_letter_code.upper(), [128 / 255] * 3),
                )
        return result
    if mode in ("rainbow", "bfactor"):
        # Native rainbow assigns one color per residue and restarts each chain.
        keys = [(a.segi, a.chain, a.resi) for a in atoms]
        if mode == "rainbow":
            chains = {}
            for key in keys:
                chains.setdefault(key[:2], {}).setdefault(key, None)
            fractions = {}
            for residues in chains.values():
                fractions.update(zip(residues, np.linspace(0, 1, len(residues))))
            values = np.array([fractions[key] for key in keys])
            palette = np.array([[0, 0, 1], [0, 1, 1], [0, 1, 0], [1, 1, 0], [1, 0, 0]])
        else:
            values = np.array([a.b for a in atoms], float)
            if not atomic:
                totals, counts = {}, {}
                for key, value in zip(keys, values):
                    totals[key] = totals.get(key, 0) + value
                    counts[key] = counts.get(key, 0) + 1
                values = np.array([totals[key] / counts[key] for key in keys])
            span = float(np.ptp(values))
            values = (
                (values - values.min()) / span if span else np.full(len(atoms), 0.5)
            )
            palette = np.array([[0, 0, 1], [1, 1, 1], [1, 0, 0]])
        stops = np.linspace(0, 1, len(palette))
        return (
            np.rint(
                np.array([np.interp(values, stops, palette[:, i]) for i in range(3)]).T
                * 255
            )
            / 255
        )
    if mode not in ("auto", "model", "element"):
        return np.tile(rgb(mode, cmd), (len(atoms), 1))
    if mode == "element" or (mode == "auto" and atomic):
        palette = {
            key: np.asarray(value, float) / 255 for key, value in ELEMENTS.items()
        }
        for i, atom in enumerate(atoms):
            key = atom.symbol.upper()
            if mode == "element" or key != "C":
                result[i] = palette.get(key, [180 / 255] * 3)
    return result
