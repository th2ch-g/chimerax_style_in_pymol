"""ChimeraX-like model, element, nucleotide, and attribute coloring."""

import colorsys

import numpy as np

MODEL_COLORS = ("#a6c8ff", "#ffb3b3", "#b2e5b2", "#ffff99", "#d9b3ff", "#99ffff")
ELEMENTS = {
    "H": "#ffffff",
    "C": "#909090",
    "N": "#3050f8",
    "O": "#ff0d0d",
    "S": "#ffff30",
    "P": "#ff8000",
    "F": "#90e050",
    "CL": "#1ff01f",
    "BR": "#a62929",
    "I": "#940094",
    "FE": "#e06633",
    "ZN": "#7d80b0",
    "MG": "#8aff00",
    "CA": "#3dff00",
    "NA": "#ab5cf2",
    "K": "#8f40d4",
}
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


def colors(model, mode, cmd, model_index=0, atomic=False):
    atoms = model.atom
    chains = list(dict.fromkeys((a.segi, a.chain) for a in atoms))
    base = rgb(MODEL_COLORS[model_index % len(MODEL_COLORS)])
    if mode not in (
        "auto",
        "model",
        "chain",
        "element",
        "keep",
        "nucleotide",
        "secondary-structure",
        "rainbow",
        "bfactor",
    ):
        return np.tile(rgb(mode, cmd), (len(atoms), 1))
    values = np.array([a.b for a in atoms])
    span = max(float(np.ptp(values)), 1e-8) if len(values) else 1
    result = []
    for i, a in enumerate(atoms):
        color = base
        if mode == "keep":
            color = cmd.get_color_tuple(a.color)
        elif mode == "chain":
            ci = chains.index((a.segi, a.chain))
            color = rgb(MODEL_COLORS[ci % len(MODEL_COLORS)])
        elif mode == "secondary-structure":
            color = rgb({"H": "#ff6961", "S": "#ffd966"}.get(a.ss, "#a6c8ff"))
        elif mode == "nucleotide":
            color = rgb(BASE_COLORS.get(a.resn.lstrip("D"), "#a6c8ff"))
        elif mode in ("rainbow", "bfactor"):
            f = (
                i / max(len(atoms) - 1, 1)
                if mode == "rainbow"
                else (a.b - values.min()) / span
            )
            color = colorsys.hsv_to_rgb((1 - f) * 0.667, 0.7, 0.95)
        if mode == "element" or (mode == "auto" and atomic and a.symbol.upper() != "C"):
            color = rgb(ELEMENTS.get(a.symbol.upper(), "#ff1493"))
        result.append(color)
    return np.asarray(result)
