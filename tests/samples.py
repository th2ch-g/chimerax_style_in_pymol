"""Deterministic local inputs shared by tests and the documentation gallery."""

import numpy as np

from chimerax_style_in_pymol.registry import GROUPS


def synthetic_protein():
    from chempy import Atom, Bond
    from chempy.models import Indexed

    model = Indexed()
    for i in range(18):
        t = i * np.deg2rad(100)
        center = np.array([2.3 * np.cos(t), 2.3 * np.sin(t), i * 1.5])
        if i >= 10:
            center = np.array([2.3 + (i - 10) * 3.2, 0, 15])
        for name, offset, symbol in (
            ("CA", [0, 0, 0], "C"),
            ("N", [-0.5, 0, -0.6], "N"),
            ("C", [0.6, 0, 0.6], "C"),
            ("O", [0.7, 1, 0.7], "O"),
        ):
            atom = Atom()
            atom.name = name
            atom.symbol = symbol
            atom.resn = "ALA"
            atom.resi = str(i + 1)
            atom.chain = "A"
            atom.coord = list(center + offset)
            atom.ss = "H" if i < 9 else "S" if i > 11 else "L"
            atom.b = 10 + i
            atom.vdw = 1.7
            model.atom.append(atom)
        for a, b in ((i * 4, i * 4 + 1), (i * 4, i * 4 + 2), (i * 4 + 2, i * 4 + 3)):
            bond = Bond()
            bond.index = [a, b]
            model.bond.append(bond)
        if i:
            bond = Bond()
            bond.index = [i * 4 - 2, i * 4 + 1]
            model.bond.append(bond)
    return model


def sugars():
    from chempy import Atom, Bond
    from chempy.models import Indexed

    model = Indexed()
    for k, resn in enumerate(("NAG", "MAN", "GAL", "FUC", "SIA", "XYS", "FRU", "KDO")):
        center = np.array([(k % 4) * 5, (k // 4) * 5, 0])
        for j in range(6):
            t = j * np.pi / 3
            a = Atom()
            a.name = f"C{j + 1}"
            a.symbol = "C"
            a.resn = resn
            a.resi = str(k + 1)
            a.chain = "G"
            a.hetatm = 1
            a.coord = list(center + [np.cos(t), np.sin(t), 0.1 * (-1) ** j])
            a.vdw = 1.7
            model.atom.append(a)
            b = Bond()
            b.index = [k * 6 + j, k * 6 + (j + 1) % 6]
            model.bond.append(b)
        if k:
            b = Bond()
            b.index = [k * 6 - 1, k * 6]
            model.bond.append(b)
    return model


def grid(segment=False):
    x, y, z = np.mgrid[-1:1:24j, -1:1:22j, -1:1:20j]
    values = 2.8 * np.exp(
        -((x + 0.3) ** 2 / 0.25 + (y - 0.1) ** 2 / 0.4 + z**2 / 0.45)
    ) + 2 * np.exp(
        -((x - 0.45) ** 2 / 0.15 + (y + 0.2) ** 2 / 0.2 + (z - 0.25) ** 2 / 0.2)
    )
    if segment:
        values = np.where(values > 1.3, 1, np.where(values > 0.7, 2, 0))
    return {"values": values, "origin": [-6, -5, -4], "spacing": [0.5, 0.5, 0.5]}


def load(cmd, style, structure=None):
    from chimerax_style_in_pymol import chimerax_style

    chimerax_style("reset", name="all", _self=cmd)
    cmd.delete("all")
    data = {}
    params = {}
    if style in GROUPS["Maps"]:
        data = grid(style == "volume-segment")
        if style == "volume-segment":
            params["segments"] = [1, 2]
        if style == "volume-slab":
            params["slab_normal"] = [0.2, 0.3, 1]
    elif style in GROUPS["Shapes"] or style in ("markers", "pseudobonds"):
        data = {
            "points": [[-3, -1, 0], [-1, 2, 0], [1, 1, 1], [3, -1, 0]],
            "links": [[0, 1], [1, 2], [2, 3]],
            "pairs": [[0, 1], [1, 2], [2, 3]],
            "faces": [[0, 1, 2], [0, 2, 3]],
        }
        if style == "shape-triangle":
            data["points"] = data["points"][:3]
            data["faces"] = [[0, 1, 2]]
        if style == "markers":
            params["radius"] = 0.5
        if style in ("shape-sphere", "shape-icosahedron", "shape-dodecahedron"):
            params["radius"] = 2.5
    elif style.startswith("nucleotides-"):
        cmd.fnab("ATGCAT", name="sample")
        if style == "nucleotides-ladder":
            data = {"pairs": [[i, 11 - i] for i in range(6)]}
    elif style == "snfg":
        cmd.load_model(sugars(), "sample")
    elif (
        style in GROUPS["Atoms"]
        or style.startswith("aniso")
        or style in ("distance", "angle", "torsion", "label", "hbonds", "contacts")
    ):
        cmd.fragment("trp", "sample")
        cmd.remove("sample and hydro")
        if style.startswith("aniso"):
            model = cmd.get_model("sample")
            for a in model.atom:
                a.u_aniso = [0.20, 0.30, 0.5, 0.04, 0.02, 0.03]
            cmd.delete("sample")
            cmd.load_model(model, "sample")
        if style in ("distance", "angle", "torsion"):
            data = {
                "indices": list(range({"distance": 2, "angle": 3, "torsion": 4}[style]))
            }
        if style == "label":
            params["label_level"] = "atom"
        if style == "hbonds":
            cmd.create("partner", "sample")
            cmd.translate([0, 3, 0], "partner")
    elif structure:
        cmd.load(str(structure), "sample")
        cmd.remove("sample and not polymer.protein")
        cmd.dss("sample")
    else:
        cmd.load_model(synthetic_protein(), "sample")
    if style == "unitcell":
        cmd.set_symmetry("sample", 35, 30, 25, 90, 90, 90, "P 1")
    cmd.bg_color("white")
    cmd.set("max_threads", 4)
    return data, params


def caption(style):
    if style in GROUPS["Maps"]:
        return "Synthetic scalar field (integer labels for segmentation)"
    if style in GROUPS["Shapes"] or style in ("markers", "pseudobonds"):
        return "Synthetic geometric inputs"
    if style.startswith("nucleotides-"):
        return "DNA built with PyMOL fnab; explicit complementary pairs for ladder"
    if style == "snfg":
        return "Synthetic glycan residues and bonds; not an experimental glycan"
    if style.startswith("aniso"):
        return "Tryptophan fragment with synthetic anisotropic displacement tensors"
    if style in GROUPS["Atoms"] or style in (
        "distance",
        "angle",
        "torsion",
        "label",
        "hbonds",
        "contacts",
    ):
        return "PyMOL tryptophan fragment; translated copy for hydrogen-bond example"
    return "Crambin (PDB 1CRN); demonstration cell for unitcell"
