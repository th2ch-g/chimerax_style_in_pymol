"""Build native PyMOL objects and CGO geometry without modifying source models."""

import numpy as np

from . import geometry, shapes, volume
from .colors import colors, rgb
from .registry import GROUPS


def recolor(cmd, obj, mode, model_index, atomic=False):
    model = cmd.get_model(obj)
    col = colors(model, mode, cmd, model_index, atomic)
    table = {}
    unique = {}
    for atom, color in zip(model.atom, col):
        key = tuple(np.round(color, 5))
        if key not in unique:
            name = "cxs_rgb_" + "".join(f"{round(v * 255):02x}" for v in color)
            cmd.set_color(name, list(color), quiet=1)
            unique[key] = cmd.get_color_index(name)
        table[atom.index] = unique[key]
    cmd.alter(obj, "color = table[index]", space={"table": table}, quiet=1)
    cmd.recolor(obj)


def atomic(cmd, obj, style, p):
    cmd.set("valence", 0, obj)
    cmd.set("stick_radius", p["stick_radius"], obj)
    cmd.set("stick_quality", 20, obj)
    cmd.set("sphere_quality", 3, obj)
    cmd.set("stick_ball", 1, obj)
    cmd.set("stick_ball_ratio", 1, obj)
    if style == "sphere":
        cmd.set("sphere_scale", 1, obj)
        cmd.show("spheres", obj)
    elif style == "ball":
        cmd.set("sphere_scale", p["ball_scale"], obj)
        cmd.show("sticks", obj)
        cmd.show("spheres", obj)
    else:
        cmd.show("sticks", obj)
    if not p["show_hydrogens"]:
        cmd.hide("everything", f"({obj}) and hydro")


def add_mesh(cmd, mesh, name, state, p, transparency, mode="solid"):
    if len(mesh.faces) > p["max_triangles"]:
        raise ValueError("Geometry exceeds max_triangles; reduce selection or quality")
    if not len(mesh.faces):
        raise ValueError("The requested representation produced no geometry")
    mesh.opacity *= 1 - transparency
    cmd.load_cgo(
        geometry.cgo(mesh, mode, stride=p["stride"]), name, state=state, zoom=0
    )
    cmd.set("cgo_lighting", 1, name)
    cmd.set("two_sided_lighting", 1, name)


def molecular(
    cmd, obj, source, style, color, quality, p, data, transparency, model_index
):
    detail = {"low": 12, "medium": 20, "high": 32}[quality]
    extras = []
    custom = None
    cmd.hide("everything", obj)
    atom_style = (
        style in GROUPS["Atoms"]
        or style in ("nucleotides-atoms", "nucleotides-fill")
        or style.startswith("aniso")
    )
    recolor(
        cmd,
        obj,
        color,
        model_index,
        atomic=atom_style
        or style in ("default", "ribbons-slabs", "cylinders-stubs", "licorice-ovals"),
    )
    for setting in (
        "transparency",
        "cartoon_transparency",
        "stick_transparency",
        "sphere_transparency",
    ):
        cmd.set(setting, transparency, obj)
    if style in ("stick", "ball", "sphere"):
        atomic(cmd, obj, style, p)
    elif style.startswith("ring-fill"):
        atomic(cmd, obj, "stick", p)
        custom = "rings"
    elif style in GROUPS["Cartoons"]:
        custom = "cartoon"
    elif style.startswith("surface") or style == "ghostly-white":
        # PyMOL quality 2 halves surface_best; match the reference 0.5 A
        # sampling scale instead of silently requesting a 0.125 A mesh.
        spacing = p["grid_spacing"] * {"high": 1, "medium": 1.5, "low": 2}[quality]
        cmd.set("surface_quality", {"high": 2, "medium": 1, "low": 0}[quality], obj)
        cmd.set("surface_best", spacing * (2 if quality == "high" else 1), obj)
        cmd.set("surface_normal", spacing, obj)
        cmd.set("solvent_radius", p["probe_radius"], obj)
        cmd.set("surface_solvent", 0, obj)
        cmd.set(
            "surface_type", {"surface-mesh": 2, "surface-dot": 1}.get(style, 0), obj
        )
        cmd.set("surface_color", -1, obj)
        cmd.show("surface", obj)
        if not p["show_hydrogens"]:
            cmd.hide("surface", f"{obj} and hydro")
    elif style == "gaussian-surface":
        custom = "gaussian"
    elif style == "nucleotides-atoms":
        atomic(cmd, obj, "stick", p)
    elif style == "nucleotides-fill":
        atomic(cmd, obj, "stick", p)
        custom = "rings"
    elif style.startswith("nucleotides-"):
        custom = "nucleotides"
        if style == "nucleotides-slab":
            atomic(cmd, obj, "stick", p)
            cmd.hide("sticks", f"{obj} and not name *'")
    elif style.startswith("aniso"):
        atomic(cmd, obj, "stick", {**p, "stick_radius": 0.055})
        custom = "aniso"
    elif style == "snfg":
        custom = "snfg"
    elif style in ("default", "ribbons-slabs", "cylinders-stubs", "licorice-ovals"):
        has_protein = cmd.count_atoms(f"{obj} and polymer.protein and name CA") > 1
        has_nucleic = cmd.count_atoms(f"{obj} and polymer.nucleic") > 0
        if not has_protein and not has_nucleic:
            atomic(cmd, obj, "stick", p)
        else:
            atomic(cmd, obj, "stick", p)
            cmd.hide("everything", f"{obj} and polymer")
            cmd.show("spheres", f"{obj} and inorganic")
            cmd.set("sphere_scale", 0.6, obj)
            custom = "preset"
    elif style in ("hbonds", "contacts", "struts"):
        atomic(cmd, obj, "stick", p)
        name = obj + "_links"
        if style == "struts":
            selection = f"{obj} and name CA+C5'"
            cmd.distance(
                name,
                selection,
                selection,
                cutoff=p["cutoff"] * 2,
                mode=0,
                label=0,
                quiet=1,
            )
        else:
            cmd.distance(
                name,
                obj,
                obj,
                cutoff=p["cutoff"],
                mode=2 if style == "hbonds" else 0,
                label=0,
                quiet=1,
            )
        cmd.set("dash_radius", p["pseudobond_radius"], name)
        cmd.set("dash_color", "forest" if style == "hbonds" else "orange", name)
        extras.append(name)
    elif style in ("distance", "angle", "torsion"):
        atomic(cmd, obj, "ball", p)
        ids = cmd.index(obj)
        count = {"distance": 2, "angle": 3, "torsion": 4}[style]
        chosen = data.get("indices", list(range(len(ids))))
        if (
            len(chosen) != count
            or len(set(chosen)) != count
            or any(i < 0 or i >= len(ids) for i in chosen)
        ):
            raise ValueError(
                f"{style} requires exactly {count} selected atoms or data.indices"
            )
        selections = [f"{obj} and index {ids[i][1]}" for i in chosen]
        name = obj + "_measurement"
        method = getattr(cmd, "dihedral" if style == "torsion" else style)
        method(name, *selections, quiet=1)
        cmd.set("label_color", "black", name)
        cmd.set("label_size", p["label_size"], name)
        extras.append(name)
    elif style == "label":
        atomic(cmd, obj, "ball", p)
        model = cmd.get_model(obj)
        labels = {}
        seen = set()
        for a in model.atom:
            key = (
                (a.segi, a.chain, a.resi)
                if p["label_level"] == "residue"
                else (a.segi, a.chain)
                if p["label_level"] == "chain"
                else a.index
            )
            if key in seen:
                continue
            seen.add(key)
            labels[a.index] = p["text"] or (
                a.name
                if p["label_level"] == "atom"
                else a.chain
                if p["label_level"] == "chain"
                else f"{a.resn} {a.resi}"
            )
        cmd.label(obj, f"{labels!r}.get(index, '')", quiet=1)
        cmd.set("label_size", p["label_size"], obj)
        cmd.set("label_color", "black", obj)
        cmd.set("label_outline_color", "white", obj)
    elif style == "unitcell":
        symmetry = cmd.get_symmetry(source)
        if symmetry is None:
            raise ValueError(
                "Unit-cell display requires crystallographic cell metadata"
            )
        cmd.set_symmetry(obj, *symmetry)
        cmd.set("cell_color", "gray30", obj)
        cmd.show("cell", obj)
    elif style in (
        "axis",
        "plane",
        "centroid",
        "inertia",
        "markers",
        "pseudobonds",
    ) or style.startswith("shape-"):
        custom = "shape"
    else:
        raise ValueError(f"Unknown molecular representation: {style}")
    if custom:
        name = obj + "_geometry"
        extras.append(name)
        for state in range(1, cmd.count_states(obj) + 1):
            model = cmd.get_model(obj, state=state)
            col = colors(model, color, cmd, model_index, atomic=atom_style)
            if custom == "rings":
                mesh = geometry.rings(
                    model,
                    col,
                    0.035 if style.endswith("thin") else 2 * p["stick_radius"],
                )
            elif custom == "cartoon":
                mesh = geometry.cartoon(model, col, style, p, detail)
            elif custom == "nucleotides":
                mesh = geometry.nucleotides(
                    model, col, style, p, detail, data.get("pairs")
                )
            elif custom == "aniso":
                mesh = geometry.aniso(model, col, style, p, detail)
            elif custom == "snfg":
                mesh = shapes.snfg(model, p, detail, data)
            elif custom == "gaussian":
                mesh = volume.gaussian(model, col, p)
            elif custom == "shape":
                mesh = shapes.geometry(
                    style, p, data, rgb(p["color"], cmd), detail, model
                )
            elif custom == "preset":
                cartoon_style = {
                    "cylinders-stubs": "helix-cylinder",
                    "licorice-ovals": "tube",
                }.get(style, "cartoon")
                mesh = geometry.cartoon(model, col, cartoon_style, p, detail)
                if has_nucleic:
                    nuc_style = {
                        "cylinders-stubs": "nucleotides-stubs",
                        "licorice-ovals": "nucleotides-ellipsoid",
                    }.get(style, "nucleotides-tube-slab")
                    mesh = geometry.nucleotides(
                        model, col, nuc_style, p, detail, data.get("pairs")
                    )
            mode = p["surface_style"] if custom == "gaussian" else "solid"
            add_mesh(cmd, mesh, name, state, p, transparency, mode)
    return [obj, *extras]
