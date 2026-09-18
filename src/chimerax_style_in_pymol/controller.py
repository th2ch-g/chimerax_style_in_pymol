"""Transactional managed views, session persistence, and source restoration."""

import re
from copy import deepcopy
from pathlib import Path
from tempfile import NamedTemporaryFile
from uuid import uuid4

from . import lighting, render, shapes, volume
from .colors import rgb
from .registry import GROUPS


def storage(cmd):
    session = cmd._pymol.session
    if not hasattr(session, "chimerax_style_state"):
        session.chimerax_style_state = {"views": {}, "sources": {}, "settings": {}}
    return session.chimerax_style_state


def source_selection(cmd, selection, name):
    cmd.select(name, f"({selection}) and not _cxs_*", quiet=1)
    rows = []
    cmd.iterate(
        name, "rows.append((model,index,reps,ID,name,resi))", space={"rows": rows}
    )
    if not rows:
        cmd.delete(name)
        raise ValueError("The selection contains no source atoms")
    return rows


def release(cmd, state, entry):
    owned = set(entry["keys"])
    other = {key for e in state["views"].values() for key in e["keys"]}
    restore = {}
    for key in owned - other:
        if key in state["sources"]:
            restore[key] = state["sources"].pop(key)
    for obj in {key[0] for key in restore}:
        if obj not in cmd.get_names("objects"):
            continue
        # Index reuse after source edits must not restore another atom's visibility.
        cmd.alter(
            obj,
            "reps = previous((model,index), ID, name, resi, reps)",
            space={
                "previous": lambda key, aid, name, resi, reps: (
                    restore[key][0]
                    if key in restore
                    and tuple(restore[key][1:]) == (aid, name, resi)
                    and reps == 0
                    else reps
                )
            },
            quiet=1,
        )
        cmd.rebuild(obj)


def reset(cmd, name):
    state = storage(cmd)
    names = list(state["views"]) if name == "all" else [name]
    for key in names:
        entry = state["views"].pop(key, None)
        if entry is None:
            continue
        for obj in entry["objects"]:
            if obj in cmd.get_names("all"):
                cmd.delete(obj)
        if entry["group"] in cmd.get_names("all"):
            cmd.delete(entry["group"])
        release(cmd, state, entry)
    if not state["views"]:
        lighting.restore(cmd, state)


def maintenance(cmd):
    state = storage(cmd)
    existing = set(cmd.get_names("all"))
    for name, entry in list(state["views"].items()):
        if entry["group"] not in existing or not any(
            o in existing for o in entry["objects"]
        ):
            reset(cmd, name)


def apply(cmd, options):
    state = storage(cmd)
    name = options["name"]
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name) or name == "all":
        raise ValueError(
            "name must start with a letter and contain only letters, digits, or underscores; all is reserved"
        )
    style = options["representation"] or options["style"]
    profile = options["lighting"]
    color = options["color"]
    opacity = options["transparency"]
    p = options["params"]
    data = options["data"]
    if style in GROUPS["Lighting"]:
        if options["lighting"] == "auto":
            profile = style
        style = "default"
    if profile == "auto":
        profile = "simple"
    outline = style == "publication"
    depth = style in ("publication-depth", "interactive")
    background = p["background"]
    if style in ("publication", "publication-depth", "interactive"):
        style = "default"
    if outline or options["style"] == "publication-depth":
        background = "white"
    if style in ("space-filling", "space-filling-single"):
        color = "chain" if style == "space-filling" else color
        if style == "space-filling-single" and color == "auto":
            color = "model"
        if options["lighting"] == "auto":
            profile = "full" if style == "space-filling" else "flat"
        style = "sphere"
    if style in ("surface-atomic", "surface-chain", "ghostly-white"):
        color = {
            "surface-atomic": "element",
            "surface-chain": "chain",
            "ghostly-white": "white",
        }[style]
        if opacity == 0:
            opacity = {"surface-atomic": 0.7, "ghostly-white": 0.8}.get(style, 0)
        style = "surface"
    if profile not in (*lighting.PROFILES, "keep"):
        raise ValueError(
            "lighting must be auto, keep, simple, full, soft, gentle, or flat"
        )
    if background != "keep":
        rgb(background, cmd)
    token = "_cxs_" + uuid4().hex[:12]
    selection = token + "_selection"
    objects = []
    rows = []
    copies = {}
    original_view = cmd.get_view()
    committed = False
    try:
        if style in GROUPS["Maps"]:
            map_data = data
            if not data:
                obj = options["selection"]
                if cmd.get_type(obj) != "object:map":
                    raise ValueError(
                        "Map styles require a map object selection or local data"
                    )
                # Session map metadata retains grid axes, including skewed unit cells.
                from .mapio import from_pymol

                grid = from_pymol(cmd, obj, max(1, options["state"]), p["max_voxels"])
            else:
                grid = volume.read(map_data, p["max_voxels"])
            mesh = volume.geometry(grid, style, p, rgb(p["color"], cmd))
            obj = token + "_map"
            objects.append(obj)
            mode = (
                "mesh"
                if style == "volume-mesh"
                else "dot"
                if style == "volume-dot"
                else "solid"
            )
            render.add_mesh(cmd, mesh, obj, 1, p, opacity, mode)
            if style in (
                "volume-plane",
                "volume-orthoplanes",
                "volume-box-faces",
                "volume-slab",
                "volume-mip",
                "volume-image",
            ):
                cmd.set("cgo_lighting", 0, obj)
        elif style.startswith("shape-") or (
            style in ("markers", "pseudobonds", "axis", "plane", "centroid", "inertia")
            and "points" in data
        ):
            detail = {"low": 12, "medium": 20, "high": 32}[options["quality"]]
            mesh = shapes.geometry(style, p, data, rgb(p["color"], cmd), detail)
            obj = token + "_shape"
            objects.append(obj)
            render.add_mesh(cmd, mesh, obj, 1, p, opacity)
        else:
            rows = source_selection(cmd, options["selection"], selection)
            sources = list(dict.fromkeys(row[0] for row in rows))
            for mi, source in enumerate(sources):
                obj = f"{token}_{mi}"
                objects.append(obj)
                count = cmd.count_states(source)
                if options["state"] > count:
                    raise ValueError(
                        f"Requested state exceeds the {count} loaded states"
                    )
                cmd.create(
                    obj,
                    f"{selection} and %{source}",
                    source_state=options["state"],
                    target_state=0 if options["state"] == 0 else 1,
                    zoom=0,
                    quiet=1,
                )
                copies[obj] = [(row[0], row[1]) for row in rows if row[0] == source]
                created = render.molecular(
                    cmd,
                    obj,
                    source,
                    style,
                    color,
                    options["quality"],
                    p,
                    data,
                    opacity,
                    mi,
                )
                objects.extend(o for o in created if o != obj)
        group = "chimerax_" + name
        old = state["views"].get(name)
        if group in cmd.get_names("all") and (old is None or old["group"] != group):
            raise ValueError(
                f"Object {group} already exists and is not owned by this command"
            )
        # Commit only after every state and every input has rendered successfully.
        new = {
            "group": group,
            "objects": objects,
            "keys": [(r[0], r[1]) for r in rows],
            "copies": copies,
            "options": deepcopy(options),
        }
        for row in rows:
            state["sources"].setdefault((row[0], row[1]), tuple(row[2:]))
        state["views"][name] = new
        if old:
            for obj in old["objects"]:
                cmd.delete(obj)
            cmd.delete(old["group"])
            release(cmd, state, old)
        for obj in objects:
            cmd.group(group, obj, quiet=1)
        if rows:
            cmd.hide("everything", selection)
        lighting.apply(cmd, state, profile, background, outline, depth)
        committed = True
        return group
    except Exception:
        # Also remove objects created inside a renderer before it raised.
        for obj in cmd.get_names("all"):
            if obj.startswith(token):
                cmd.delete(obj)
        raise
    finally:
        cmd.delete(selection)
        restored_view = list(original_view)
        if committed:
            restored_view[17] = abs(restored_view[17]) * (
                1 if cmd.get_setting_int("orthoscopic") else -1
            )
        cmd.set_view(restored_view)


def refresh(cmd, name):
    state = storage(cmd)
    names = list(state["views"]) if name == "all" else [name]
    if any(n not in state["views"] for n in names):
        raise ValueError("No managed view with that name")
    return [apply(cmd, deepcopy(state["views"][n]["options"])) for n in names]


def select_source(cmd, selection, target="sele"):
    state = storage(cmd)
    selected = cmd.index(selection)
    original = []
    for obj, index in selected:
        for entry in state["views"].values():
            if obj in entry["copies"]:
                ids = [i for _, i in cmd.index(obj)]
                if index in ids:
                    original.append(entry["copies"][obj][ids.index(index)])
    if not original:
        raise ValueError("Select atoms in a managed native representation first")
    temporary = "_cxs_pick_" + uuid4().hex[:8]
    cmd.select(target, "none", quiet=1)
    try:
        for obj in dict(original):
            cmd.select_list(
                temporary,
                obj,
                [i for o, i in original if o == obj],
                mode="index",
                quiet=1,
            )
            cmd.select(target, temporary, merge=1, quiet=1)
    finally:
        cmd.delete(temporary)
    return len(original)


def export(cmd, filename, width, height, ray):
    if not filename or Path(filename).suffix.lower() != ".png":
        raise ValueError("filename must name a PNG output file")
    if width < 1 or height < 1:
        raise ValueError("Image dimensions must be positive integers")
    path = Path(filename).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=path.parent, suffix=".png", delete=False) as stream:
        scratch = Path(stream.name)
    try:
        cmd.png(str(scratch), width=width, height=height, ray=int(ray), quiet=1)
        cmd.sync()
        with scratch.open("rb") as stream:
            valid = stream.read(8) == b"\x89PNG\r\n\x1a\n"
        if not valid:
            raise RuntimeError(
                "PyMOL did not write a PNG; GPU export requires an active OpenGL window"
            )
        scratch.replace(path)
    finally:
        scratch.unlink(missing_ok=True)
    return str(path)
