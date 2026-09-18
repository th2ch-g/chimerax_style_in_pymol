"""Independent cartoon, nucleotide, ring, and thermal-displacement geometry."""

from collections import defaultdict

import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline

from .helix import fit_helix
from .mesh import Mesh, merge, unit
from .primitives import convex, cylinder, sphere_template, sweep


def residues(model):
    groups = {}
    for i, a in enumerate(model.atom):
        groups.setdefault((a.segi, a.chain, a.resi, a.resn), {})[a.name] = i
    return list(groups.values())


def coordinates(model):
    return np.array([a.coord for a in model.atom], float).reshape(-1, 3)


def ellipsoid(center, radii, axes, color, owner=0, detail=24):
    v, f = sphere_template(detail)
    radii = np.maximum(np.asarray(radii), 1e-5)
    return Mesh(
        v * radii @ axes + center,
        unit(v / radii @ axes),
        np.tile(color, (len(v), 1)),
        f,
        np.full(len(v), owner),
    )


def polymer_runs(model):
    xyz = coordinates(model)
    runs, run = [], []
    previous = None
    for res in residues(model):
        idx = (
            res.get("CA")
            if "CA" in res and model.atom[res["CA"]].symbol == "C"
            else res.get("C5'")
        )
        if idx is None:
            continue
        a = model.atom[idx]
        key = (a.segi, a.chain, "CA" in res)
        if previous is not None and (
            key != previous
            or np.linalg.norm(xyz[idx] - xyz[run[-1][0]])
            > (4.8 if "CA" in res else 9.0)
        ):
            if len(run) > 1:
                runs.append(run)
            run = []
        run.append((idx, res))
        previous = key
    if len(run) > 1:
        runs.append(run)
    return runs


def cartoon(model, col, style, p, detail):
    xyz = coordinates(model)
    pieces = []
    for run in polymer_runs(model):
        ids = np.array([r[0] for r in run])
        points = xyz[ids].copy()
        ss = np.array([model.atom[i].ss if "CA" in r else "N" for i, r in run])
        hints = np.zeros_like(points)
        delta = np.diff(points, axis=0)
        for k, (idx, res) in enumerate(run):
            if ss[k] == "N" and "C1'" in res:
                hints[k] = xyz[res["C1'"]] - points[k]
            elif ss[k] == "H" and 0 < k < len(run) - 1:
                hints[k] = np.cross(delta[k - 1], delta[k])
            elif "O" in res and "C" in res:
                hints[k] = xyz[res["O"]] - xyz[res["C"]]
            else:
                hints[k] = [0, 1, 0]
        for k in range(1, len(hints)):
            if np.dot(hints[k], hints[k - 1]) < 0:
                hints[k] *= -1
        for k in range(1, len(points) - 1):
            if ss[k] == "S":
                points[k] = (1 - p["smooth"]) * points[k] + p["smooth"] * (
                    xyz[ids[k - 1]] + 2 * xyz[ids[k]] + xyz[ids[k + 1]]
                ) / 4
        mode = p["helix_mode"]
        if style.startswith("helix-"):
            mode = style[6:]
        override = (
            style.removeprefix("cartoon-")
            if style.startswith("cartoon-")
            else p["xsection"]
        )
        # Fit each helix independently and connect coil segments to its endpoints.
        fitted_helices = np.zeros(len(points), dtype=bool)
        helix = np.flatnonzero(ss == "H")
        for group in np.split(helix, np.flatnonzero(np.diff(helix) > 1) + 1):
            if len(group) < 3 or mode == "default":
                continue
            h = points[group]
            midline, axes, radius = fit_helix(h, curved=mode != "cylinder")
            radius = radius if p["helix_radius"] is None else p["helix_radius"]
            if mode == "wrap":
                points[group] = midline + unit(h - midline) * radius
                hints[group] = axes
            else:
                points[group] = midline
                fitted_helices[group] = True
                # Keep the fitted tube independent of the loop spline. Otherwise
                # loop curvature bends and flares the cylinder's terminal rings.
                spacing = np.linalg.norm(np.diff(midline, axis=0), axis=1).mean()
                knots = np.arange(len(group))
                samples = np.r_[
                    -0.3, np.arange(0, len(group) - 0.5, 0.5), len(group) - 0.7
                ]
                tube_path = CubicHermiteSpline(knots, midline, axes * spacing)(samples)
                tube_hints = np.tile(unit(h[0] - midline[0]), (len(samples), 1))
                nearest_helix = np.clip(
                    np.floor(samples + 0.5).astype(int), 0, len(group) - 1
                )
                tube_ids = ids[group[nearest_helix]]
                pieces.append(
                    sweep(
                        tube_path,
                        radius,
                        radius,
                        tube_hints,
                        col[tube_ids],
                        tube_ids,
                        "ellipse",
                        detail,
                    )
                )

        sample = np.linspace(
            0, len(points) - 1, (len(points) - 1) * max(8, detail // 2) + 1
        )
        path = CubicSpline(np.arange(len(points)), points, axis=0, bc_type="natural")(
            sample
        )
        hint = CubicSpline(
            np.arange(len(points)), unit(hints), axis=0, bc_type="natural"
        )(sample)
        nearest = np.clip(np.floor(sample + 0.5).astype(int), 0, len(points) - 1)
        kinds = ss[nearest]
        width = np.where(
            np.isin(kinds, ["H", "S", "N"]), p["width"] / 2, p["coil_radius"]
        )
        thick = np.full(len(sample), p["thickness"] / 2)
        thick[~np.isin(kinds, ["H", "S", "N"])] = p["coil_radius"]
        if mode in ("tube", "cylinder"):
            # A thin connecting path remains inside the separately capped tubes.
            mask = fitted_helices[nearest]
            width[mask] = thick[mask] = p["coil_radius"]
        if style == "tube":
            width[:] = thick[:] = p["coil_radius"] * 2.5
        if style == "worm":
            values = np.array(
                [getattr(model.atom[i], p["worm_attribute"]) for i in ids]
            )
            if not np.isfinite(values).all():
                raise ValueError("Worm attributes must be finite")
            mapped = p["worm_min"] + (values - values.min()) / max(
                float(np.ptp(values)), 1e-12
            ) * (p["worm_max"] - p["worm_min"])
            width = thick = np.interp(sample, np.arange(len(points)), mapped)
        if style not in ("tube", "worm"):
            for k in range(len(points)):
                if ss[k] == "S" and (k == len(points) - 1 or ss[k + 1] != "S"):
                    mask = (sample >= max(0, k - 0.85)) & (sample <= k + 0.1)
                    width[mask] = np.maximum(
                        0.025,
                        p["width"]
                        / 2
                        * p["arrow_scale"]
                        * (k + 0.1 - sample[mask])
                        / 0.95,
                    )
        sections = np.array(
            ["rectangle" if k in ("S", "N") else "ellipse" for k in kinds], object
        )
        if override != "auto":
            sections[np.isin(kinds, ["H", "S", "N"])] = {
                "oval": "ellipse",
                "rectangle": "rectangle",
                "barbell": "fancy",
            }[override]
        if style in ("tube", "worm"):
            sections[:] = "ellipse"
        if mode in ("tube", "cylinder"):
            sections[kinds == "H"] = "ellipse"
        # Adjacent sections share an endpoint to avoid cracks in mixed profiles.
        boundaries = np.r_[
            0, np.flatnonzero(sections[1:] != sections[:-1]) + 1, len(path) - 1
        ]
        for start, end in zip(boundaries[:-1], boundaries[1:]):
            sl = slice(start, end + 1)
            pieces.append(
                sweep(
                    path[sl],
                    width[sl],
                    thick[sl],
                    hint[sl],
                    col[ids[nearest[sl]]],
                    ids[nearest[sl]],
                    sections[start],
                    detail,
                )
            )
    if not pieces:
        raise ValueError("Cartoon requires at least two connected CA or C5' atoms")
    return merge(pieces)


def ring_cycles(model):
    graph = defaultdict(set)
    for bond in model.bond:
        i, j = bond.index
        if (
            model.atom[i].resi == model.atom[j].resi
            and model.atom[i].chain == model.atom[j].chain
        ):
            graph[i].add(j)
            graph[j].add(i)
    cycles = set()
    for start in graph:

        def walk(path):
            last = path[-1]
            for nxt in graph[last]:
                if nxt == start and 3 <= len(path) <= 6:
                    cycle = tuple(path)
                    cycles.add(min(cycle, (cycle[0], *cycle[:0:-1])))
                elif nxt > start and nxt not in path and len(path) < 6:
                    walk(path + [nxt])

        walk([start])
    return sorted(cycles)


def rings(model, col, thickness, sugar_only=False):
    xyz = coordinates(model)
    pieces = []
    for ids in ring_cycles(model):
        if sugar_only and not all("'" in model.atom[i].name for i in ids):
            continue
        points = xyz[list(ids)]
        normal = np.linalg.svd(points - points.mean(axis=0), full_matrices=False)[2][-1]
        pieces.append(
            convex(
                np.r_[points + normal * thickness / 2, points - normal * thickness / 2],
                col[ids[0]],
                ids[0],
            )
        )
    if not pieces:
        raise ValueError("No bonded rings of three to six atoms in the selection")
    return merge(pieces)


def nucleotides(model, col, style, p, detail, pairs=None):
    xyz = coordinates(model)
    pieces = [cartoon(model, col, "cartoon", p, detail)]
    entries = []
    for res in residues(model):
        if "C3'" not in res or "C1'" not in res:
            continue
        base = [
            i
            for name, i in res.items()
            if "'" not in name
            and name not in ("P", "OP1", "OP2", "OP3", "O1P", "O2P")
            and model.atom[i].symbol != "H"
        ]
        if len(base) < 3:
            continue
        anchor = xyz[res["C3'"]]
        end = xyz[res.get("N1" if "N9" in res else "N3", base[0])]
        entries.append((res, base, anchor, end))
    if not entries:
        raise ValueError("Nucleotide representations require ribose and base atoms")
    joined = set()
    if style.endswith("ladder") and pairs:
        # Pairs are explicit residue indices; no distance-only base-pair inference.
        for i, j in pairs:
            if not (0 <= i < len(entries) and 0 <= j < len(entries)) or i == j:
                raise ValueError("Nucleotide pair indices are outside the residue list")
            a, b = entries[i][2], entries[j][2]
            mid = (a + b) / 2
            pieces.extend(
                [
                    cylinder(
                        a, mid, p["rung_radius"], col[entries[i][1][0]], detail=detail
                    ),
                    cylinder(
                        mid, b, p["rung_radius"], col[entries[j][1][0]], detail=detail
                    ),
                ]
            )
            joined.update((i, j))
    for k, (res, base, anchor, end) in enumerate(entries):
        owner = base[0]
        color = col[owner]
        if style.endswith(("ladder", "stubs")):
            if k not in joined:
                pieces.append(
                    cylinder(anchor, end, p["rung_radius"], color, owner, detail)
                )
            continue
        points = xyz[base]
        center = points.mean(axis=0)
        normal = np.linalg.svd(points - center, full_matrices=False)[2][-1]
        along = unit(center - xyz[res["C1'"]])
        across = unit(np.cross(normal, along))
        along = unit(np.cross(across, normal))
        axes = np.array([along, across, normal])
        local = (points - center) @ axes.T
        radius = np.maximum(np.max(abs(local), axis=0), [1, 1, 0])
        radius[:2] += 0.2
        radius[2] = p["slab_thickness"] / 2
        shape = (
            "ellipsoid"
            if style.endswith("ellipsoid")
            else "muffler"
            if style.endswith("muffler")
            else p["nucleotide_shape"]
        )
        if shape == "ellipsoid":
            pieces.append(ellipsoid(center, radius, axes, color, owner, detail))
        elif shape == "muffler":
            pieces.append(
                sweep(
                    np.array([center - along * radius[0], center + along * radius[0]]),
                    radius[1],
                    radius[2],
                    np.array([across, across]),
                    [color, color],
                    [owner, owner],
                    "ellipse",
                    detail,
                )
            )
        else:
            vertices = np.array(
                [[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
            )
            pieces.append(convex(vertices * radius @ axes + center, color, owner))
        pieces.append(
            cylinder(anchor, xyz[res["C1'"]], p["sugar_radius"], color, owner, detail)
        )
        pieces.append(
            cylinder(xyz[res["C1'"]], center, p["sugar_radius"], color, owner, detail)
        )
        if p["show_orientation"]:
            # A small contrasting pip marks the positive face of the base plane.
            pieces.append(
                ellipsoid(
                    center + normal * (radius[2] + 0.025),
                    [0.18, 0.18, 0.06],
                    axes,
                    color * 0.45,
                    owner,
                    detail,
                )
            )
    if style == "nucleotides-slab":
        pieces.append(rings(model, col, 0.15, sugar_only=True))
    return merge(pieces)


def aniso(model, col, style, p, detail):
    pieces = []
    scale = p["aniso_scale"]
    if p["probability"] is not None:
        from scipy.stats import chi2

        scale *= np.sqrt(chi2.ppf(float(p["probability"]), 3))
    for i, a in enumerate(model.atom):
        u = getattr(a, "u_aniso", ())
        if len(u) != 6 or not np.any(u):
            continue
        xx, yy, zz, xy, xz, yz = u
        eigen, axes = np.linalg.eigh([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
        if not np.isfinite(eigen).all() or np.min(eigen) < -1e-7:
            raise ValueError("Anisotropic tensors must be positive semidefinite")
        radii = np.sqrt(np.maximum(eigen, 1e-8)) * scale
        axes = axes.T
        center = np.asarray(a.coord)
        if style == "aniso":
            pieces.append(ellipsoid(center, radii, axes, col[i], i, detail))
        elif style == "aniso-axes":
            for r, axis in zip(radii, axes):
                pieces.append(
                    cylinder(
                        center - r * axis,
                        center + r * axis,
                        p["axis_radius"],
                        col[i],
                        i,
                        detail,
                    )
                )
        else:
            t = np.linspace(0, 2 * np.pi, detail * 2 + 1)
            for j, k in ((0, 1), (0, 2), (1, 2)):
                path = (
                    center
                    + np.cos(t)[:, None] * radii[j] * axes[j]
                    + np.sin(t)[:, None] * radii[k] * axes[k]
                )
                pieces.extend(
                    cylinder(a, b, p["ellipse_radius"], col[i], i, 8)
                    for a, b in zip(path[:-1], path[1:])
                )
    if not pieces:
        raise ValueError(
            "No ANISOU tensors are present; isotropic B factors are not substituted"
        )
    return merge(pieces)


def cgo(mesh, mode="solid", radius=0.04, stride=1):
    from pymol import cgo as cg

    result = [cg.ALPHA, mesh.opacity]
    if mode == "mesh":
        edges = np.unique(
            np.sort(
                np.r_[
                    mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]]
                ],
                axis=1,
            ),
            axis=0,
        )
        pairs = edges[::stride]
        values = np.empty((len(pairs), 14), np.float32)
        values[:, 0] = cg.CYLINDER
        values[:, 1:4] = mesh.vertices[pairs[:, 0]]
        values[:, 4:7] = mesh.vertices[pairs[:, 1]]
        values[:, 7] = radius
        values[:, 8:11] = mesh.colors[pairs[:, 0]]
        values[:, 11:14] = mesh.colors[pairs[:, 1]]
        result.extend(values.ravel().tolist())
    elif mode == "dot":
        ids = np.arange(0, len(mesh.vertices), stride)
        values = np.empty((len(ids), 9), np.float32)
        values[:, 0] = cg.COLOR
        values[:, 1:4] = mesh.colors[ids]
        values[:, 4] = cg.SPHERE
        values[:, 5:8] = mesh.vertices[ids]
        values[:, 8] = radius
        result.extend(values.ravel().tolist())
    else:
        # PyMOL applies alpha to a BEGIN/END block, not individual vertices.
        alpha = np.round(mesh.alphas[mesh.faces].mean(axis=1) * mesh.opacity, 3)
        for opacity in np.unique(alpha):
            if opacity <= 0:
                continue
            result.extend([cg.ALPHA, float(opacity), cg.BEGIN, cg.TRIANGLES])
            ids = mesh.faces[alpha == opacity].ravel()
            values = np.empty((len(ids), 12), np.float32)
            values[:, 0] = cg.COLOR
            values[:, 1:4] = mesh.colors[ids]
            values[:, 4] = cg.NORMAL
            values[:, 5:8] = mesh.normals[ids]
            values[:, 8] = cg.VERTEX
            values[:, 9:12] = mesh.vertices[ids]
            result.extend(values.ravel().tolist())
            result.append(cg.END)
    return result
