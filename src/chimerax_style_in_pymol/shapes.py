"""Geometric annotations, molecular principal axes, and SNFG glycan symbols."""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial import ConvexHull

from .colors import rgb
from .geometry import coordinates, ellipsoid, residues
from .mesh import merge, unit
from .primitives import arrow, convex, cylinder, dashed, indexed, sphere, sweep

# Common CCD aliases. Additional symbols can be supplied explicitly in data.
SUGARS = {
    "GLC": ("sphere", "#0090bc"),
    "BGC": ("sphere", "#0090bc"),
    "MAN": ("sphere", "#00a651"),
    "BMA": ("sphere", "#00a651"),
    "GAL": ("sphere", "#ffd400"),
    "GLA": ("sphere", "#ffd400"),
    "NAG": ("box", "#0090bc"),
    "NDG": ("box", "#0090bc"),
    "NGA": ("box", "#ffd400"),
    "A2G": ("box", "#ffd400"),
    "BM3": ("box", "#00a651"),
    "FUC": ("cone", "#ed1c24"),
    "FUL": ("cone", "#ed1c24"),
    "RAM": ("cone", "#00a651"),
    "SIA": ("diamond", "#a54399"),
    "SLB": ("diamond", "#a54399"),
    "NGC": ("diamond", "#00a651"),
    "GCU": ("diamond", "#0090bc"),
    "IDR": ("diamond", "#a17a4d"),
    "XYS": ("star", "#f47920"),
    "XYP": ("star", "#f47920"),
    "ARA": ("star", "#00a651"),
    "ARB": ("star", "#00a651"),
    "FRU": ("pentagon", "#00a651"),
    "KDO": ("hexagon", "#ffd400"),
    "MUR": ("hexagon", "#a54399"),
    "TYV": ("rectangle", "#00a651"),
    "ABE": ("rectangle", "#f47920"),
}


def symbol(center, shape, radius, color, detail=24):
    if shape == "sphere":
        return sphere(center, radius, color, 0, detail)
    if shape == "cone":
        return cylinder(
            np.asarray(center) - [0, radius, 0],
            np.asarray(center) + [0, radius, 0],
            radius,
            color,
            detail=detail,
            radius_b=0.001,
        )
    if shape in ("box", "rectangle"):
        vertices = (
            np.array(
                [[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)], float
            )
            * 0.7
        )
        if shape == "rectangle":
            vertices[:, 1:] *= 0.5
    elif shape == "diamond":
        vertices = np.r_[np.eye(3), -np.eye(3)] * 1.3
    elif shape in ("star", "pentagon", "hexagon"):
        n = {"star": 10, "pentagon": 5, "hexagon": 6}[shape]
        t = np.arange(n) * 2 * np.pi / n
        r = np.where(np.arange(n) % 2, 0.42, 1) if shape == "star" else np.ones(n)
        rim = np.c_[r * np.cos(t), r * np.sin(t)]
        vertices = np.array([[x, y, z] for z in (-0.25, 0.25) for x, y in rim])
        if shape == "star":
            vertices = np.r_[vertices, [[0, 0, -0.25], [0, 0, 0.25]]]
            faces = []
            for i in range(n):
                j = (i + 1) % n
                faces.extend(
                    [
                        [2 * n, j, i],
                        [2 * n + 1, i + n, j + n],
                        [i, j, j + n],
                        [i, j + n, i + n],
                    ]
                )
            return indexed(vertices * radius + center, faces, color)
    else:
        raise ValueError(f"Unknown SNFG symbol: {shape}")
    return convex(vertices * radius + center, color)


def snfg(model, p, detail, data):
    xyz = coordinates(model)
    pieces = []
    centers = {}
    custom = data.get("sugars", {})
    for res in residues(model):
        ids = list(res.values())
        a = model.atom[ids[0]]
        info = custom.get(a.resn, SUGARS.get(a.resn))
        if not info:
            continue
        center = xyz[ids].mean(axis=0)
        pieces.append(symbol(center, info[0], p["symbol_radius"], rgb(info[1]), detail))
        centers.update({i: (center, ids[0]) for i in ids})
    links = set()
    for b in model.bond:
        i, j = b.index
        if i not in centers and j not in centers:
            continue
        a, ai = centers.get(i, (xyz[i], i))
        b, bi = centers.get(j, (xyz[j], j))
        key = tuple(sorted((ai, bi)))
        if ai != bi and key not in links:
            links.add(key)
            pieces.append(cylinder(a, b, 0.18, [0.55] * 3, detail=detail))
    if not pieces:
        raise ValueError(
            "No recognized carbohydrate residues; provide data.sugars for additional CCD aliases"
        )
    return merge(pieces)


def points(data, model=None):
    result = np.asarray(
        data.get("points", coordinates(model) if model is not None else []), float
    )
    if (
        result.ndim != 2
        or result.shape[1] != 3
        or not len(result)
        or not np.isfinite(result).all()
    ):
        raise ValueError(
            "data.points must be a nonempty array of finite XYZ coordinates"
        )
    return result


def geometry(style, p, data, color, detail, model=None):
    center = np.asarray(p["center"], float)
    size = np.asarray(p["size"], float)
    if (
        center.shape != (3,)
        or size.shape != (3,)
        or not np.isfinite(center).all()
        or not np.isfinite(size).all()
        or np.any(size <= 0)
    ):
        raise ValueError(
            "center and size require finite XYZ triples and positive sizes"
        )
    radius = p["radius"]
    if style in ("shape-sphere", "shape-ellipsoid"):
        return ellipsoid(
            center,
            np.full(3, radius) if style == "shape-sphere" else size / 2,
            np.eye(3),
            color,
            detail=detail,
        )
    if style in ("shape-cylinder", "shape-cone"):
        return cylinder(
            center - [0, p["height"] / 2, 0],
            center + [0, p["height"] / 2, 0],
            radius,
            color,
            detail=detail,
            radius_b=0.001 if style.endswith("cone") else None,
        )
    if style == "shape-box":
        return convex(
            center
            + np.array([[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)])
            * size
            / 2,
            color,
        )
    if style == "shape-rectangle":
        return indexed(
            center
            + np.array([[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]]) * size / 2,
            [[0, 1, 2], [0, 2, 3]],
            color,
        )
    if style in ("shape-icosahedron", "shape-dodecahedron"):
        phi = (1 + np.sqrt(5)) / 2
        v = np.array(
            [(0, a, b * phi) for a in (-1, 1) for b in (-1, 1)]
            + [(a, b * phi, 0) for a in (-1, 1) for b in (-1, 1)]
            + [(b * phi, 0, a) for a in (-1, 1) for b in (-1, 1)]
        )
        if style.endswith("dodecahedron"):
            v = v[ConvexHull(v).simplices].mean(axis=1)
        return convex(center + unit(v) * radius, color)
    if style in ("shape-triangle", "shape-mesh"):
        v = points(data)
        if style == "shape-triangle" and len(v) != 3:
            raise ValueError("A triangle requires exactly three points")
        faces = data.get("faces", [[0, 1, 2]] if style == "shape-triangle" else [])
        if not len(faces):
            raise ValueError("A custom mesh requires triangle faces")
        return indexed(v, faces, color)
    xyz = points(data, model)
    if style in ("shape-tube", "shape-ribbon"):
        if len(xyz) < 2 or np.any(np.linalg.norm(np.diff(xyz, axis=0), axis=1) < 1e-8):
            raise ValueError(
                "Tube/ribbon needs at least two distinct successive points"
            )
        path = CubicSpline(np.arange(len(xyz)), xyz, axis=0)(
            np.linspace(0, len(xyz) - 1, (len(xyz) - 1) * detail)
        )
        return sweep(
            path,
            np.full(len(path), radius),
            np.full(
                len(path), radius if style.endswith("tube") else p["thickness"] / 2
            ),
            np.tile([0, 1, 0], (len(path), 1)),
            np.tile(color, (len(path), 1)),
            np.zeros(len(path), int),
            "ellipse" if style.endswith("tube") else "rectangle",
            detail,
        )
    if style == "markers":
        pieces = [sphere(x, radius, color, 0, detail) for x in xyz]
        for i, j in data.get("links", []):
            if not 0 <= i < len(xyz) or not 0 <= j < len(xyz):
                raise ValueError("Marker link indices are out of bounds")
            pieces.append(
                cylinder(xyz[i], xyz[j], p["stick_radius"], color, detail=detail)
            )
        return merge(pieces)
    if style == "pseudobonds":
        pairs = data.get("pairs")
        if not pairs:
            raise ValueError("Pseudobonds require explicit data.pairs of point indices")
        pieces = []
        for i, j in pairs:
            if not 0 <= i < len(xyz) or not 0 <= j < len(xyz) or i == j:
                raise ValueError("Invalid pseudobond endpoints")
            pieces.append(
                dashed(
                    xyz[i],
                    xyz[j],
                    p["pseudobond_radius"],
                    color,
                    detail=8,
                    count=p["dashes"],
                )
            )
        return merge(pieces)
    center = xyz.mean(axis=0)
    if style == "centroid":
        return sphere(center, radius, color, 0, detail)
    if len(xyz) < 3:
        raise ValueError("Principal geometry requires at least three points")
    _, singular, axes = np.linalg.svd(xyz - center, full_matrices=False)
    lengths = np.maximum(np.max(abs((xyz - center) @ axes.T), axis=0), 0.1)
    if style == "axis":
        return arrow(
            center - lengths[0] * axes[0],
            center + lengths[0] * axes[0],
            p["stick_radius"],
            color,
            detail=detail,
        )
    if style == "plane":
        v = (
            np.array([[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]]) * lengths @ axes
            + center
        )
        return indexed(v, [[0, 1, 2], [0, 2, 3]], color)
    if style == "inertia":
        # Equal-weight covariance ellipsoid; not a mass-weighted inertia tensor.
        return ellipsoid(
            center,
            np.maximum(singular / np.sqrt(len(xyz)) * np.sqrt(5), 0.05),
            axes,
            color,
            detail=detail,
        )
    raise ValueError(f"Unknown shape: {style}")
