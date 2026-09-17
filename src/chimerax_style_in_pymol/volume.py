"""Local scalar fields, molecular envelopes, contours, and sampled images."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree

from .colors import rgb
from .mesh import Mesh, merge, unit
from .primitives import indexed


@dataclass
class Grid:
    values: np.ndarray
    transform: np.ndarray

    def world(self, points):
        return np.asarray(points) @ self.transform[:3, :3].T + self.transform[:3, 3]


def read(data, budget):
    if "file" in data:
        path = Path(data["file"])
        if path.suffix.lower() == ".npz":
            with np.load(path, allow_pickle=False) as archive:
                data = {k: archive[k] for k in archive.files}
        elif path.suffix.lower() in (".mrc", ".map", ".ccp4"):
            import gemmi

            ccp4 = gemmi.read_ccp4_map(str(path), setup=False)
            sampling = np.array([ccp4.header_i32(k) for k in (8, 9, 10)], float)
            order = np.array([ccp4.header_i32(k) for k in (17, 18, 19)]) - 1
            if np.any(sampling <= 0) or sorted(order) != [0, 1, 2]:
                raise ValueError("MRC/CCP4 header has invalid lattice sampling or axes")
            starts = np.zeros(3)
            starts[order] = [ccp4.header_i32(k) for k in (5, 6, 7)]
            origin = np.array([ccp4.header_float(k) for k in (50, 51, 52)])
            ccp4.setup(0, gemmi.MapSetup.ReorderOnly)
            grid = ccp4.grid
            values = np.array(grid, copy=True)
            transform = np.eye(4)
            transform[:3, :3] = np.asarray(grid.unit_cell.orth.mat) / sampling[None, :]
            transform[:3, 3] = origin + transform[:3, :3] @ starts
            data = {"values": values, "transform": transform}
        else:
            raise ValueError("Maps must be local NPZ, MRC, MAP, or CCP4 files")
    if "values" not in data:
        raise ValueError(
            "A map requires data with values and origin/spacing or transform"
        )
    values = np.asarray(data["values"], np.float32)
    if values.ndim != 3 or min(values.shape) < 2 or not np.isfinite(values).all():
        raise ValueError(
            "Map values must be a finite 3D array with at least two samples per axis"
        )
    if values.size > budget:
        raise ValueError("Map exceeds max_voxels; crop or downsample the input")
    if "transform" in data:
        transform = np.asarray(data["transform"], float)
    else:
        transform = np.eye(4)
        transform[:3, :3] = np.diag(np.broadcast_to(data.get("spacing", 1.0), (3,)))
        transform[:3, 3] = data.get("origin", [0, 0, 0])
    if (
        transform.shape != (4, 4)
        or not np.isfinite(transform).all()
        or not np.allclose(transform[3], [0, 0, 0, 1])
        or abs(np.linalg.det(transform[:3, :3])) < 1e-10
    ):
        raise ValueError("Map transform must be a finite nonsingular affine 4x4 matrix")
    return Grid(values, transform)


def contour(grid, level, color):
    from skimage.measure import marching_cubes

    if not grid.values.min() < level < grid.values.max():
        raise ValueError(f"Contour level {level:g} must be inside the map data range")
    v, f, n, _ = marching_cubes(
        grid.values, level, allow_degenerate=False, gradient_direction="ascent"
    )
    normals = unit(n @ np.linalg.inv(grid.transform[:3, :3]))
    vertices = grid.world(v)
    # Enforce outward winding against the returned gradient normals.
    flip = (
        np.einsum(
            "ij,ij->i",
            np.cross(
                vertices[f[:, 1]] - vertices[f[:, 0]],
                vertices[f[:, 2]] - vertices[f[:, 0]],
            ),
            normals[f].mean(axis=1),
        )
        < 0
    )
    f[flip] = f[flip, ::-1]
    return Mesh(
        vertices, normals, np.tile(color, (len(v), 1)), f, np.zeros(len(v), int)
    )


def gaussian(model, colors, p):
    xyz = np.array([a.coord for a in model.atom])
    sigma = p["resolution"] / (2 * np.sqrt(2 * np.log(2)))
    spacing = p["grid_spacing"]
    low = xyz.min(axis=0) - 4 * sigma
    shape = np.ceil((xyz.max(axis=0) + 4 * sigma - low) / spacing).astype(int) + 1
    if np.prod(shape, dtype=np.int64) > p["max_voxels"]:
        raise ValueError("Gaussian grid exceeds max_voxels; increase grid_spacing")
    values = np.zeros(shape, np.float32)
    for a, center in zip(model.atom, xyz):
        if a.symbol == "H" and not p["show_hydrogens"]:
            continue
        lo = np.maximum(0, np.floor((center - 4 * sigma - low) / spacing).astype(int))
        hi = np.minimum(
            shape, np.ceil((center + 4 * sigma - low) / spacing).astype(int) + 1
        )
        axes = np.meshgrid(
            *(low[k] + np.arange(lo[k], hi[k]) * spacing for k in range(3)),
            indexing="ij",
            sparse=True,
        )
        d2 = sum((x - c) ** 2 for x, c in zip(axes, center))
        values[tuple(slice(a, b) for a, b in zip(lo, hi))] += np.exp(
            -d2 / (2 * sigma**2)
        )
    transform = np.diag([spacing] * 3 + [1.0])
    transform[:3, 3] = low
    mesh = contour(Grid(values, transform), p["level"], colors[0])
    _, ids = cKDTree(xyz).query(mesh.vertices)
    mesh.colors = colors[ids].astype(np.float32)
    mesh.owners = ids
    return mesh


def transfer(grid, p):
    if p["transfer"] is None:
        low, high = float(grid.values.min()), float(grid.values.max())
        if low == high:
            raise ValueError("Image rendering requires a nonconstant map")
        rows = np.array(
            [
                [low, 0.06, 0.1, 0.2, 0],
                [low + (high - low) * 0.25, 0.1, 0.3, 0.65, 0.01],
                [low + (high - low) * 0.6, 0.25, 0.7, 0.9, 0.15],
                [high, 1, 0.9, 0.5, 0.7],
            ]
        )
    else:
        rows = np.array([[float(v), *rgb(c), float(a)] for v, c, a in p["transfer"]])
    if (
        rows.ndim != 2
        or rows.shape[1] != 5
        or len(rows) < 2
        or not np.isfinite(rows).all()
        or np.any(np.diff(rows[:, 0]) <= 0)
        or np.any((rows[:, 1:] < 0) | (rows[:, 1:] > 1))
    ):
        raise ValueError(
            "transfer requires increasing [value, RGB color, opacity] entries"
        )
    return rows


def plane(grid, axis, position, p, projection=False, topography=False, alpha=False):
    other = [k for k in range(3) if k != axis]
    stride = p["stride"]
    a, b = np.meshgrid(
        np.arange(0, grid.values.shape[other[0]], stride),
        np.arange(0, grid.values.shape[other[1]], stride),
        indexing="ij",
    )
    pts = np.zeros(a.shape + (3,))
    pts[..., axis] = position
    pts[..., other[0]] = a
    pts[..., other[1]] = b
    if projection:
        values = grid.values.max(axis=axis)[::stride, ::stride]
    else:
        values = map_coordinates(grid.values, pts.reshape(-1, 3).T, order=1).reshape(
            a.shape
        )
    if topography:
        pts[..., axis] += (
            (values - values.min()) / max(float(np.ptp(values)), 1e-8) * p["height"]
        )
    rows = transfer(grid, p)
    colors = np.stack(
        [np.interp(values, rows[:, 0], rows[:, i]) for i in (1, 2, 3)], axis=-1
    )
    indices = np.arange(a.size).reshape(a.shape)
    i = indices[:-1, :-1].ravel()
    j = indices[1:, :-1].ravel()
    k = indices[1:, 1:].ravel()
    last = indices[:-1, 1:].ravel()
    mesh = indexed(
        grid.world(pts.reshape(-1, 3)),
        np.r_[np.c_[i, j, k], np.c_[i, k, last]],
        colors.reshape(-1, 3),
    )
    if alpha:
        mesh.alphas = np.interp(values.ravel(), rows[:, 0], rows[:, 4]).astype(
            np.float32
        )
    return mesh


def geometry(grid, style, p, color):
    if style in ("volume-surface", "volume-mesh", "volume-dot"):
        return contour(grid, p["level"], color)
    if style == "volume-segment":
        if not np.allclose(grid.values, np.round(grid.values)):
            raise ValueError("Segmentation requires an integer label map")
        pieces = []
        for k, label in enumerate(p["segments"]):
            mask = (grid.values == label).astype(np.float32)
            if mask.min() == mask.max():
                raise ValueError(f"Segment {label} has no boundary in the map")
            pieces.append(
                contour(Grid(mask, grid.transform), 0.5, np.roll(color, k % 3))
            )
        return merge(pieces)
    axis = "xyz".index(p["axis"])
    position = (grid.values.shape[axis] - 1) * p["position"]
    if style in ("volume-plane", "volume-mip", "volume-topography"):
        return plane(
            grid,
            axis,
            position,
            p,
            projection=style.endswith("mip"),
            topography=style.endswith("topography"),
        )
    if style == "volume-orthoplanes":
        return merge(
            [
                plane(grid, k, (grid.values.shape[k] - 1) * p["position"], p)
                for k in range(3)
            ]
        )
    if style == "volume-box-faces":
        return merge(
            [
                plane(grid, k, pos, p)
                for k in range(3)
                for pos in (0, grid.values.shape[k] - 1)
            ]
        )
    if style == "volume-slab":
        normal = unit(p["slab_normal"])
        if np.linalg.norm(normal) < 0.5:
            raise ValueError("slab_normal cannot be zero")
        side = unit(np.cross(normal, np.eye(3)[np.argmin(abs(normal))]))
        up = np.cross(normal, side)
        shape = np.array(grid.values.shape) - 1
        center = grid.world(shape * p["position"])
        extent = np.linalg.norm(grid.transform[:3, :3] @ shape) / 2
        count = min(180, max(grid.values.shape))
        u, v = np.meshgrid(
            np.linspace(-extent, extent, count),
            np.linspace(-extent, extent, count),
            indexing="ij",
        )
        world = center + u[..., None] * side + v[..., None] * up
        inv = np.linalg.inv(grid.transform)
        coords = (world - grid.transform[:3, 3]) @ inv[:3, :3].T
        grid_normal = normal @ inv[:3, :3].T
        values = np.mean(
            [
                map_coordinates(
                    grid.values,
                    (coords + grid_normal * offset).reshape(-1, 3).T,
                    order=1,
                    mode="constant",
                ).reshape(count, count)
                for offset in np.linspace(
                    -p["slab_depth"] / 2, p["slab_depth"] / 2, p["slab_depth"]
                )
            ],
            axis=0,
        )
        # Reuse the plane tessellator with an explicitly oriented local grid.
        affine = np.eye(4)
        affine[:3, 0] = side * 2 * extent / (count - 1)
        affine[:3, 1] = up * 2 * extent / (count - 1)
        affine[:3, 2] = normal
        affine[:3, 3] = center - extent * (side + up)
        return plane(Grid(np.repeat(values[:, :, None], 2, axis=2), affine), 2, 0, p)
    if style == "volume-image":
        spacing = np.linalg.norm(grid.transform[:3, axis])
        step = p["volume_step"] / spacing * p["stride"]
        positions = np.arange(0, grid.values.shape[axis] - 1 + 1e-6, step)
        other = [k for k in range(3) if k != axis]
        estimate = (
            2
            * len(positions)
            * np.prod([(grid.values.shape[k] - 1) // p["stride"] for k in other])
        )
        if estimate > p["max_triangles"]:
            raise ValueError(
                "Volume sampling exceeds max_triangles; increase stride or volume_step"
            )
        pieces = []
        for k in positions:
            mesh = plane(grid, axis, k, p, alpha=True)
            mesh.alphas = 1 - (1 - mesh.alphas) ** step
            pieces.append(mesh)
        return merge(pieces)
    raise ValueError(f"Unknown map representation: {style}")
