"""Geometry utilities adapted from the sibling MIT projects; see NOTICE."""

from dataclasses import dataclass

import numpy as np


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-12)


@dataclass
class Mesh:
    vertices: np.ndarray
    normals: np.ndarray
    colors: np.ndarray
    faces: np.ndarray
    owners: np.ndarray
    opacity: float = 1.0
    alphas: object = None

    def __post_init__(self):
        self.vertices = np.ascontiguousarray(self.vertices, dtype=np.float32).reshape(
            -1, 3
        )
        self.normals = np.ascontiguousarray(self.normals, dtype=np.float32).reshape(
            -1, 3
        )
        self.colors = np.ascontiguousarray(self.colors, dtype=np.float32).reshape(-1, 3)
        self.faces = np.ascontiguousarray(self.faces, dtype=np.uint32).reshape(-1, 3)
        self.owners = np.asarray(self.owners, dtype=np.int32)
        self.alphas = (
            np.ones(len(self.vertices), np.float32)
            if self.alphas is None
            else np.asarray(self.alphas, np.float32)
        )
        if (
            self.alphas.shape != (len(self.vertices),)
            or not np.isfinite(self.alphas).all()
            or np.any((self.alphas < 0) | (self.alphas > 1))
        ):
            raise ValueError("Invalid vertex opacity")
        if not np.isfinite(self.vertices).all() or not np.isfinite(self.normals).all():
            raise ValueError("Non-finite mesh coordinates or normals")

    @property
    def nbytes(self):
        return sum(
            a.nbytes
            for a in (
                self.vertices,
                self.normals,
                self.colors,
                self.faces,
                self.owners,
                self.alphas,
            )
        )

    def edge_data(self):
        """Weld section seams before finding boundary and crease edges."""
        _, welded = np.unique(np.round(self.vertices, 5), axis=0, return_inverse=True)
        faces = self.faces.astype(int)
        area = np.cross(
            self.vertices[faces[:, 1]] - self.vertices[faces[:, 0]],
            self.vertices[faces[:, 2]] - self.vertices[faces[:, 0]],
        )
        faces = faces[np.linalg.norm(area, axis=1) > 1e-8]
        normals = unit(
            np.cross(
                self.vertices[faces[:, 1]] - self.vertices[faces[:, 0]],
                self.vertices[faces[:, 2]] - self.vertices[faces[:, 0]],
            )
        )
        if not len(faces):
            return np.empty((0, 12), np.float32)
        edges = np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2)
        keys = np.sort(welded[edges], axis=1)
        order = np.lexsort((keys[:, 1], keys[:, 0]))
        keys, edges = keys[order], edges[order]
        fn = np.repeat(normals, 3, axis=0)[order]
        first = np.r_[0, np.flatnonzero(np.any(np.diff(keys, axis=0), axis=1)) + 1]
        last = np.r_[first[1:] - 1, len(keys) - 1]
        n0 = fn[first]
        n1 = np.where((last == first)[:, None], -n0, fn[last])
        keep = (np.sum(n0 * n1, axis=1) < 0.9999) & (keys[first, 0] != keys[first, 1])
        return np.ascontiguousarray(
            np.c_[
                self.vertices[edges[first, 0]], self.vertices[edges[first, 1]], n0, n1
            ][keep],
            dtype=np.float32,
        )


def merge(meshes, opacity=1.0):
    meshes = [m for m in meshes if len(m.faces)]
    if not meshes:
        return Mesh([], [], [], [], [], opacity)
    offsets = np.cumsum([0] + [len(m.vertices) for m in meshes[:-1]])
    return Mesh(
        np.concatenate([m.vertices for m in meshes]),
        np.concatenate([m.normals for m in meshes]),
        np.concatenate([m.colors for m in meshes]),
        np.concatenate([m.faces + o for m, o in zip(meshes, offsets)]),
        np.concatenate([m.owners for m in meshes]),
        opacity,
        np.concatenate([m.alphas * m.opacity for m in meshes]),
    )
