"""Read PyMOL map lattice geometry without assuming an orthogonal unit cell."""

import numpy as np

from .volume import read


def from_pymol(cmd, name, state, budget):
    values = np.asarray(cmd.get_volume_field(name, state=state, copy=1))
    if values.size > budget:
        raise ValueError("Map exceeds max_voxels")
    session = cmd.get_session(name, partial=1, quiet=1)
    entry = next(n for n in session["names"] if n and n[0] == name)
    try:
        corners = np.array(entry[5][2][state - 1][6], float).reshape(8, 3)
    except (IndexError, TypeError, ValueError) as exc:
        raise ValueError(
            "Cannot read this PyMOL map's lattice; use a local MRC/CCP4/NPZ input"
        ) from exc
    transform = np.eye(4)
    transform[:3, :3] = np.column_stack(
        [corners[k] - corners[0] for k in (1, 2, 4)]
    ) / (np.array(values.shape) - 1)
    transform[:3, 3] = corners[0]
    matrix = cmd.get_object_matrix(name, state=state)
    if matrix is not None:
        transform = np.array(matrix).reshape(4, 4) @ transform
    return read({"values": values, "transform": transform}, budget)
