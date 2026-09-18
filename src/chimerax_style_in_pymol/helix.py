"""Fitted circular and straight helix axes with analytic optimization gradients."""

import numpy as np
from scipy.optimize import minimize

from .mesh import unit

# Published ChimeraX ideal-helix calibration coordinates (CA, angstroms).
_IDEAL = np.array(
    [
        [-4.543, -1.381, -5.088],
        [-4.871, -2.280, -1.408],
        [-1.332, -3.668, -1.376],
        [-0.007, -0.471, -2.952],
        [-1.782, 1.624, -0.322],
        [-0.267, -0.483, 2.456],
        [3.207, -0.031, 0.978],
        [2.712, 3.735, 0.827],
        [1.634, 3.794, 4.473],
    ]
)


def fit_helix(points, curved=True):
    """Use a line below 13 residues and a fitted circular arc for longer helices."""
    points = np.asarray(points, float)
    center = points.mean(axis=0)
    axis = np.linalg.svd(points - center, full_matrices=False)[2][0]
    if len(points) <= len(_IDEAL):
        ideal = _IDEAL[: len(points)]
        a, _, b = np.linalg.svd((ideal - ideal.mean(axis=0)).T @ (points - center))
        rotation = a @ np.diag([1, 1, np.linalg.det(a @ b)]) @ b
        center = (
            np.array([-0.395, -0.049, -0.215]) - ideal.mean(axis=0)
        ) @ rotation + center
        axis = np.array([0.613, 0.501, 0.610]) @ rotation
    guess = np.r_[center, axis]
    arc = False
    if curved and len(points) >= 13:
        middle = len(points) // 2
        centroids = np.array(
            [
                points[1:4].mean(axis=0),
                points[middle - 1 : middle + 2].mean(axis=0),
                points[-4:-1].mean(axis=0),
            ]
        )
        offsets = centroids[1:] - centroids[0]
        normal = np.cross(*offsets)
        if np.dot(normal, normal) >= 1e-8:
            axis = unit(normal)
            offset = np.linalg.solve(
                np.vstack([2 * offsets, axis]),
                np.r_[(offsets * offsets).sum(axis=1), 0],
            )
            center = centroids[0] + offset
            guess = np.r_[center, axis, np.linalg.norm(offset)]
            arc = True

    def objective(parameters):
        c, u = parameters[:3], unit(parameters[3:6])
        delta = points - c
        along = delta @ u
        radial = delta - along[:, None] * u
        lengths = np.maximum(np.linalg.norm(radial, axis=1), 1e-12)
        directions = radial / lengths[:, None]
        if arc:
            radius = parameters[6]
            residual = lengths - radius
            value = max(float(np.sqrt(np.sum(residual**2 + along**2))), 1e-12)
            center_gradient = (
                (radius - lengths)[:, None] * directions - along[:, None] * u
            ).sum(axis=0) / value
            axis_gradient = (delta * (along * radius / lengths)[:, None]).sum(
                axis=0
            ) / value
            radius_gradient = -residual.sum() / value
        else:
            residual = lengths - lengths.mean()
            value = float(np.mean(residual**2))
            center_gradient = -2 * (directions * residual[:, None]).mean(axis=0)
            axis_gradient = -2 * (delta * (along * residual / lengths)[:, None]).mean(
                axis=0
            )
        axis_gradient = (axis_gradient - u * np.dot(axis_gradient, u)) / max(
            np.linalg.norm(parameters[3:6]), 1e-12
        )
        gradient = np.r_[center_gradient, axis_gradient]
        if arc:
            gradient = np.r_[gradient, radius_gradient]
        return value, gradient

    fitted = minimize(
        objective, guess, jac=True, method="BFGS", tol=0.1, options={"maxiter": 5}
    ).x
    center, axis = fitted[:3], unit(fitted[3:6])
    delta = points - center
    if arc:
        radial = delta - (delta @ axis)[:, None] * axis
        direction = unit(radial)
        # Sort by unwrapped polar angle so projected residues cannot double back.
        basis = direction[0]
        angles = np.unwrap(
            np.arctan2(direction @ np.cross(axis, basis), direction @ basis)
        )
        centers = center + fitted[6] * direction[np.argsort(angles)]
        tangents = unit(np.cross(axis, centers - center))
    else:
        if np.dot(axis, points[-1] - points[0]) < 0:
            axis = -axis
        centers = center + np.sort(delta @ axis)[:, None] * axis
        tangents = np.tile(axis, (len(points), 1))
    radius = min(2.5, float(np.linalg.norm(points - centers, axis=1).mean()))
    return centers, tangents, max(radius, 0.05)
