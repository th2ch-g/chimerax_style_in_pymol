"""Reversible PyMOL lighting approximations for ChimeraX's built-in presets."""

from .colors import rgb

BASE = {
    "orthoscopic": 1,
    "ambient": 0.4,
    "direct": 0.6,
    "reflect": 0.45,
    "specular": 0.5,
    "shininess": 30,
    "spec_power": 30,
    "spec_reflect": 0.8,
    "light_count": 2,
    "light": [-1, 1, 1],
    "light2": [0.2, 0.2, 0.959],
    "two_sided_lighting": 1,
    "depth_cue": 0,
    "ray_trace_fog": 0,
    "ray_shadows": 0,
    "ambient_occlusion_mode": 0,
    "ray_trace_mode": 0,
    "ray_opaque_background": 1,
    "antialias": 2,
    "gamma": 1.0,
    "transparency_mode": 1,
    "ray_max_passes": 128,
    "ray_transparency_specular": 0,
}
PROFILES = {
    "simple": {},
    "full": {
        "ray_shadows": 1,
        "ambient_occlusion_mode": 2,
        "ambient": 0.35,
        "direct": 0.65,
    },
    "soft": {
        "ray_shadows": 1,
        "ambient_occlusion_mode": 2,
        "ambient": 0.65,
        "direct": 0.2,
        "reflect": 0.15,
        "specular": 0.1,
    },
    "gentle": {
        "ray_shadows": 1,
        "ambient_occlusion_mode": 1,
        "ambient": 0.75,
        "direct": 0.18,
        "reflect": 0.1,
        "specular": 0.05,
    },
    "flat": {
        "ambient": 1,
        "direct": 0,
        "reflect": 0,
        "specular": 0,
        "ray_trace_mode": 1,
    },
}


def apply(cmd, storage, profile, background="keep", outline=False, depth=False):
    if profile == "keep":
        return
    if profile not in PROFILES:
        raise ValueError("lighting must be simple, full, soft, gentle, flat, or keep")
    settings = {**BASE, **PROFILES[profile]}
    if outline:
        settings["ray_trace_mode"] = 1
    if depth:
        settings.update(depth_cue=1, ray_trace_fog=1)
    if background != "keep":
        settings["bg_rgb"] = tuple(rgb(background, cmd))
    old = storage.setdefault("settings", {})
    for key, value in settings.items():
        if key not in old:
            old[key] = [cmd.get_setting_tuple(key), None]
        cmd.set(key, value, quiet=1)
        old[key][1] = cmd.get_setting_tuple(key)


def restore(cmd, storage):
    for key, (original, last) in storage.pop("settings", {}).items():
        if cmd.get_setting_tuple(key) == last:
            values = original[1]
            cmd.set(key, values if len(values) > 1 else values[0], quiet=1)
