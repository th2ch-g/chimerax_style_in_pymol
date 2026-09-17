"""Public representation vocabulary and validated rendering parameters."""

import json
import math
from pathlib import Path

GROUPS = {
    "Atoms": ("stick", "ball", "sphere", "ring-fill", "ring-fill-thin"),
    "Cartoons": (
        "cartoon",
        "cartoon-oval",
        "cartoon-rectangle",
        "cartoon-barbell",
        "helix-tube",
        "helix-cylinder",
        "helix-wrap",
        "tube",
        "worm",
    ),
    "Surfaces": ("surface", "surface-mesh", "surface-dot", "gaussian-surface"),
    "Nucleotides": (
        "nucleotides-atoms",
        "nucleotides-fill",
        "nucleotides-slab",
        "nucleotides-tube-slab",
        "nucleotides-muffler",
        "nucleotides-ellipsoid",
        "nucleotides-ladder",
        "nucleotides-stubs",
    ),
    "Special molecular views": ("snfg", "aniso", "aniso-axes", "aniso-ellipses"),
    "Maps": (
        "volume-surface",
        "volume-mesh",
        "volume-dot",
        "volume-image",
        "volume-plane",
        "volume-orthoplanes",
        "volume-box-faces",
        "volume-slab",
        "volume-mip",
        "volume-topography",
        "volume-segment",
    ),
    "Annotations": (
        "pseudobonds",
        "distance",
        "angle",
        "torsion",
        "hbonds",
        "contacts",
        "struts",
        "label",
        "markers",
        "axis",
        "plane",
        "centroid",
        "inertia",
        "unitcell",
    ),
    "Shapes": (
        "shape-sphere",
        "shape-ellipsoid",
        "shape-cylinder",
        "shape-cone",
        "shape-box",
        "shape-rectangle",
        "shape-triangle",
        "shape-icosahedron",
        "shape-dodecahedron",
        "shape-tube",
        "shape-ribbon",
        "shape-mesh",
    ),
    "Presets": (
        "default",
        "ribbons-slabs",
        "cylinders-stubs",
        "licorice-ovals",
        "space-filling",
        "space-filling-single",
        "surface-atomic",
        "surface-chain",
        "ghostly-white",
        "publication",
        "publication-depth",
        "interactive",
    ),
    "Lighting": ("simple", "full", "soft", "gentle", "flat"),
}
STYLES = tuple(s for names in GROUPS.values() for s in names)
ALIASES = {
    "ribbon": "cartoon",
    "sticks": "stick",
    "ballstick": "ball",
    "ball-and-stick": "ball",
    "spheres": "sphere",
    "spacefill": "sphere",
    "cpk": "sphere",
    "licorice": "tube",
    "putty": "worm",
    "nucleic": "nucleotides-tube-slab",
    "slab": "nucleotides-slab",
    "tube/slab": "nucleotides-tube-slab",
    "ladder": "nucleotides-ladder",
    "stubs": "nucleotides-stubs",
    "ellipsoid": "aniso",
    "carbohydrate": "snfg",
    "mesh": "surface-mesh",
    "dot": "surface-dot",
    "isosurface": "volume-surface",
    "volume": "volume-image",
    "slice": "volume-plane",
    "dihedral": "torsion",
    "original-look": "default",
    "silhouette": "publication",
}
DEFAULTS = {
    "stick_radius": 0.2,
    "ball_scale": 0.3,
    "width": 2.0,
    "thickness": 0.4,
    "coil_radius": 0.2,
    "arrow_scale": 2.0,
    "helix_radius": 2.0,
    "smooth": 1.0,
    "worm_min": 0.25,
    "worm_max": 2.0,
    "worm_attribute": "b",
    "xsection": "auto",
    "helix_mode": "default",
    "probe_radius": 1.4,
    "grid_spacing": 0.5,
    "resolution": 4.0,
    "level": 1.0,
    "slab_thickness": 0.5,
    "nucleotide_shape": "box",
    "rung_radius": 0.45,
    "sugar_radius": 0.2,
    "symbol_radius": 1.75,
    "aniso_scale": 1.0,
    "probability": None,
    "axis_radius": 0.025,
    "ellipse_radius": 0.02,
    "axis": "z",
    "position": 0.5,
    "slab_depth": 3,
    "slab_normal": [0, 0, 1],
    "segments": [1],
    "radius": 1.0,
    "height": 4.0,
    "size": [4, 3, 2],
    "center": [0, 0, 0],
    "color": "#6495ed",
    "dashes": 8,
    "pseudobond_radius": 0.075,
    "cutoff": 3.5,
    "label_size": 16.0,
    "label_level": "residue",
    "text": "",
    "show_hydrogens": False,
    "show_orientation": True,
    "transfer": None,
    "max_voxels": 8000000,
    "max_triangles": 2000000,
    "stride": 1,
    "volume_step": 0.2,
    "surface_style": "solid",
    "background": "keep",
}


def canonical(style):
    name = str(style).strip().lower().replace("_", "-")
    return ALIASES.get(name, name)


def local_data(value):
    if value is None or value == "":
        return {}
    if isinstance(value, dict):
        return dict(value)
    path = Path(str(value)).expanduser()
    if not path.is_file():
        raise ValueError(f"Local input file does not exist: {path.name}")
    if path.suffix.lower() == ".json":
        with path.open() as stream:
            result = json.load(stream)
        if not isinstance(result, dict):
            raise ValueError("JSON input must contain an object")
        return result
    return {"file": str(path)}


def parameters(value):
    extra = local_data(value)
    unknown = set(extra) - DEFAULTS.keys()
    if unknown:
        raise ValueError("Unknown parameters: " + ", ".join(sorted(unknown)))
    result = {**DEFAULTS, **extra}
    for key, default in DEFAULTS.items():
        v = result[key]
        if isinstance(default, bool):
            if not isinstance(v, bool):
                raise ValueError(f"{key} must be a JSON boolean")
        elif isinstance(default, (int, float)):
            if not isinstance(v, (int, float)) or not math.isfinite(v):
                raise ValueError(f"{key} must be a finite number")
            if key not in ("level", "position", "smooth", "dashes") and v <= 0:
                raise ValueError(f"{key} must be positive")
            if key == "dashes" and v < 0:
                raise ValueError("dashes must be nonnegative")
    for key in ("position", "smooth"):
        if not 0 <= result[key] <= 1:
            raise ValueError(f"{key} must lie between 0 and 1")
    for key in ("stride", "dashes", "slab_depth", "max_voxels", "max_triangles"):
        if int(result[key]) != result[key]:
            raise ValueError(f"{key} must be an integer")
        result[key] = int(result[key])
    if result["probability"] is not None:
        if not 0 < float(result["probability"]) < 1:
            raise ValueError("probability must lie strictly between 0 and 1")
    for key, choices in {
        "xsection": ("auto", "oval", "rectangle", "barbell"),
        "helix_mode": ("default", "tube", "cylinder", "wrap"),
        "nucleotide_shape": ("box", "muffler", "ellipsoid"),
        "axis": ("x", "y", "z"),
        "surface_style": ("solid", "mesh", "dot"),
        "label_level": ("atom", "residue", "chain"),
        "worm_attribute": ("b", "q", "partial_charge"),
    }.items():
        if result[key] not in choices:
            raise ValueError(f"{key} must be one of: {', '.join(choices)}")
    return result
