"""ChimeraX-inspired representations; importing does not start or alter PyMOL."""

__version__ = "0.1.0"


def chimerax_style(
    style="default",
    selection="all",
    representation=None,
    color="auto",
    quality="high",
    name="chimerax",
    transparency=0,
    lighting="auto",
    state=0,
    data=None,
    params=None,
    filename="",
    width=1200,
    height=900,
    quiet=0,
    _self=None,
):
    """
    DESCRIPTION
        Apply ChimeraX-inspired molecular and map representations in PyMOL.

    USAGE
        chimerax_style [style [, selection [, representation [, color]]]]
        chimerax_style cartoon, selection=chain A
        chimerax_style soft, representation=surface
        chimerax_style nucleotides-tube-slab, color=nucleotide
        chimerax_style volume-mesh, selection=density, params=map.json
        chimerax_style list
        chimerax_style refresh, name=all
        chimerax_style reset, name=all
        chimerax_style select, selection=pk1
        chimerax_style ray, filename=figure.png, width=1600, height=1200

    ARGUMENTS
        style: named representation, preset, or lighting profile (see list).
        selection: a PyMOL atom selection or one loaded map object.
        representation: geometry override when applying a lighting profile.
        color: auto, keep, model, chain, element, nucleotide, secondary-structure,
            rainbow, bfactor, a PyMOL color name, or #rrggbb.
        quality: low, medium, high.
        name: managed view name; applying again replaces this view transactionally.
        transparency: fraction from 0 (opaque) to 1 (invisible).
        lighting: auto, keep, simple, full, soft, gentle, flat.
        state: 0 prepares every loaded state; a positive integer selects one state.
        data: local JSON/NPZ/MRC/CCP4 path or Python dictionary.
        params: local JSON path or Python dictionary of rendering parameters.

    NOTES
        Source atoms, colors, bonds, and coordinates are preserved. Native copies
        and CGOs work in the PyMOL viewport and standard ray export. Refresh after
        source edits. Reset restores source visibility and owned global lighting.
        See docs/pymol_chimerax.md for data schemas and approximation boundaries.
    """
    from .registry import GROUPS, STYLES, canonical, local_data, parameters

    style = canonical(style)
    if style == "help":
        print(chimerax_style.__doc__)
        return chimerax_style.__doc__
    if style == "list":
        for group, names in GROUPS.items():
            print(group + ": " + ", ".join(names))
        return dict(GROUPS)
    if _self is None:
        from pymol import cmd as _self
    from . import controller

    if style in ("reset", "refresh"):
        return getattr(controller, style)(_self, str(name))
    controller.maintenance(_self)
    if style in ("png", "ray"):
        return controller.export(
            _self, str(filename), int(width), int(height), style == "ray"
        )
    if style == "select":
        return controller.select_source(_self, str(selection))
    representation = canonical(representation) if representation else None
    if style not in STYLES or (
        representation is not None and representation not in STYLES
    ):
        raise ValueError("Unknown style or representation; use chimerax_style list")
    quality = str(quality).lower()
    if quality not in ("low", "medium", "high"):
        raise ValueError("quality must be low, medium, or high")
    transparency = float(transparency)
    if not 0 <= transparency <= 1:
        raise ValueError("transparency must be between 0 and 1")
    original_state = float(state)
    state = int(original_state)
    if state < 0 or state != original_state:
        raise ValueError("state must be a nonnegative integer")
    profile = str(lighting).lower()
    if style in GROUPS["Lighting"] and profile == "auto":
        profile = style
    options = dict(
        style=style,
        selection=str(selection),
        representation=representation,
        color=str(color),
        quality=quality,
        name=str(name),
        transparency=transparency,
        lighting=profile,
        state=state,
        data=local_data(data),
        params=parameters(params),
    )
    result = controller.apply(_self, options)
    if not int(quiet):
        print(f"ChimeraX-style {style}: {result}")
    return result


def __init_plugin__(app=None, _self=None):
    if _self is None:
        from pymol import cmd as _self
    _self.extend("chimerax_style", chimerax_style)
    from pymol import CmdException

    from .registry import STYLES

    # PyMOL's completion table is optional in embedded/headless instances.
    try:
        _self.auto_arg[0]["chimerax_style"] = [
            _self.Shortcut(
                [*STYLES, "list", "help", "reset", "refresh", "png", "ray", "select"]
            ),
            "style",
            ",",
        ]
    except (AttributeError, KeyError, CmdException):
        pass
