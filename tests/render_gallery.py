"""Render every public style in an isolated PyMOL process and audit its images."""

import argparse
import hashlib
import json
import sys
from pathlib import Path
from time import monotonic, sleep
from urllib.request import urlopen

import numpy as np
from PIL import Image
from samples import caption, load

from chimerax_style_in_pymol import chimerax_style
from chimerax_style_in_pymol.registry import GROUPS, STYLES

ROOT = Path(__file__).resolve().parents[1]


def validate(path, width, height):
    with Image.open(path) as source:
        source.load()
        rgba = source.convert("RGBA")
    canvas = Image.new("RGBA", rgba.size, "white")
    canvas.alpha_composite(rgba)
    pixels = np.asarray(canvas.convert("RGB"))
    if pixels.shape != (height, width, 3):
        raise ValueError(f"Wrong image dimensions: {path.name}")
    mask = pixels.min(axis=2) < 238
    if np.count_nonzero(mask) < 150:
        raise ValueError(f"Empty image: {path.name}")
    border = np.r_[mask[0], mask[-1], mask[:, 0], mask[:, -1]]
    if np.mean(border) > 0.03:
        raise ValueError(f"Clipped image: {path.name}")
    canvas.convert("RGB").save(path, optimize=True)
    return {
        "visible_pixels": int(np.count_nonzero(mask)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def gallery_markdown():
    lines = [
        "# Gallery",
        "",
        "[README](../README.md) | [Representation and fidelity guide](pymol_chimerax.md)",
        "",
        "All 84 named styles rendered in real PyMOL at 1200 x 900. GPU and ray columns use identical cameras and geometry. These are PyMOL outputs, not ChimeraX screenshots. Lighting, shadows, outlines, and translucent volume compositing can differ between backends.",
        "",
    ]
    for group, styles in GROUPS.items():
        lines.extend(
            ["## " + group, "", "| Style / input | GPU | Ray |", "| --- | --- | --- |"]
        )
        for style in styles:
            lines.append(
                f"| `{style}`<br>{caption(style)} | ![{style} GPU](gallery/{style}-gpu.png) | ![{style} ray](gallery/{style}-ray.png) |"
            )
        lines.append("")
    (ROOT / "docs/gallery.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Capture both real OpenGL and ray output in an isolated Qt window",
    )
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="Download 1CRN to the ignored cache if it is missing",
    )
    parser.add_argument("--structure", type=Path, default=ROOT / ".cache/1CRN.pdb")
    parser.add_argument("--output", type=Path, default=ROOT / "docs/gallery")
    parser.add_argument("--only", nargs="+", choices=STYLES)
    parser.add_argument("--width", type=int, default=1200)
    parser.add_argument("--height", type=int, default=900)
    args = parser.parse_args()
    if args.fetch and not args.structure.is_file():
        args.structure.parent.mkdir(parents=True, exist_ok=True)
        with urlopen("https://files.rcsb.org/download/1CRN.pdb", timeout=60) as stream:
            args.structure.write_bytes(stream.read())
    if not args.structure.is_file():
        parser.error(
            "Provide a local --structure, or use --fetch for the documented 1CRN example"
        )
    args.output.mkdir(parents=True, exist_ok=True)
    styles = args.only or STYLES
    app = window = pm = None
    if args.gui:
        import pymol

        pymol.invocation.options.show_splash = 0
        from pmg_qt.pymol_qt_gui import PyMOLApplication, PyMOLQtGUI

        app = PyMOLApplication(["chimerax-gallery"])
        window = PyMOLQtGUI()
        window.show()
        cmd = window.pymolwidget.cmd
        window.resize(1300, 1050)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

        def pump():
            end = monotonic() + 0.15
            while monotonic() < end:
                app.processEvents()
                sleep(0.002)

        pump()
        for setting in (
            "internal_gui",
            "internal_feedback",
            "internal_prompt",
            "seq_view",
            "movie_panel",
        ):
            cmd.set(setting, 0)
        widget = window.pymolwidget
        ratio = widget.devicePixelRatioF()
        widget.setFixedSize(round(args.width / ratio), round(args.height / ratio))
        window.adjustSize()
        pump()
        with widget:
            widget.resizeGL(widget.width(), widget.height())
    else:
        import pymol2

        pm = pymol2.SingletonPyMOL()
        pm.start()
        cmd = pm.cmd
        cmd.viewport(args.width, args.height)

        def pump():
            pass

    reports = []
    try:
        for style in styles:
            started = monotonic()
            data, params = load(cmd, style, args.structure)
            color = "nucleotide" if style.startswith("nucleotides-") else "auto"
            if style in ("ribbons-slabs", "cylinders-stubs", "licorice-ovals"):
                cmd.fnab("ATGC", name="dna")
                cmd.translate([30, 0, 0], "dna")
            chimerax_style(
                style,
                data=data,
                params=params,
                color=color,
                quality="high",
                quiet=1,
                _self=cmd,
            )
            cmd.reset()
            if cmd.count_atoms("all"):
                cmd.orient("all")
            cmd.turn("y", 20)
            cmd.turn("z", -15)
            cmd.zoom("chimerax_chimerax", buffer=3, complete=1)
            if style == "unitcell":
                cmd.zoom("center", buffer=32)
            if style in (
                "shape-rectangle",
                "shape-triangle",
                "shape-mesh",
                "plane",
                "volume-plane",
                "volume-mip",
                "volume-topography",
            ):
                cmd.turn("x", 20)
            cmd.clip("slab", 200)
            pump()
            renderers = ("gpu", "ray") if args.gui else ("ray",)
            for renderer in renderers:
                path = args.output / f"{style}-{renderer}.png"
                cmd.draw(
                    args.width, args.height, antialias=2
                ) if renderer == "gpu" else None
                pump()
                chimerax_style(
                    "png" if renderer == "gpu" else "ray",
                    filename=path,
                    width=args.width,
                    height=args.height,
                    quiet=1,
                    _self=cmd,
                )
                try:
                    metrics = validate(path, args.width, args.height)
                except ValueError:
                    cmd.save(str(ROOT / ".cache/gallery-failure.pse"))
                    print(
                        {
                            "style": style,
                            "objects": cmd.get_names("all"),
                            "enabled": cmd.get_names("all", enabled_only=1),
                            "view": cmd.get_view(),
                            "state": cmd.get_state(),
                        },
                        flush=True,
                    )
                    cmd.png(
                        str(ROOT / ".cache/gallery-failure-ray.png"),
                        width=args.width,
                        height=args.height,
                        ray=1,
                    )
                    raise
                reports.append({"style": style, "renderer": renderer, **metrics})
            print(f"{style}: {monotonic() - started:.2f}s", flush=True)
        if args.gui and not args.only:
            gallery_markdown()
        audit = ROOT / ".cache/gallery-audit.json"
        audit.parent.mkdir(parents=True, exist_ok=True)
        audit.write_text(json.dumps(reports, indent=2) + "\n")
    finally:
        if pm:
            pm.stop()
        if window:
            window.hide()


if __name__ == "__main__":
    main()
