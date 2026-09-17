"""Portable source-tree loader: run chimerax_style.py in the PyMOL console."""

import sys
from pathlib import Path

source = str(Path(__file__).resolve().parent / "src")
if source not in sys.path:
    sys.path.insert(0, source)

from chimerax_style_in_pymol import __init_plugin__  # noqa: E402

__init_plugin__()
