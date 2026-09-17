# chimerax_style_in_pymol

ChimeraX-inspired molecular and volumetric representations in PyMOL.
The standalone `chimerax_style` command provides 84 named styles: atoms,
cartoons, surfaces, nucleotides, glycans, thermal ellipsoids, density maps,
annotations, shapes, presets, and lighting. Source atoms are preserved;
managed views support state playback, refresh/reset, and session save/load.

[Full guide](docs/pymol_chimerax.md) · [日本語](docs/ja/pymol_chimerax.md) ·
[GPU/ray comparison gallery](docs/gallery.md)

## Install

The reproducible environment supplies Python 3.10, PyMOL 3.1, and dependencies:

```sh
pixi install --locked
pixi run pymol
```

For an existing PyMOL environment (Python 3.10+):

```sh
uv pip install --python <pymol-python> .
```

Replace `<pymol-python>` with the interpreter used by PyMOL. PyMOL itself is
provided by conda/pixi, not pip. Runtime Python dependencies are NumPy, SciPy,
scikit-image, Pillow, and Gemmi. ChimeraX and the sibling CueMol/Mol* packages
are not required. Standard PyMOL commands are unchanged.

Register from PyMOL's Python console or `.pymolrc.py`:

```python
from chimerax_style_in_pymol import __init_plugin__
__init_plugin__()
```

Alternatively, with dependencies already installed, run `run chimerax_style.py`
from this checkout in the PyMOL command line. Importing the package alone does
not start PyMOL or change scene settings.

## Gallery

All 84 named styles, rendered in real PyMOL at **1200 × 900** with `quality=high`.
These previews show the interactive GPU output. The [full gallery](docs/gallery.md)
compares each with ray output and identifies its inputs. Protein examples use
[crambin, PDB 1CRN](https://www.rcsb.org/structure/1CRN); DNA and amino acids use
PyMOL builders. Map, glycan, shape, and displacement-tensor examples use clearly
identified synthetic data. Protein-only examples use a consistent camera.

### Atoms

<table>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/stick-gpu.png"><img src="docs/gallery/stick-gpu.png" alt="stick: PyMOL GPU rendering" width="260"></a><br><code>stick</code></td>
    <td width="33%" align="center"><a href="docs/gallery/ball-gpu.png"><img src="docs/gallery/ball-gpu.png" alt="ball: PyMOL GPU rendering" width="260"></a><br><code>ball</code></td>
    <td width="33%" align="center"><a href="docs/gallery/sphere-gpu.png"><img src="docs/gallery/sphere-gpu.png" alt="sphere: PyMOL GPU rendering" width="260"></a><br><code>sphere</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/ring-fill-gpu.png"><img src="docs/gallery/ring-fill-gpu.png" alt="ring-fill: PyMOL GPU rendering" width="260"></a><br><code>ring-fill</code></td>
    <td width="33%" align="center"><a href="docs/gallery/ring-fill-thin-gpu.png"><img src="docs/gallery/ring-fill-thin-gpu.png" alt="ring-fill-thin: PyMOL GPU rendering" width="260"></a><br><code>ring-fill-thin</code></td>
    <td width="33%"></td>
  </tr>
</table>

### Cartoons

<table>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/cartoon-gpu.png"><img src="docs/gallery/cartoon-gpu.png" alt="cartoon: PyMOL GPU rendering" width="260"></a><br><code>cartoon</code></td>
    <td width="33%" align="center"><a href="docs/gallery/cartoon-oval-gpu.png"><img src="docs/gallery/cartoon-oval-gpu.png" alt="cartoon-oval: PyMOL GPU rendering" width="260"></a><br><code>cartoon-oval</code></td>
    <td width="33%" align="center"><a href="docs/gallery/cartoon-rectangle-gpu.png"><img src="docs/gallery/cartoon-rectangle-gpu.png" alt="cartoon-rectangle: PyMOL GPU rendering" width="260"></a><br><code>cartoon-rectangle</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/cartoon-barbell-gpu.png"><img src="docs/gallery/cartoon-barbell-gpu.png" alt="cartoon-barbell: PyMOL GPU rendering" width="260"></a><br><code>cartoon-barbell</code></td>
    <td width="33%" align="center"><a href="docs/gallery/helix-tube-gpu.png"><img src="docs/gallery/helix-tube-gpu.png" alt="helix-tube: PyMOL GPU rendering" width="260"></a><br><code>helix-tube</code></td>
    <td width="33%" align="center"><a href="docs/gallery/helix-cylinder-gpu.png"><img src="docs/gallery/helix-cylinder-gpu.png" alt="helix-cylinder: PyMOL GPU rendering" width="260"></a><br><code>helix-cylinder</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/helix-wrap-gpu.png"><img src="docs/gallery/helix-wrap-gpu.png" alt="helix-wrap: PyMOL GPU rendering" width="260"></a><br><code>helix-wrap</code></td>
    <td width="33%" align="center"><a href="docs/gallery/tube-gpu.png"><img src="docs/gallery/tube-gpu.png" alt="tube: PyMOL GPU rendering" width="260"></a><br><code>tube</code></td>
    <td width="33%" align="center"><a href="docs/gallery/worm-gpu.png"><img src="docs/gallery/worm-gpu.png" alt="worm: PyMOL GPU rendering" width="260"></a><br><code>worm</code></td>
  </tr>
</table>

### Surfaces

<table>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/surface-gpu.png"><img src="docs/gallery/surface-gpu.png" alt="surface: PyMOL GPU rendering" width="260"></a><br><code>surface</code></td>
    <td width="33%" align="center"><a href="docs/gallery/surface-mesh-gpu.png"><img src="docs/gallery/surface-mesh-gpu.png" alt="surface-mesh: PyMOL GPU rendering" width="260"></a><br><code>surface-mesh</code></td>
    <td width="33%" align="center"><a href="docs/gallery/surface-dot-gpu.png"><img src="docs/gallery/surface-dot-gpu.png" alt="surface-dot: PyMOL GPU rendering" width="260"></a><br><code>surface-dot</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/gaussian-surface-gpu.png"><img src="docs/gallery/gaussian-surface-gpu.png" alt="gaussian-surface: PyMOL GPU rendering" width="260"></a><br><code>gaussian-surface</code></td>
    <td width="33%"></td>
    <td width="33%"></td>
  </tr>
</table>

### Nucleotides

<table>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/nucleotides-atoms-gpu.png"><img src="docs/gallery/nucleotides-atoms-gpu.png" alt="nucleotides-atoms: PyMOL GPU rendering" width="260"></a><br><code>nucleotides-atoms</code></td>
    <td width="33%" align="center"><a href="docs/gallery/nucleotides-fill-gpu.png"><img src="docs/gallery/nucleotides-fill-gpu.png" alt="nucleotides-fill: PyMOL GPU rendering" width="260"></a><br><code>nucleotides-fill</code></td>
    <td width="33%" align="center"><a href="docs/gallery/nucleotides-slab-gpu.png"><img src="docs/gallery/nucleotides-slab-gpu.png" alt="nucleotides-slab: PyMOL GPU rendering" width="260"></a><br><code>nucleotides-slab</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/nucleotides-tube-slab-gpu.png"><img src="docs/gallery/nucleotides-tube-slab-gpu.png" alt="nucleotides-tube-slab: PyMOL GPU rendering" width="260"></a><br><code>nucleotides-tube-slab</code></td>
    <td width="33%" align="center"><a href="docs/gallery/nucleotides-muffler-gpu.png"><img src="docs/gallery/nucleotides-muffler-gpu.png" alt="nucleotides-muffler: PyMOL GPU rendering" width="260"></a><br><code>nucleotides-muffler</code></td>
    <td width="33%" align="center"><a href="docs/gallery/nucleotides-ellipsoid-gpu.png"><img src="docs/gallery/nucleotides-ellipsoid-gpu.png" alt="nucleotides-ellipsoid: PyMOL GPU rendering" width="260"></a><br><code>nucleotides-ellipsoid</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/nucleotides-ladder-gpu.png"><img src="docs/gallery/nucleotides-ladder-gpu.png" alt="nucleotides-ladder: PyMOL GPU rendering" width="260"></a><br><code>nucleotides-ladder</code></td>
    <td width="33%" align="center"><a href="docs/gallery/nucleotides-stubs-gpu.png"><img src="docs/gallery/nucleotides-stubs-gpu.png" alt="nucleotides-stubs: PyMOL GPU rendering" width="260"></a><br><code>nucleotides-stubs</code></td>
    <td width="33%"></td>
  </tr>
</table>

### Special molecular views

<table>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/snfg-gpu.png"><img src="docs/gallery/snfg-gpu.png" alt="snfg: PyMOL GPU rendering" width="260"></a><br><code>snfg</code></td>
    <td width="33%" align="center"><a href="docs/gallery/aniso-gpu.png"><img src="docs/gallery/aniso-gpu.png" alt="aniso: PyMOL GPU rendering" width="260"></a><br><code>aniso</code></td>
    <td width="33%" align="center"><a href="docs/gallery/aniso-axes-gpu.png"><img src="docs/gallery/aniso-axes-gpu.png" alt="aniso-axes: PyMOL GPU rendering" width="260"></a><br><code>aniso-axes</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/aniso-ellipses-gpu.png"><img src="docs/gallery/aniso-ellipses-gpu.png" alt="aniso-ellipses: PyMOL GPU rendering" width="260"></a><br><code>aniso-ellipses</code></td>
    <td width="33%"></td>
    <td width="33%"></td>
  </tr>
</table>

### Maps

<table>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/volume-surface-gpu.png"><img src="docs/gallery/volume-surface-gpu.png" alt="volume-surface: PyMOL GPU rendering" width="260"></a><br><code>volume-surface</code></td>
    <td width="33%" align="center"><a href="docs/gallery/volume-mesh-gpu.png"><img src="docs/gallery/volume-mesh-gpu.png" alt="volume-mesh: PyMOL GPU rendering" width="260"></a><br><code>volume-mesh</code></td>
    <td width="33%" align="center"><a href="docs/gallery/volume-dot-gpu.png"><img src="docs/gallery/volume-dot-gpu.png" alt="volume-dot: PyMOL GPU rendering" width="260"></a><br><code>volume-dot</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/volume-image-gpu.png"><img src="docs/gallery/volume-image-gpu.png" alt="volume-image: PyMOL GPU rendering" width="260"></a><br><code>volume-image</code></td>
    <td width="33%" align="center"><a href="docs/gallery/volume-plane-gpu.png"><img src="docs/gallery/volume-plane-gpu.png" alt="volume-plane: PyMOL GPU rendering" width="260"></a><br><code>volume-plane</code></td>
    <td width="33%" align="center"><a href="docs/gallery/volume-orthoplanes-gpu.png"><img src="docs/gallery/volume-orthoplanes-gpu.png" alt="volume-orthoplanes: PyMOL GPU rendering" width="260"></a><br><code>volume-orthoplanes</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/volume-box-faces-gpu.png"><img src="docs/gallery/volume-box-faces-gpu.png" alt="volume-box-faces: PyMOL GPU rendering" width="260"></a><br><code>volume-box-faces</code></td>
    <td width="33%" align="center"><a href="docs/gallery/volume-slab-gpu.png"><img src="docs/gallery/volume-slab-gpu.png" alt="volume-slab: PyMOL GPU rendering" width="260"></a><br><code>volume-slab</code></td>
    <td width="33%" align="center"><a href="docs/gallery/volume-mip-gpu.png"><img src="docs/gallery/volume-mip-gpu.png" alt="volume-mip: PyMOL GPU rendering" width="260"></a><br><code>volume-mip</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/volume-topography-gpu.png"><img src="docs/gallery/volume-topography-gpu.png" alt="volume-topography: PyMOL GPU rendering" width="260"></a><br><code>volume-topography</code></td>
    <td width="33%" align="center"><a href="docs/gallery/volume-segment-gpu.png"><img src="docs/gallery/volume-segment-gpu.png" alt="volume-segment: PyMOL GPU rendering" width="260"></a><br><code>volume-segment</code></td>
    <td width="33%"></td>
  </tr>
</table>

### Annotations

<table>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/pseudobonds-gpu.png"><img src="docs/gallery/pseudobonds-gpu.png" alt="pseudobonds: PyMOL GPU rendering" width="260"></a><br><code>pseudobonds</code></td>
    <td width="33%" align="center"><a href="docs/gallery/distance-gpu.png"><img src="docs/gallery/distance-gpu.png" alt="distance: PyMOL GPU rendering" width="260"></a><br><code>distance</code></td>
    <td width="33%" align="center"><a href="docs/gallery/angle-gpu.png"><img src="docs/gallery/angle-gpu.png" alt="angle: PyMOL GPU rendering" width="260"></a><br><code>angle</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/torsion-gpu.png"><img src="docs/gallery/torsion-gpu.png" alt="torsion: PyMOL GPU rendering" width="260"></a><br><code>torsion</code></td>
    <td width="33%" align="center"><a href="docs/gallery/hbonds-gpu.png"><img src="docs/gallery/hbonds-gpu.png" alt="hbonds: PyMOL GPU rendering" width="260"></a><br><code>hbonds</code></td>
    <td width="33%" align="center"><a href="docs/gallery/contacts-gpu.png"><img src="docs/gallery/contacts-gpu.png" alt="contacts: PyMOL GPU rendering" width="260"></a><br><code>contacts</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/struts-gpu.png"><img src="docs/gallery/struts-gpu.png" alt="struts: PyMOL GPU rendering" width="260"></a><br><code>struts</code></td>
    <td width="33%" align="center"><a href="docs/gallery/label-gpu.png"><img src="docs/gallery/label-gpu.png" alt="label: PyMOL GPU rendering" width="260"></a><br><code>label</code></td>
    <td width="33%" align="center"><a href="docs/gallery/markers-gpu.png"><img src="docs/gallery/markers-gpu.png" alt="markers: PyMOL GPU rendering" width="260"></a><br><code>markers</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/axis-gpu.png"><img src="docs/gallery/axis-gpu.png" alt="axis: PyMOL GPU rendering" width="260"></a><br><code>axis</code></td>
    <td width="33%" align="center"><a href="docs/gallery/plane-gpu.png"><img src="docs/gallery/plane-gpu.png" alt="plane: PyMOL GPU rendering" width="260"></a><br><code>plane</code></td>
    <td width="33%" align="center"><a href="docs/gallery/centroid-gpu.png"><img src="docs/gallery/centroid-gpu.png" alt="centroid: PyMOL GPU rendering" width="260"></a><br><code>centroid</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/inertia-gpu.png"><img src="docs/gallery/inertia-gpu.png" alt="inertia: PyMOL GPU rendering" width="260"></a><br><code>inertia</code></td>
    <td width="33%" align="center"><a href="docs/gallery/unitcell-gpu.png"><img src="docs/gallery/unitcell-gpu.png" alt="unitcell: PyMOL GPU rendering" width="260"></a><br><code>unitcell</code></td>
    <td width="33%"></td>
  </tr>
</table>

### Shapes

<table>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/shape-sphere-gpu.png"><img src="docs/gallery/shape-sphere-gpu.png" alt="shape-sphere: PyMOL GPU rendering" width="260"></a><br><code>shape-sphere</code></td>
    <td width="33%" align="center"><a href="docs/gallery/shape-ellipsoid-gpu.png"><img src="docs/gallery/shape-ellipsoid-gpu.png" alt="shape-ellipsoid: PyMOL GPU rendering" width="260"></a><br><code>shape-ellipsoid</code></td>
    <td width="33%" align="center"><a href="docs/gallery/shape-cylinder-gpu.png"><img src="docs/gallery/shape-cylinder-gpu.png" alt="shape-cylinder: PyMOL GPU rendering" width="260"></a><br><code>shape-cylinder</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/shape-cone-gpu.png"><img src="docs/gallery/shape-cone-gpu.png" alt="shape-cone: PyMOL GPU rendering" width="260"></a><br><code>shape-cone</code></td>
    <td width="33%" align="center"><a href="docs/gallery/shape-box-gpu.png"><img src="docs/gallery/shape-box-gpu.png" alt="shape-box: PyMOL GPU rendering" width="260"></a><br><code>shape-box</code></td>
    <td width="33%" align="center"><a href="docs/gallery/shape-rectangle-gpu.png"><img src="docs/gallery/shape-rectangle-gpu.png" alt="shape-rectangle: PyMOL GPU rendering" width="260"></a><br><code>shape-rectangle</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/shape-triangle-gpu.png"><img src="docs/gallery/shape-triangle-gpu.png" alt="shape-triangle: PyMOL GPU rendering" width="260"></a><br><code>shape-triangle</code></td>
    <td width="33%" align="center"><a href="docs/gallery/shape-icosahedron-gpu.png"><img src="docs/gallery/shape-icosahedron-gpu.png" alt="shape-icosahedron: PyMOL GPU rendering" width="260"></a><br><code>shape-icosahedron</code></td>
    <td width="33%" align="center"><a href="docs/gallery/shape-dodecahedron-gpu.png"><img src="docs/gallery/shape-dodecahedron-gpu.png" alt="shape-dodecahedron: PyMOL GPU rendering" width="260"></a><br><code>shape-dodecahedron</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/shape-tube-gpu.png"><img src="docs/gallery/shape-tube-gpu.png" alt="shape-tube: PyMOL GPU rendering" width="260"></a><br><code>shape-tube</code></td>
    <td width="33%" align="center"><a href="docs/gallery/shape-ribbon-gpu.png"><img src="docs/gallery/shape-ribbon-gpu.png" alt="shape-ribbon: PyMOL GPU rendering" width="260"></a><br><code>shape-ribbon</code></td>
    <td width="33%" align="center"><a href="docs/gallery/shape-mesh-gpu.png"><img src="docs/gallery/shape-mesh-gpu.png" alt="shape-mesh: PyMOL GPU rendering" width="260"></a><br><code>shape-mesh</code></td>
  </tr>
</table>

### Presets

<table>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/default-gpu.png"><img src="docs/gallery/default-gpu.png" alt="default: PyMOL GPU rendering" width="260"></a><br><code>default</code></td>
    <td width="33%" align="center"><a href="docs/gallery/ribbons-slabs-gpu.png"><img src="docs/gallery/ribbons-slabs-gpu.png" alt="ribbons-slabs: PyMOL GPU rendering" width="260"></a><br><code>ribbons-slabs</code></td>
    <td width="33%" align="center"><a href="docs/gallery/cylinders-stubs-gpu.png"><img src="docs/gallery/cylinders-stubs-gpu.png" alt="cylinders-stubs: PyMOL GPU rendering" width="260"></a><br><code>cylinders-stubs</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/licorice-ovals-gpu.png"><img src="docs/gallery/licorice-ovals-gpu.png" alt="licorice-ovals: PyMOL GPU rendering" width="260"></a><br><code>licorice-ovals</code></td>
    <td width="33%" align="center"><a href="docs/gallery/space-filling-gpu.png"><img src="docs/gallery/space-filling-gpu.png" alt="space-filling: PyMOL GPU rendering" width="260"></a><br><code>space-filling</code></td>
    <td width="33%" align="center"><a href="docs/gallery/space-filling-single-gpu.png"><img src="docs/gallery/space-filling-single-gpu.png" alt="space-filling-single: PyMOL GPU rendering" width="260"></a><br><code>space-filling-single</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/surface-atomic-gpu.png"><img src="docs/gallery/surface-atomic-gpu.png" alt="surface-atomic: PyMOL GPU rendering" width="260"></a><br><code>surface-atomic</code></td>
    <td width="33%" align="center"><a href="docs/gallery/surface-chain-gpu.png"><img src="docs/gallery/surface-chain-gpu.png" alt="surface-chain: PyMOL GPU rendering" width="260"></a><br><code>surface-chain</code></td>
    <td width="33%" align="center"><a href="docs/gallery/ghostly-white-gpu.png"><img src="docs/gallery/ghostly-white-gpu.png" alt="ghostly-white: PyMOL GPU rendering" width="260"></a><br><code>ghostly-white</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/publication-gpu.png"><img src="docs/gallery/publication-gpu.png" alt="publication: PyMOL GPU rendering" width="260"></a><br><code>publication</code></td>
    <td width="33%" align="center"><a href="docs/gallery/publication-depth-gpu.png"><img src="docs/gallery/publication-depth-gpu.png" alt="publication-depth: PyMOL GPU rendering" width="260"></a><br><code>publication-depth</code></td>
    <td width="33%" align="center"><a href="docs/gallery/interactive-gpu.png"><img src="docs/gallery/interactive-gpu.png" alt="interactive: PyMOL GPU rendering" width="260"></a><br><code>interactive</code></td>
  </tr>
</table>

### Lighting

<table>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/simple-gpu.png"><img src="docs/gallery/simple-gpu.png" alt="simple: PyMOL GPU rendering" width="260"></a><br><code>simple</code></td>
    <td width="33%" align="center"><a href="docs/gallery/full-gpu.png"><img src="docs/gallery/full-gpu.png" alt="full: PyMOL GPU rendering" width="260"></a><br><code>full</code></td>
    <td width="33%" align="center"><a href="docs/gallery/soft-gpu.png"><img src="docs/gallery/soft-gpu.png" alt="soft: PyMOL GPU rendering" width="260"></a><br><code>soft</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="docs/gallery/gentle-gpu.png"><img src="docs/gallery/gentle-gpu.png" alt="gentle: PyMOL GPU rendering" width="260"></a><br><code>gentle</code></td>
    <td width="33%" align="center"><a href="docs/gallery/flat-gpu.png"><img src="docs/gallery/flat-gpu.png" alt="flat: PyMOL GPU rendering" width="260"></a><br><code>flat</code></td>
    <td width="33%"></td>
  </tr>
</table>

## Use

After loading a structure, use the PyMOL command line:

```text
chimerax_style
chimerax_style cartoon-barbell, selection=chain A, color=chain
chimerax_style full, representation=surface
chimerax_style nucleotides-tube-slab, color=nucleotide
chimerax_style aniso
chimerax_style ball, selection=organic, name=ligand, lighting=keep
chimerax_style volume-mesh, selection=density, params=contour.json
chimerax_style ray, filename=figure.png, width=2400, height=1800
chimerax_style refresh, name=all
chimerax_style reset, name=all
chimerax_style list
chimerax_style help
```

`selection` uses PyMOL syntax. Reusing `name` replaces a view after successful
preparation; distinct names compose views. `state=0` prepares all molecular
states. Source coordinates, bonds, colors and attributes remain unchanged.
`reset` restores source representations. `color=keep` retains source colors;
`lighting=keep` retains scene lighting. Global lighting affects unrelated
visible objects too. The existing background is preserved except by explicit
publication/background settings.

Additional geometry and scientific inputs use local JSON/NPZ/MRC/CCP4 files or
Python dictionaries. Unknown parameters and missing required data raise errors.
For example, `aniso` needs ANISOU tensors and ladder rungs need explicit base
pairs; absent data is not replaced with invented scientific annotations.

```python
from chimerax_style_in_pymol import chimerax_style

chimerax_style("cartoon", params={"width": 2.0, "thickness": 0.4})
chimerax_style("aniso", params={"probability": 0.5})
chimerax_style("volume-surface", data="density.npz", params={"level": 1.2})
```

## Fidelity and scope

The [representation matrix](docs/pymol_chimerax.md#representation-coverage)
covers ChimeraX's built-in atomic, cartoon, surface, nucleotide, SNFG, thermal,
map, marker, and shape display families. It lists every approximation and
input requirement. This is a PyMOL renderer, not a ChimeraX command interpreter
or an implementation of arbitrary Toolshed extensions and analysis tools.

Numeric defaults include 0.2 Å sticks, VDW-scaled 0.3 balls, 2.0 × 0.4 Å cartoons,
0.5 Å nucleotide slabs, and 0.45 Å rungs. Smooth highlights, pastel model colors,
element/base colors, and distinct lighting profiles aim for the ChimeraX look.
Spline fitting, VDW radii, surface tessellation, and shading differ.

`publication` silhouettes and shadow/occlusion effects are strongest in ray
output; native GPU shading differs. `volume-image` uses sampled transparent
planes and can show slicing artifacts. There is no pixel-equivalence claim or
comparison against a running ChimeraX renderer. See the paired gallery to
assess the actual outputs.

## Develop and regenerate

```sh
pixi install --locked
pixi run check
pixi run test
pixi run uv build
pixi run uv run --no-project --python .pixi/envs/default/bin/python python tests/render_gallery.py --gui --fetch
```

Gallery generation requires a desktop OpenGL context. `--fetch` downloads the
1CRN example only if missing. A local `--structure` can be used instead. Omit
`--gui` for headless ray-only generation, or use `--only cartoon surface` for
specific entries. Gallery rendering is local; CI runs lint, tests, and builds.

Published `docs/gallery/*.png` files are versioned so GitHub can show them.
Downloaded structures, caches, environments, build output, and temporary renders
are ignored. Geometry utilities adapted from the sibling projects retain their
MIT attribution in [NOTICE](NOTICE). Official ChimeraX references are listed in
the [guide](docs/pymol_chimerax.md#reference-and-verification).
