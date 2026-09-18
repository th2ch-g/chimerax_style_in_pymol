# ChimeraX-inspired representations in PyMOL

[README](../README.md) | [Gallery](gallery.md) | [日本語](ja/pymol_chimerax.md)

The `chimerax_style` command provides 84 named representations, presets, and
lighting profiles. It runs in PyMOL 3.1 without ChimeraX or either sibling
package. It uses native molecular objects for atoms and solvent-excluded
surfaces, and independently constructed CGO meshes for other geometry.
Both the interactive viewport and PyMOL ray renderer can display the geometry.

## Interface

```text
chimerax_style [style], selection=all, representation=, color=auto, quality=high, name=chimerax, transparency=0, lighting=auto, state=0
```

`style` is a name in `chimerax_style list`. `representation` overrides geometry
when selecting a lighting profile. `selection` uses **PyMOL selection syntax**,
not ChimeraX model IDs or atom specifications. Maps accept one map object name.
`quality` is `low`, `medium`, or `high`. `transparency` is a fraction in [0, 1].
`state=0` prepares all molecular states; a positive state selects one. Lighting
is scene-wide: `lighting=keep` leaves the existing scene settings unchanged.
The background is preserved unless a publication preset or `params.background`
explicitly changes it.

```text
chimerax_style cartoon, selection=polymer, color=chain
chimerax_style soft, representation=surface, selection=chain A
chimerax_style nucleotides-tube-slab, selection=polymer.nucleic, color=nucleotide
chimerax_style ball, selection=organic, name=ligand, lighting=keep
chimerax_style volume-mesh, selection=density, params=contour.json
chimerax_style refresh, name=all
chimerax_style reset, name=all
chimerax_style ray, filename=figure.png, width=2400, height=1800
chimerax_style png, filename=viewport.png, width=1200, height=900
```

`data` and `params` accept a Python dictionary or a local file. Use JSON files
from the PyMOL command line because commas delimit command arguments. No style
downloads structures or annotations. The gallery generator's explicit `--fetch`
option is the only network operation in this package.

```python
from chimerax_style_in_pymol import chimerax_style

chimerax_style("cartoon-barbell", "chain A", params={
    "width": 2.0, "thickness": 0.4, "arrow_scale": 2.0,
})
chimerax_style("worm", params={"worm_attribute": "b", "worm_min": 0.25, "worm_max": 2.0})
chimerax_style("aniso", params={"probability": 0.5})
```

## Representation coverage

This covers the built-in atomic, cartoon, molecular-surface, nucleotide, SNFG,
thermal-ellipsoid, map, marker, and geometric-shape **display families**.
It is not a ChimeraX command interpreter. Arbitrary Toolshed extensions,
GUI tools, data processing commands, VR, medical image viewers, custom shaders,
and file readers are outside this display interface. For each family, the table
below identifies the implemented rendering and its limits; it does not claim
pixel-identical output or identical scientific analysis.

| ChimeraX family | `chimerax_style` names | Implementation / limits |
| --- | --- | --- |
| Atom/bond styles | `stick`, `ball`, `sphere` | Rounded sticks, VDW-scaled balls, full VDW spheres. Stick radius 0.2 Å and ball scale 0.3. PyMOL supplies bond perception and VDW radii. |
| Ring fill | `ring-fill`, `ring-fill-thin` | Bonded 3–6 member rings, including fused rings. Convex extrusions; use `stick` to remove fills. |
| Cartoon cross sections | `cartoon`, `cartoon-oval`, `cartoon-rectangle`, `cartoon-barbell` | Default oval helices/coils and rectangular strands/nucleic backbones. Width 2 Å, thickness 0.4 Å, arrow scale 2. Barbell profile and spline smoothing are independent approximations. |
| Helix modes | `helix-tube`, `helix-cylinder`, `helix-wrap` | Fitted circular arc, straight cylinder, or ribbon projected onto that fitted surface. Short helices use ideal-helix alignment; tube helices with at least 13 residues use arc fitting. Radius is fitted (maximum 2.5 Å); `params.helix_radius` overrides it. Spline joins, end caps, and triangulation remain approximations. |
| Uniform / attribute cartoons | `tube`, `worm` | Uniform tube or a radius interpolated from `b`, `q`, or `partial_charge`. Chain changes and large coordinate gaps split the path. |
| Molecular surfaces | `surface`, `surface-mesh`, `surface-dot` | Native PyMOL solvent-excluded surface, probe 1.4 Å; solid, mesh, and dot views of that surface. PyMOL's surface triangulation differs. |
| Gaussian envelope | `gaussian-surface` | Sum of atom-centered Gaussians; `resolution` is FWHM in Å. Absolute contour level, no scientific density normalization. |
| Nucleotide atoms / fill | `nucleotides-atoms`, `nucleotides-fill` | Native sticks, optionally with ring fill. |
| Nucleotide slabs | `nucleotides-slab`, `nucleotides-tube-slab`, `nucleotides-muffler`, `nucleotides-ellipsoid` | Base-aligned boxes, oval-section mufflers, or ellipsoids; slab mode retains ribose sticks/fill. Tube/slab uses connectors. Base thickness 0.5 Å; extents are fitted to observed base coordinates. |
| Ladder / half-rungs | `nucleotides-ladder`, `nucleotides-stubs` | Radius 0.45 Å. Explicit residue pairs create full rungs; unpaired residues remain stubs. Distance alone never establishes a base pair. |
| SNFG | `snfg` | All eight shape classes: sphere, cube, diamond, cone, star, rectangle, pentagon, hexagon. Common CCD names and covalent links are recognized; supply additional aliases via `data.sugars`. Unknown residues are not guessed. |
| Thermal ellipsoids | `aniso`, `aniso-axes`, `aniso-ellipses` | Eigenvectors and square roots of ANISOU eigenvalues. Default scale 1.0 is RMS displacement; optional probability uses a 3D chi-square quantile. Missing/invalid tensors raise errors. |
| Map contours | `volume-surface`, `volume-mesh`, `volume-dot` | Marching cubes with full affine coordinates. Levels are absolute values, not sigma unless the input already is sigma-normalized. |
| Map images | `volume-image`, `volume-plane`, `volume-orthoplanes`, `volume-box-faces`, `volume-slab`, `volume-mip` | Transfer-colored density planes, orthoplanes, six faces, tilted slab average, and axis-aligned maximum-intensity projection. `volume-image` samples a stack of translucent planes; it is not ChimeraX's volume ray marcher. Rotation can expose sampling artifacts; choose `axis` for the viewing direction. |
| Topography / segmentation | `volume-topography`, `volume-segment` | Height surface from one plane, or boundaries of explicit integer labels. No automatic segmentation or density interpretation. |
| Pseudobonds / measurements | `pseudobonds`, `distance`, `angle`, `torsion`, `hbonds`, `contacts`, `struts` | Explicit pairs or PyMOL measurements. `hbonds` uses PyMOL's polar-contact detector; `contacts` uses a distance cutoff; `struts` connects nearby backbone points. These do not reproduce ChimeraX's detection/strut optimization algorithms. |
| Labels / markers | `label`, `markers` | Native atom/residue/chain labels, spheres and explicit links. Scene labels, not a screen-layout/2D-label editor. |
| Principal geometry | `axis`, `plane`, `centroid`, `inertia`, `unitcell` | Principal axis/plane, centroid, covariance ellipsoid, and native crystallographic cell. Cell metadata is required. |
| General shapes | `shape-sphere`, `shape-ellipsoid`, `shape-cylinder`, `shape-cone`, `shape-box`, `shape-rectangle`, `shape-triangle`, `shape-icosahedron`, `shape-dodecahedron`, `shape-tube`, `shape-ribbon`, `shape-mesh` | Capped solids, coordinate-following tubes/ribbons, and explicit indexed meshes. |
| Composition presets | `default`, `ribbons-slabs`, `cylinders-stubs`, `licorice-ovals`, `space-filling`, `space-filling-single`, `surface-atomic`, `surface-chain`, `ghostly-white` | Cartoon/ligand compositions, paired protein/nucleotide depictions, space filling, or surfaces with preset color/opacity. Default chooses cartoon when backbone atoms exist, otherwise sticks. ChimeraX's chain-count heuristics are not duplicated. |
| Overall look | `publication`, `publication-depth`, `interactive` | White background with ray silhouettes, white depth-cued display, or depth-cued current background. PyMOL silhouettes are a ray effect; the GPU viewport has no matching screen-space silhouette pass. |
| Lighting | `simple`, `full`, `soft`, `gentle`, `flat` | Scene-wide PyMOL lighting approximations. Key/fill directions and smooth highlights target the ChimeraX appearance. PyMOL ray occlusion and shadows differ from ChimeraX's multi-direction shadows; the GPU viewport uses native PyMOL shading. |

Aliases include `ribbon`, `sticks`, `ball-and-stick`, `cpk`, `spacefill`, `putty`,
`nucleic`, `tube/slab`, `ladder`, `stubs`, `ellipsoid`, `carbohydrate`, `isosurface`,
`slice`, and `dihedral`. No `cuemol_style` alias is registered.

## Data schemas

Map data uses a finite XYZ array (axis order **x, y, z**, not z, y, x), with
`origin` and positive `spacing`, or a nonsingular affine `transform` mapping
voxel indices to Cartesian Å. NPZ input is loaded with `allow_pickle=False`.
MRC/CCP4 files and already loaded PyMOL map objects are also accepted.

```python
import numpy as np

np.savez("density.npz", values=values, origin=[0, 0, 0], spacing=[1, 1, 1])
chimerax_style("volume-surface", data="density.npz", params={"level": 1.2})
chimerax_style("volume-plane", data="density.npz", params={
    "axis": "z", "position": 0.5,
    "transfer": [[0, "#102040", 0], [1, "#4080ff", 0.1], [3, "#ffee80", 0.8]],
})
chimerax_style("volume-segment", data="labels.npz", params={"segments": [1, 2]})
```

Transfer entries are increasing `[value, color, opacity]` rows. Plane position
is a fraction from 0 to 1. `slab_normal` is a Cartesian vector. `slab_depth`
sets the sampled slab thickness in Å and sample count. `stride` subsamples
planes/dots/mesh edges. `volume_step` sets image-plane spacing in Å (default
0.2), with optical-depth correction when the step changes. `grid_spacing`
controls Gaussian surface sampling.
`max_voxels` and `max_triangles` bound map and geometry allocation; invalid
parameters are rejected before replacing the previous view.

```python
chimerax_style("nucleotides-ladder", data={"pairs": [[0, 7], [1, 6], [2, 5]]})
chimerax_style("snfg", data={"sugars": {"XYZ": ["sphere", "#00a651"]}})
chimerax_style("markers", data={
    "points": [[0, 0, 0], [2, 1, 0], [4, 0, 1]], "links": [[0, 1], [1, 2]],
})
chimerax_style("shape-mesh", data={
    "points": [[0, 0, 0], [2, 0, 0], [0, 2, 0]], "faces": [[0, 1, 2]],
})
chimerax_style("pseudobonds", data={
    "points": [[0, 0, 0], [3, 2, 1]], "pairs": [[0, 1]],
})
chimerax_style("distance", "ligand", data={"indices": [0, 5]})
```

Indices are zero-based. Nucleotide pairs refer to nucleotide residue order
within each selected molecular object. Point pairs index `data.points`;
measurement indices index the atoms of each selected object. A measurement
without explicit indices requires exactly 2, 3, or 4 selected atoms.
Shape dimensions use `center`, `size` (XYZ full extents), `radius`, and `height`.
Custom colors use `params.color` for shapes/maps and `color` for molecules.

## Lifecycle and sessions

Importing the module does not import/start PyMOL or alter settings. Register
with `__init_plugin__()`; standard PyMOL commands remain unchanged. Managed
views appear as `chimerax_<name>` groups. Source molecular objects remain in
the session with their original coordinates, properties, bonds and colors;
only selected source representations are hidden while a view owns them.
Reset restores their original visibility. Overlapping views retain ownership
until the last view is removed. Unrelated selections and objects are untouched.

Reapplying a name builds all replacement states first. Invalid input leaves the
previous view intact. `refresh` rebuilds from current source coordinates and
colors; automatic tracking of subsequent edits is not performed. All loaded
molecular states are prepared for native PyMOL frame playback. Map rendering
uses one selected map state. Display objects and plain session metadata survive
PSE save/load when this package is installed. Re-register the command after
starting a new PyMOL session.

Pick a native atom representation, then run `chimerax_style select,
selection=pk1` to select its original source atom. Custom CGO surfaces and
ribbons do not provide per-atom picking. If a source object is edited so that
its atom indices change, rebuild the view; reset avoids restoring visibility
onto a different atom with a reused index. Deleting a managed group is noticed
at the next style command; use `reset` for immediate cleanup.

Lighting affects all visible objects because PyMOL lighting is global. The
last applied view determines it. Resetting the final managed view restores
settings owned by this command, preserving subsequent user changes. GPU `png`
needs an active OpenGL window. Ray export also works headlessly. Export validates the rendered PNG before
atomically replacing the destination; a failed render preserves an existing file.

Custom solid, mesh, and dot CGO payloads are assembled in NumPy batches. This
reduces preparation overhead for large cartoons, nucleotide shapes, and maps;
triangle order, normals, colors, and opacity groups are preserved. All states
are still prepared before replacing a managed view.

## Reference and verification

The behavior and numeric defaults were checked against the official ChimeraX
documentation for [style](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/style.html),
[size](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/size.html),
[cartoon](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/cartoon.html),
[nucleotides](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/nucleotides.html),
[surface](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/surface.html),
[aniso](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/aniso.html),
[volume](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/volume.html),
[SNFG](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/snfg.html),
[lighting](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/lighting.html), and
[presets](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/preset.html).
Model, chain, element and nucleotide colors are checked against the official
[atomic color definitions](https://github.com/RBVI/ChimeraX/blob/8b1067a516a6d439537bcd355c3d1bdbea0f1292/src/bundles/atomic/src/colors.py).
Models begin tan, sky blue, plum and light green; chains use their case-insensitive
IDs rather than selection order. The element table includes all 109 published
entries and the gray fallback. Rainbow runs per residue within each chain;
B-factor uses blue-white-red and averages residues for cartoons. The optional
`secondary-structure` mode remains a plugin-specific red/gold/blue palette.

Helix fitting was compared numerically with the original `sse.py` on 1CRN/1GGG:
36 straight/curved fits have maximum center/radius errors of 0.000036/0.0000012 Å.
Analytic gradients avoid repeated finite-difference objective evaluations.
No ChimeraX application image comparison has been performed; native lighting,
spline joins, surfaces and other approximations above remain distinct.

The NumPy 2.5.3 audit renders all 84 named styles in GPU and ray. For 20,000
atoms, chain-color preparation takes 2.2 ms versus 73.6 ms before this update;
B-factor preparation takes 7.9 ms versus 48.9 ms. Chain lookup and B-factor
normalization no longer repeat linear scans per atom.

The sibling CueMol/Mol* projects informed package structure, managed-view
ownership, explicit scientific input handling, and paired GPU/ray galleries.
Generic mesh primitives were adapted under MIT; attribution is in [NOTICE](../NOTICE).
Tests exercise every style in real PyMOL, exact displacement-tensor geometry,
affine grids, source preservation, overlapping views, failed replacements,
state playback data, session round trips, parser registration, and PNG export.

Publication backgrounds use Python float components at the PyMOL setting boundary,
fixing rejected NumPy 2.x tuple representations. Fitted helix tubes are capped
separately from loop interpolation to prevent flared ends.
