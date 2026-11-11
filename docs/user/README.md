# Ommatidia: user guide

Cover any mesh (a whole head, or just the eye region of a sculpt) with insect
ommatidia: a regular hexagonal lattice of domed corneal lenses with
interommatidial bristles, modelled on SEM images of the *Drosophila* eye.
Spacing is even *on the surface*, so facets never stretch, whatever the mesh's
UVs or triangle sizes.

Status: lattice, domed facets, bristles and an SEM-look material work.
Bristle sockets and corneal nipples are still to come.

## In Blender (Blender 4.2 or newer)

1. *Edit ▸ Preferences ▸ Get Extensions ▸ ⌄ ▸ Install from Disk…* and pick
   `ommatidia-<version>.zip`. To build the zip, run this from the repository
   root (it writes `dist/ommatidia-<version>.zip`):

   ```powershell
   # Windows PowerShell: Blender is usually not on PATH, so call it by full path
   New-Item -ItemType Directory -Force dist
   & "C:\Program Files\Blender Foundation\Blender 4.5\blender.exe" --factory-startup --command extension build --source-dir src/ommatidia --output-dir dist
   ```

   ```sh
   # macOS / Linux (or Windows with Blender on PATH)
   mkdir -p dist
   blender --factory-startup --command extension build --source-dir src/ommatidia --output-dir dist
   ```

   The output folder must exist first. `--factory-startup` keeps your other
   add-ons from loading (and printing noise) during the build.
2. Select the mesh to cover. Open the sidebar (**N**) ▸ **Ommatidia**.
3. Choose the density:
   - **Count** (default 750, a female *Drosophila* eye) works at any scale.
   - **Spacing** in real µm (default 16 µm) needs **Model scale**: Blender
     units per real metre. Example: a 400 µm eye sculpted 2 units wide has
     model scale 2 / 400e-6 = 5000. At model scale 1 the mesh is real-size.
4. **Generate Ommatidia.** A new object `<mesh>_Ommatidia` appears, parented
   to your mesh: one connected mesh of domed facets and bristles. (The lattice
   itself is kept, hidden, as `<mesh>_OmmatidiaLattice`; it stores the settings
   and the Diagnostics numbers.)
5. Now play with the settings. **Facets** and **bristle** shape settings change
   the result immediately. **Lattice** and density settings regenerate about
   half a second after you stop changing them (switch *Live update* off for
   very large eyes, and press Generate yourself).

### Where the lattice is regular

A hexagonal lattice cannot cover a curved surface without some 5- and 7-sided
facets (on a closed sphere at least twelve). The add-on grows a **perfect
patch** from a centre and pushes all those defects outward, into seams far
from it, where you can hide them (behind the head, at the edge of the eye).

- **Lattice ▸ Centre** picks where the patch starts:
  *Auto* (an open sheet: the point farthest from its edge; a closed mesh: the
  side facing −Y, i.e. toward Blender's front view), *Cursor* (the 3D cursor,
  placed with Shift + right-click on the mesh), *Object* (e.g. an Empty), or
  *Point*. **Use 3D Cursor** copies the cursor position into a fixed point.
- **Rotation** turns the rows of facets around the centre.
- How far it stays perfect depends on curvature: on a sphere, about 60° around
  the centre; on an eye-shaped cap, all the way to the last row or two at the
  rim. Cells shrink slightly (up to ~10%) toward the seams; that is the price
  of keeping them all hexagonal.
- **Seed** only changes the seams. **Irregularity** adds random jitter
  everywhere, because real eyes are not perfect either.

### Facets

- **Inflation**: dome height. 0 = flat, 1 = semi-dome (globular lenses).
- **Facet shape**: 0 = round domes (*Drosophila*), 1 = flat-topped hexagonal
  facets (many wasps and bees).
- **Pokeyness**: 0 = round lens, 1 = pointed cone.
- **Crease sharpness**: how narrow the groove between facets is.
- **Rim fade**: on an open sheet, facets flatten over this many spacings at the edge.
- **Resolution**: mesh detail per facet (−1 = automatic). Higher is smoother
  and heavier; it regenerates.

### Bristles

- **Bristles**: *Alternating* (every other facet corner: one bristle per
  ommatidium, as in *Drosophila*), *Random* (each corner with probability
  *Frequency*), or *None*.
- **Length / Radius** in real µm. In *Count* mode the size of your mesh is
  unknown, so one facet spacing counts as the *Spacing* value (16 µm): a 10 µm
  bristle is then about 0.6 facets long on any mesh.
- **Tilt**: angle from the surface normal; **Comb direction**: the direction
  bristles lean toward (world axes).

### Look

The ommatidia object gets an "Ommatidia SEM" material in both slots (facets,
bristles): surfaces seen edge-on glow, as in an electron micrograph. Use
*Material Preview* or *Rendered* view to see it. Your own materials in those
slots survive regeneration.

### Tips

- An eye made with *Sculpt ▸ Mask Extract* is usually a thin solid ("Extract
  as Solid" is on by default): two skins a hair apart. The add-on notices and
  covers only the outer skin (Diagnostics says so). Unticking "Extract as
  Solid" gives a single sheet, which is cleaner still.
- Smooth meshes work best. Hard creases (a raw cube) break some facet walls;
  add a Subdivision Surface or Bevel modifier. Diagnostics warns about it.
- On an open sheet the outermost half-facets at the rim are dropped, so you get
  ~10% fewer than the count. Closed meshes get a few % more (curvature packs
  the patch a little tighter).
- Same settings give the same result.

## Headless (batch) runs with main.py

```
blender -b input.blend --python main.py -- path/to/config.json
```

- `input.blend` should contain the mesh to cover. `main.py` uses the active
  mesh, or the only mesh in the file; set `TARGET_NAME` near the bottom of
  `main.py` to pick one by name.
- `config.json` is optional. It may contain any subset of the settings (see
  `src/ommatidia/config.py`); missing ones use the *Drosophila* defaults.
  Lengths are in metres, e.g. `{"density_mode": "SPACING", "spacing": 1.6e-05,
  "model_scale": 5000, "inflation": 0.6, "lattice_origin": [0, -0.0002, 0]}`.
  Unknown keys are rejected, so typos don't pass silently.
- Each run writes `00_results/<date-time>_seed<seed>/`:
  `scene.blend` (with the ommatidia and lattice objects), `config.json` (the
  full settings used) and `metrics.json` (count, spacing CV, hexagon fraction,
  defects, patches, bristles, timings).

What `main.py` does, in order: read the config → read the mesh as a surface →
grow the ommatidium lattice → derive facet outlines and bristle sites → build
the domed facets and bristles → write the objects and material → save the
results.

How the code is organised: [modules.md](modules.md).
