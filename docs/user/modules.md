# How the add-on is organised

```
                   main.py  (batch)        ui/  (Blender panel and button)
                        └──────────┬──────────┘
                                   ▼
   blender_io ──SURFACE──► lattice ──centres──► topology ──outlines──► facets ──domes, bristles──┐
   (reads your mesh;       (grows a perfect     (neighbours,            (the visible             │
    converts units)         hexagonal patch)     corners, walls)         composite mesh)         │
        ▲                                                                                       │
        └──── writes <mesh>_Ommatidia (visible) and <mesh>_OmmatidiaLattice (hidden) ◄──────────┘
                                   │
                    material: SEM look (corneal nipples later)
```

| Module | What it does | Details |
|---|---|---|
| `config.py` | Every setting, with defaults and checks; saved next to every result | [config.py](../../src/ommatidia/config.py) |
| `blender_io/` | Reads your mesh as a surface in real units; writes the lattice object | [blender_io.md](../../src/ommatidia/blender_io/blender_io.md) |
| `lattice/` | Places the ommatidium centres: a perfect patch, defects in the seams | [lattice.md](../../src/ommatidia/lattice/lattice.md) |
| `topology/` | Connects neighbours into facets; picks bristle corners | [topology.md](../../src/ommatidia/topology/topology.md) |
| `facets/` | Builds the domed facets and bristles, and reshapes them instantly | [facets.md](../../src/ommatidia/facets/facets.md) |
| `material/` | The SEM-look material | [material.md](../../src/ommatidia/material/material.md) |
| `ui/` | The sidebar panels, Generate button and live updates | [ui.md](../../src/ommatidia/ui/ui.md) |

`lattice/`, `topology/` and `facets/` are plain Python and numpy. They never
touch Blender, which is why they can be tested on ideal shapes (spheres,
ellipsoids, caps) outside Blender.
