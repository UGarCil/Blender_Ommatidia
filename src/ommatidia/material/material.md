# module name: material

**purpose:** make renders of the ommatidia read like scanning electron
micrographs.

Code: [material.py](material.py) · used by: [pipeline.py](../pipeline.py)

## Introduction

In an SEM image, brightness comes mostly from the *edge effect*: a surface seen
at a grazing angle releases more secondary electrons than one seen head-on. So
facet centres look dark, their slopes and outlines glow, and thin bristles look
bright. Blender's Layer Weight node gives that ramp directly ("Facing"), so
the material is Layer Weight → Color Ramp (dark → bright) → Emission, mixed
with a little diffuse shading so shapes still read under lights.

The ommatidia object has two material slots: 0 for the facets and 1 for the
bristles. Both start with the "Ommatidia SEM" material. Slots the user has
filled are kept when the object is regenerated, and an existing "Ommatidia SEM"
material is reused, so edits to it survive too.

Sockets are linked by index because Blender 5 renamed some ("Fac" → "Factor").

Corneal nipples (a sub-micron bump pattern) are planned as a shader bump (M6).

## Atomic algorithms

- **sem_material** ([sem_look.py](atomic_algorithm_scripts/sem_look.py)): get or
  build the SEM material.

## Compound algorithm

`assign_materials(obj)`: make sure there are two slots; fill the empty ones with
`sem_material()`.
