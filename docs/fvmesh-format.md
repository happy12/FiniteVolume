# FVMESH 1.0 mesh format

A plain-text, purpose-built 2D mesh format for this solver, designed to
natively support control-volume cells with any number of faces (unlike
Gmsh's legacy triangle/quad-only element types). Read by
`MeshReader::read_fvmesh` ([src/MeshReader.cpp](../src/MeshReader.cpp)),
written by `FvMeshWriter::write` ([src/FvMeshWriter.cpp](../src/FvMeshWriter.cpp)).

Selected by the `.fvmesh` file extension (see `MeshReader::read`); `.msh`
files continue to go through the existing Gmsh reader.

## Grammar

```
# FVMESH 1.0
POINTS <n_nodes>
x0 y0
x1 y1
...

CELLS <n_cells>
<n_verts> v0 v1 ... v_{n_verts-1}
...

BOUNDARY <n_boundary_edges>
<node_a> <node_b> <patch_name>
...
```

- The first line must be `# FVMESH <major>.<minor>`, e.g. `# FVMESH 1.0`.
  It's version-checked, not exact-string-matched (the same way Gmsh's
  `$MeshFormat` version is): a version this reader doesn't recognize fails
  with a clear "unsupported FVMESH version" error rather than a generic
  "missing header" one. Only 1.0 exists today.
- **POINTS**: `n_nodes` lines follow, each `x y` (double precision, mesh
  length units). Node index is 0-based, in file order. `FvMeshWriter::write`
  takes a `precision` parameter (1-17 significant digits, default 6) for how
  these are formatted; `read_fvmesh` accepts any precision on input.
- **CELLS**: `n_cells` lines follow, each a control-volume cell as
  `n_verts v0 v1 ... v_{n_verts-1}` -- an ordered, counter-clockwise (CCW)
  list of node indices forming a simple polygon. `n_verts` may be any value
  >= 3; there is no fixed limit (a cell may have 3, 4, 20, or more faces).
- **BOUNDARY**: `n_boundary_edges` lines follow, each `node_a node_b
  patch_name`, naming an edge of the mesh (the node pair, in either order,
  must match exactly one edge of exactly one cell) and the boundary patch it
  belongs to. Patches are created on demand, keyed by `patch_name` --
  there is no separate numeric id/name mapping section, unlike Gmsh's
  physical groups. Interior edges (shared by two cells) are never listed
  here; they are derived automatically by matching shared node pairs
  between cells.
- Sections may appear in any order, but each must be complete (its declared
  count of lines) before the next section begins. Blank lines and `#`
  comment lines are ignored between sections.

## Non-goals

- No support for 3D volume cells (this is, and remains, a 2D solver).
- No attempt at compatibility with third-party VTK/Gmsh tooling -- this
  format is only meant to be read/written by this project and by an
  eventual in-house mesh generator targeting it.
- No embedded solver/boundary-condition parameters. Physics, boundary
  condition types/values, and run parameters are assigned separately by
  patch name in the case file (see [include/CaseInput.h](../include/CaseInput.h));
  the mesh file only supplies geometry and patch names.
