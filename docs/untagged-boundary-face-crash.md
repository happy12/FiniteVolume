# Crash: untagged boundary face causes an out-of-bounds `bcs[]` access

## Symptom

Running a real case (Euler equation, mesh loaded from a `.fvmesh` file with
`airfoil`/`inlet`/`outlet` boundary patches) crashes with:

```
.../bits/stl_vector.h:1272: std::vector<_Tp, _Alloc>::const_reference
std::vector<_Tp, _Alloc>::operator[](size_type) const
[with _Tp = EulerBoundaryCondition; ...]: Assertion '__n < this->size()' failed.
```

This happens deep into a run (hundreds/thousands of steps in), not at load
time, making it hard to diagnose from the crash site alone.

## Root cause

`Face::patch_id` defaults to `-1`, documented as "internal face"
(`include/UnstructuredMesh.h:31`).

When a mesh is loaded, only boundary edges that actually appear in the
`BOUNDARY` section (fvmesh) / as tagged line elements (gmsh) get their face's
`patch_id` set — see the matching loop in
`src/MeshReader.cpp:150-171` (`build_cells_faces_patches`). Any face that is a
genuine domain boundary (`cell_right == -1`) but whose edge never appears in
that tagged list is silently left at `patch_id == -1`.

`EulerFVMSolver::ghost_state()` then does an **unchecked** access:

```cpp
// src/EulerFVMSolver.cpp:175
const EulerBoundaryCondition& bc = bcs[face.patch_id];
```

`bcs` is a `std::vector<EulerBoundaryCondition>` sized to the number of
*known* patches (`main.cpp:333`, `run_euler()`). When `face.patch_id == -1`
and `size_type` is unsigned, `bcs[-1]` becomes `bcs[SIZE_MAX]` — the
out-of-bounds access that trips the assertion.

### Confirmed real-world trigger

The sibling `AirfoilMesher` project's own docs state its `.fvmesh` writer
never populates `top`/`bottom` boundary tags (only `airfoil`/`inlet`/`outlet`
are actually written) — this is a known, documented gap on that side. That
means a real structured/multi-block mesh from AirfoilMesher can contain
genuine domain-boundary faces that never get any patch tag at all, and this
crash was reproduced end-to-end with exactly such a mesh.

Meshes where every boundary edge is fully tagged (e.g. Bravo's Gmsh-generated
`.msh` meshes, which tag the entire outer loop) do not trigger this.

## Recommended fix

1. **Validate at load time, not deep into the run.** Right after
   `MeshReader::read()` succeeds (or inside `build_cells_faces_patches`
   itself), scan every face with `cell_right == -1` and confirm
   `patch_id >= 0`. If any aren't, fail immediately with a clear, actionable
   message, e.g.:

   ```
   Error: mesh has N boundary faces with no assigned patch (untagged
   boundary edges). Every boundary face must belong to a tagged patch.
   ```

   This mirrors the existing pattern in `MeshReader.cpp` for other malformed
   meshes (e.g. the out-of-range node index check a few lines above the
   `BOUNDARY`-parsing loop) — fail loudly and early, don't let a bad mesh
   silently reach the solver.

2. **Defense in depth in `ghost_state()`.** Replace the raw
   `bcs[face.patch_id]` with a bounds-checked access (e.g. `.at()`, or an
   explicit check with a clear error) so that even if some other code path
   ever bypasses the load-time validation, the failure mode is a message,
   not a `std::vector` assertion trace deep in the solve loop.

### Why "hard error" rather than defaulting to Wall

The codebase already has a precedent for a *named*-but-unconfigured patch:
it defaults to Wall and prints a warning (`main.cpp:347-348`), reasoning that
a wall can't leak/inject mass or energy through an unconfigured boundary.

That reasoning doesn't extend to a face with **no name at all** — there's no
way to know what that boundary is actually supposed to be (it could just as
easily be a farfield or symmetry plane as a wall), so silently picking Wall
risks quietly changing the physics without the user ever finding out. Failing
at load time, with a face count in the message, gives the user something
concrete to act on (fix the mesh generator's boundary tagging) instead of a
plausible-looking but possibly-wrong result.

## Not in scope for this note

Fixing AirfoilMesher's own `top`/`bottom` tagging gap is a separate task in
that project. This note only covers making FiniteVolume itself fail safely
and clearly when handed a mesh with untagged boundary faces, regardless of
which mesh generator produced it.
