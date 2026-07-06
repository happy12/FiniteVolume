# Wall diagnostics (skin friction, drag, boundary-layer thickness) -- implementation plan

Discussed 2026-07-05, following a factual audit confirming none of the
following exist in this codebase today: a skin-friction coefficient, any
force/drag integration, or a boundary-layer-thickness value reachable
outside the single `--verify-flat-plate` CLI harness (where `tau_wall`/
`u_tau`/`delta_99` are local variables in
[run_verify_flat_plate_boundary_layer()](../src/main.cpp) -- not returned,
not written anywhere, not part of any struct). Planning only -- no
implementation started or agreed to. This document is meant to be
self-contained enough to hand to an agent (or human) picking up the actual
work later.

## Scope

- **In scope**: per-patch skin-friction coefficient (`Cf`), pressure
  coefficient (`Cp`), raw lift/drag forces and their coefficients (`Cl`,
  `Cd`), a moment coefficient (`Cm`) about a case-file-specified reference
  point, wall-normal `y+`, and boundary-layer thickness estimates
  (`delta_99`, displacement thickness `delta*`, momentum thickness `theta`,
  shape factor `H = delta*/theta`) for `NavierStokesFVMSolver`, with the
  design kept reusable by `RANSFVMSolver`. Every per-location quantity
  (`Cf`, `Cp`, `y+`, and the boundary-layer thickness estimates) is reported
  **per wall mesh node**, not per boundary face/panel -- see "Face-to-node
  averaging" below for why and how.
- **Out of scope / explicit non-goals**:
  - 3D anything -- this is a 2D solver; drag/lift are per-unit-span
    quantities, consistent with the rest of the codebase's 2D framing (see
    [docs/MANUAL.md](MANUAL.md)'s "Known limitations").
  - Wiring `RANSFVMSolver` into `CaseInput`/`main.cpp`'s case-file/CLI path
    -- that is a separate, larger, pre-existing gap (see
    [docs/archive/rans-spalart-allmaras-tracker.md](archive/rans-spalart-allmaras-tracker.md)).
    This plan designs the wall-diagnostics computation as a shared,
    solver-agnostic free function specifically so `RANSFVMSolver` can adopt
    it with one call site once/if that wiring ever happens, but does not
    attempt that wiring itself.
  - A fully general, point-location/ray-casting boundary-layer profile
    extractor that works on an arbitrary (e.g. Delaunay-triangulated)
    near-wall mesh. See "Boundary-layer thickness" below for why, and the
    deferred alternative it's traded off against.

## Why a new shared module, not solver methods

`NavierStokesFVMSolver::step()` already computes per-face `tau_xx`/`tau_yy`/
`tau_xy` inline
([NavierStokesFVMSolver.cpp:178-183](../src/NavierStokesFVMSolver.cpp)), and
`RANSFVMSolver::step()` duplicates the identical pattern with `mu_eff` in
place of `mu`
([RANSFVMSolver.cpp:236-238](../src/RANSFVMSolver.cpp)) -- consistent with
this project's existing "one class per equation set, but share the
underlying numerics module" precedent (`GradientReconstruction.h`,
`WallDistance.h`). Neither retains that per-face value past the hot loop,
and recomputing it there would touch the hottest per-step OpenMP pass for
every internal face just to serve an occasional, wall-faces-only
diagnostic. Instead, wall diagnostics should be a **separate,
on-demand post-processing pass** over only the (small) set of `NoSlipWall`
boundary faces, following the exact precedent
`NavierStokesFVMSolver::compute_resolution_diagnostics()` already sets: it
also recomputes fresh gradients on demand rather than reusing anything
retained from `step()`
([NavierStokesFVMSolver.cpp:301-344](../src/NavierStokesFVMSolver.cpp)).

New header/source, no dependency on either solver class (same reasoning as
`WallDistance.h`'s "depends only on mesh geometry and a face list, not on
which solver produced the field" design):
[include/WallTraction.h](../include/WallTraction.h) /
[src/WallTraction.cpp](../src/WallTraction.cpp).

```cpp
struct WallFaceSample {
    int face_index;
    double x_mid, y_mid;         // face midpoint
    double p;                    // wall pressure (see "Wall pressure" below)
    double rho;                  // wall-adjacent cell's density (needed for y+)
    double effective_viscosity;  // mu (NS) or mu + rho*nu_t (RANS) at this face
    double tau_wall;             // wall-tangential shear stress magnitude (signed: positive = same sense as the tangent direction below)
    double tx, ty;                // unit tangent direction used for tau_wall's sign (rotate (nx, ny) by +90 deg)
    double y_wall_normal;         // first-cell wall-normal distance: face_normal_distance() from this face to its owning cell's centroid
};

// Computes the above for every NoSlipWall boundary face in 'wall_faces'
// (caller-supplied, e.g. WallDistance.h's own precedent for taking a face
// list rather than discovering walls itself). 'effective_viscosity' is a
// per-cell array so RANSFVMSolver can pass mu + rho*nu_t while
// NavierStokesFVMSolver passes a uniform mu. Deliberately knows nothing
// about case-file reference quantities -- Cf/Cp/y+ are all derived
// downstream from this struct's raw fields (see below), keeping this
// module as reference-agnostic as GradientReconstruction.h/WallDistance.h.
std::vector<WallFaceSample> compute_wall_traction(
    const UnstructuredMesh& mesh, const std::vector<int>& wall_faces,
    const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& p, const std::vector<double>& rho,
    const std::vector<double>& effective_viscosity,
    const std::vector<Gradient2>& grad_u, const std::vector<Gradient2>& grad_v,
    GradientScheme gradient_scheme);
```

Internally this reuses `face_gradient()` and a promoted (no longer
anonymous-namespace-local) `corrected_face_gradient_vector()` -- currently
private to `NavierStokesFVMSolver.cpp`
([NavierStokesFVMSolver.cpp:33-37](../src/NavierStokesFVMSolver.cpp)) -- to
get `tau_xx`/`tau_yy`/`tau_xy` at each wall face exactly as `step()` does,
then projects onto the tangent direction for `tau_wall` and reports the raw
tensor's normal-direction trace for pressure separately.

### Wall pressure

Neither solver stores a separate boundary pressure state -- Euler-family
walls are inviscid-mirror ghost states, not a reconstructed boundary value
(see `ghost_state()`,
[NavierStokesFVMSolver.cpp:279-298](../src/NavierStokesFVMSolver.cpp)).
Consistent with that existing convention, wall pressure for the force
integral is the owning cell's own pressure, `p(U[face.cell_left])` --
a first-order approximation, not a boundary-extrapolated one. This should
be stated plainly in the eventual code comment and MANUAL.md entry, the
same way this project documents every other first-order approximation
(e.g. `AdvectionDiffusionFVMSolver`'s Neumann-boundary gradient stencil
value).

## Skin friction coefficient

`Cf = tau_wall / (0.5 * rho_ref * V_ref^2)`, derived from `WallFaceSample`
at output time using the reference quantities below (not stored on the
struct itself, per the reference-agnostic design above). Reported per wall
node (see "Face-to-node averaging").

## Pressure coefficient (Cp)

`Cp = (p - p_ref) / (0.5 * rho_ref * V_ref^2)`, `p_ref` being the reference
static pressure -- add one more reference key, `reference_pressure`
(default: taken from the same `ns_farfield` state used for
`reference_density`/`reference_velocity_*` if unset, same fallback rule).
Cheapest of every quantity in this plan: no viscosity, no gradients, just
`WallFaceSample::p` and the reference values. Reported per wall node
alongside `Cf`.

## Wall-normal y+

`y+ = rho * u_tau * y_wall_normal / effective_viscosity`, where
`u_tau = sqrt(tau_wall / rho)` (friction velocity) and `y_wall_normal` is
the **first-cell wall-normal distance** already captured on
`WallFaceSample` -- `face_normal_distance()`
([GradientReconstruction.h](../include/GradientReconstruction.h)) from the
wall face to its owning cell's centroid, the same quantity
`compute_dt()`'s viscous stability term already uses
(Phase 6 of [docs/navier-stokes-tracker.md](navier-stokes-tracker.md)), so
no new geometry primitive is needed. Every term (`rho`, `tau_wall`,
`effective_viscosity`, `y_wall_normal`) is already on `WallFaceSample`,
so `y+` needs no reference quantities at all, unlike `Cf`/`Cp`.

This closes a real, previously-flagged gap: the RANS tracker's flat-plate
phase had to size its mesh for `y+ ~ 1` *a priori* via a Schlichting
estimate, with no way to check afterward whether the actual run landed
there (see
[docs/archive/rans-spalart-allmaras-tracker.md](archive/rans-spalart-allmaras-tracker.md)
Phase 4's "known setup difficulty"). Reported per wall node alongside
`Cf`/`Cp`.

## Drag, lift, and moment

Per named wall patch, sum over that patch's boundary faces:

```
friction_force = sum( face.area * tau_wall * (tx, ty) )
pressure_force = sum( face.area * -p * (nx, ny) )   // outward normal; force ON the wall from the fluid
total_force    = friction_force + pressure_force
```

then project `total_force` (and the friction-only/pressure-only parts) onto
the reference flow direction for drag, and onto its perpendicular for lift.
Report **per patch** (an airfoil patch and a tunnel-wall patch should not be
silently summed together) and a **domain total** across every `NoSlipWall`
patch. Nondimensionalize as `Cd = drag / (0.5 * rho_ref * V_ref^2 *
reference_length)`, `Cl` likewise.

**Moment**, about a new case-file-specified reference point
`(moment_reference_x, moment_reference_y)` (e.g. an airfoil's quarter-chord),
accumulated in the same per-face loop at near-zero extra cost once the
per-face force vector already exists:

```
r  = (face.x_mid - moment_reference_x, face.y_mid - moment_reference_y)
dF = face.area * (tau_wall * (tx, ty) + -p * (nx, ny))    // this face's own force contribution
M += r.x * dF.y - r.y * dF.x                              // 2D out-of-plane (z) moment, right-hand rule
```

reported per patch and as a domain total, same as force. Nondimensionalize
as `Cm = M / (0.5 * rho_ref * V_ref^2 * reference_length^2)`.

### New reference quantities (case-file keys)

Needed for both `Cf` and `Cd`/`Cl` -- dynamic pressure and a flow direction
don't exist as first-class concepts anywhere in this codebase today (the
closest is a `Farfield`/`ns_farfield` boundary's prescribed state, which is
already exactly "the freestream" for the common case of one farfield patch).

| Key | Default | Meaning |
|---|---|---|
| `reference_density` | auto | Freestream `rho` for dynamic pressure; if unset, taken from the first `ns_farfield` patch found |
| `reference_velocity_x`, `reference_velocity_y` | auto | Freestream velocity components; drag direction = this vector normalized; if unset, taken from the first `ns_farfield` patch |
| `reference_pressure` | auto | Freestream `p` for `Cp`; if unset, taken from the first `ns_farfield` patch found |
| `reference_length` | 1.0 | Length scale for `Cd`/`Cl`/`Cm` (e.g. chord); `1.0` degrades gracefully to a per-unit-span coefficient, documented as such |
| `moment_reference_x`, `moment_reference_y` | 0.0, 0.0 | Point `Cm` is taken about (e.g. an airfoil's quarter-chord) |

If no `ns_farfield` patch exists and these keys are also unset, fail loudly
at load time (`CaseInput::load` returning `false`) rather than silently
dividing by a zero dynamic pressure -- consistent with this project's
existing "fail loudly on ambiguous/missing config" precedent (e.g. the
untagged-boundary-face fail-fast recommendation in
[docs/untagged-boundary-face-crash.md](untagged-boundary-face-crash.md)).
If *multiple* `ns_farfield` patches with different states exist, warn and
use the first (in patch declaration order) rather than guessing which one
is "the" freestream.

## Boundary-layer thickness

This is the hard one, because it needs a wall-normal *profile* through the
mesh, not just a wall-face-local quantity -- fundamentally different data
access than the force integral above.

**Approach: cell-column marching, not point-location ray casting.**
`--verify-flat-plate`'s existing `delta_99` computation
([main.cpp](../src/main.cpp), inside `run_verify_flat_plate_boundary_layer()`)
implicitly exploits that its generator mesh
(`build_flat_plate_mesh()`/`build_structured_verification_mesh()`) has cells
stacked in a known wall-normal column order, and simply walks that known
index structure. Generalizing this to *any* mesh means walking the mesh's
actual face-adjacency graph instead of a known index scheme:

1. Start at a `NoSlipWall` boundary face's owning cell.
2. Among that cell's other faces (excluding the one just crossed), pick the
   neighbor across whichever face's outward normal is most nearly parallel
   to the current wall-normal marching direction (initialized to the wall
   face's own outward normal). This is the same "most-aligned-neighbor"
   idea `face_normal_distance()` already uses for a single step
   ([GradientReconstruction.h](../include/GradientReconstruction.h)),
   repeated cell-to-cell.
3. Accumulate wall-normal distance via `face_normal_distance()` at each
   step; sample the tangential velocity component (tangent = wall face's
   own tangent, held fixed for the whole march, not re-derived per cell) at
   each visited cell centroid.
4. Stop after a max cell count / max distance, or once the sampled velocity
   plateaus near `reference_velocity`'s magnitude.
5. From the resulting `(distance, u_tangential)` samples: `delta_99` by
   linear interpolation where `u/U_e` crosses `0.99` (matching
   `--verify-flat-plate`'s existing interpolation approach exactly, just
   generalized to an arbitrary sample list instead of a known column);
   `delta*` and `theta` by trapezoidal integration of `(1 - u/U_e)` and
   `(u/U_e)(1 - u/U_e)` respectively, using the wall's own prescribed
   `(wall_u, wall_v)` (tangential component) as the first sample point at
   distance 0; `H = delta*/theta`.

**Explicit, disclosed limitation (not deferred to be discovered later):**
this heuristic assumes a locally wall-normal-ish structured cell stacking
near the wall -- true for boundary-layer-clustered meshes (AirfoilMesher's
own style, per
[docs/ns-cfl-margin-and-farfield-bc-findings.md](ns-cfl-margin-and-farfield-bc-findings.md)'s
real Bravo integration case) but not guaranteed on an arbitrary
unstructured triangulation right up to a wall, where "most-aligned-neighbor"
could wander laterally along the wall instead of marching outward. This
plan deliberately does not attempt the fully general alternative (a
point-location/ray-casting profile sampler, needing a new spatial index or
brute-force point-in-cell query that doesn't exist anywhere in this
codebase today) -- consistent with this project's recurring "simple first,
document the gap" stance (see `WallDistance.h`'s own brute-force,
no-acceleration-structure precedent). Validate directly against a real
unstructured-near-wall mesh (see "Validation" below) rather than assuming
the heuristic degrades gracefully; if it doesn't, the point-location
alternative becomes a concrete follow-up, not a hypothetical one.

New free function alongside `compute_wall_traction()`:

```cpp
struct BoundaryLayerProfile {
    int wall_face_index;
    double x_mid, y_mid;
    double delta_99, displacement_thickness, momentum_thickness, shape_factor;
    int n_cells_marched;   // diagnostic: did the march hit its cap before plateauing?
};

std::vector<BoundaryLayerProfile> compute_boundary_layer_profiles(
    const UnstructuredMesh& mesh, const std::vector<int>& wall_faces,
    const std::vector<double>& u, const std::vector<double>& v,
    double u_edge, int max_cells_per_march);
```

## Face-to-node averaging

`compute_wall_traction()` and `compute_boundary_layer_profiles()` both
produce one sample **per wall face** internally, because a face is the
natural anchor for both computations: it has a single well-defined owning
cell, an outward normal, and a tangent direction, none of which a bare
mesh node has on its own (a wall node is shared by, generically, two wall
faces). But the *reported* `wall_profile_file` output should be **per wall
node**, not per face/panel, so it lines up with actual mesh-vertex
coordinates for plotting a surface distribution against arc length or
`x`-coordinate, rather than a set of face-midpoint samples offset from
them.

A small reduction step, run once after `compute_wall_traction()`/
`compute_boundary_layer_profiles()`, averages every quantity (`Cf`, `Cp`,
`y+`, `delta_99`, `delta*`, `theta`, `H`) from a node's adjacent wall
face(s) onto that node:

- An interior wall node (shared by exactly two wall faces along the same
  patch) gets the plain average of its two adjacent faces' values.
- An endpoint node of an open wall patch (exactly one adjacent wall face)
  just takes that one face's value directly -- no averaging needed or
  possible.
- A node shared by two *different* wall patches (a corner) is reported
  once per patch it belongs to (same convention as the force/moment totals
  being kept per patch), not silently merged across patches.

This mirrors `MeshReader.cpp`'s own node/edge topology (`Face::node1`/
`node2`), so identifying each wall node's adjacent wall face(s) is a direct
lookup, not a new spatial search.

## Output wiring

Follows `resolution_report_file`/`resolution_report_interval`'s exact
existing precedent
([CaseInput.h:255-256](../include/CaseInput.h),
wiring in `run_navier_stokes()` at
[main.cpp:1731](../src/main.cpp)/[main.cpp:1763](../src/main.cpp)):
`ensure_parent_directory` at startup, `resolve_output_path` rebasing under
`scratch_dir`, opened once, a row (or block of rows) appended every N steps.

| Key | Default | Meaning |
|---|---|---|
| `wall_forces_file` | *(disabled)* | CSV: one row per `(step, patch)` -- `friction_drag`, `pressure_drag`, `total_drag`, `Cd_friction`, `Cd_pressure`, `Cd_total`, `lift`, `Cl`, `moment`, `Cm` |
| `wall_forces_interval` | 1 | Write a `wall_forces_file` block every N steps |
| `wall_profile_file` | *(disabled)* | CSV: one row per wall mesh node at the time of writing (see "Face-to-node averaging") -- `x`, `y`, `patch_name`, `tau_wall`, `Cf`, `Cp`, `y_plus`, `delta_99`, `displacement_thickness`, `momentum_thickness`, `shape_factor` |
| `wall_profile_interval` | 0 (write once, at the run's natural end, like a `write_interval = 0` VTK) | Write a `wall_profile_file` snapshot every N steps if > 0 |

Both are Navier-Stokes-only keys, same table placement as
`resolution_report_file` in [docs/MANUAL.md](MANUAL.md).

## Validation plan

1. **Skin friction / drag, exact analytic case: reuse `--verify-couette`'s
   existing 1x16 sheared mesh.** Planar Couette flow has an exact,
   closed-form wall shear stress `tau_wall = mu * U / H` at both walls and
   *zero* pressure-gradient-driven force along the flow direction (uniform
   `p` by construction) -- this is a clean, already-available, zero-new-mesh
   validation target for the force/`Cf` plumbing specifically (not
   boundary-layer thickness, which has no meaning in fully-developed
   Couette flow). New `--verify-wall-forces` CLI gate: run the existing
   Couette setup, call `compute_wall_traction()` on both wall patches,
   compare `tau_wall`/`Cf` and the resulting friction drag against the exact
   value to machine precision (same bar as `--verify-couette`'s own
   `~1e-15` result).
2. **`Cp`/`Cm`/`y+`, exact/closed-form checks on the same Couette setup.**
   Couette flow has uniform pressure by construction, so `Cp` must come out
   exactly `0` everywhere -- a good degenerate check for the `Cp` plumbing
   specifically (its own force contribution is zero, so this doesn't
   overlap with the `tau_wall`/`Cd` check above). `y+` is exactly computable
   by hand on this mesh (`tau_wall`, `rho`, `mu`, and the first-cell
   wall-normal distance are all known/uniform), so compare against that.
   `Cm` also has a closed form here: both walls carry equal-and-opposite
   pure shear forces offset by `+-H/2` from the domain's vertical center, so
   taking `(moment_reference_x, moment_reference_y)` at that center gives an
   exact, nonzero analytic moment to check against -- a real check, not a
   degenerate zero, since a `Cm` implementation with a sign error would
   still pass a same-value-both-sides check but fail this one.
3. **Boundary-layer thickness, approximate analytic case: extend
   `--verify-flat-plate`.** That harness already computes `delta_99`,
   `tau_wall`, `u_tau` by hand for its one hardcoded structured mesh
   ([main.cpp](../src/main.cpp)) and compares the *velocity profile shape*
   against the Pohlhausen approximation. Once the general
   `compute_boundary_layer_profiles()` exists, rerun that same case through
   it and confirm the general code reproduces today's hardcoded `delta_99`
   result (regression check), then add a Blasius laminar skin-friction
   comparison (`Cf(x) = 0.664/sqrt(Re_x)`) at the exit station, and Pohlhausen's
   own closed-form `delta*`/`theta` ratios as an additional check beyond the
   velocity-profile-only comparison that exists today.
4. **Structured-assumption stress test.** Build (or reuse, if one already
   exists) a mesh with an unstructured triangulation near a wall and run
   `compute_boundary_layer_profiles()` against it specifically to observe
   and document the marching heuristic's actual failure mode, rather than
   asserting one without evidence.
5. **Regression.** Rerun every existing `--verify-*` gate to confirm the new
   post-processing pass (called only when its output file is configured, and
   never from inside `step()`) changes nothing about existing solver
   behavior.

## Phasing

1. **Phase 1 -- forces/Cf/Cp/Cm/y+**: `WallTraction.h`/`.cpp`,
   reference-quantity case-file keys (including `moment_reference_x/y`),
   `wall_forces_file`/`wall_profile_file` wiring, face-to-node averaging,
   `--verify-wall-forces` against Couette's exact `tau_wall`/`Cp`/`y+`/`Cm`.
2. **Phase 2 -- boundary-layer thickness**: cell-column marching,
   `wall_profile_file` wiring, `--verify-flat-plate` extended with the
   Blasius/Pohlhausen `Cf`/`delta*`/`theta` checks, the unstructured-mesh
   stress test.
3. **Phase 3 -- deferred/optional**: wire into `RANSFVMSolver` if/when it
   ever gets a `CaseInput`/`main.cpp` case-file path; revisit the
   point-location/ray-casting alternative only if Phase 2's stress test
   shows the marching heuristic actually fails on a realistic mesh.

## Phase 2 stress-test finding (implemented 2026-07-05)

Phases 1 and 2 are both implemented (`include/WallTraction.h` /
`src/WallTraction.cpp`, case-file keys, `wall_forces_file`/`wall_profile_file`
wiring, `--verify-wall-forces`, `--verify-flat-plate`'s Blasius/Pohlhausen
extension). The Validation step 4 stress test
(`--verify-bl-marching-unstructured`, `run_verify_bl_marching_unstructured()`
in `src/main.cpp`) **fails as designed** on a checkerboard-triangulated
near-wall mesh: a manufactured exact planar-Couette field (`u = U*y/H`, exact
`delta_99 = 0.99*H` regardless of triangulation) comes back with a
systematic ~20-25% `delta_99` error, essentially identical across every wall
face despite their differing local diagonal orientation.

**Root cause, not just a symptom**: the failure is NOT lateral wandering --
the "most-aligned-neighbor" cell selection itself consistently picks the
cell stacked above, even across the diagonal cuts (fixed by requiring a
strictly positive alignment before accepting a candidate face, without which
the march would erroneously wander sideways once it reaches the top of a
column and its only remaining interior neighbor is lateral -- a real bug
found and fixed during this same implementation pass). The remaining error
is a **distance-accounting bias**: `compute_boundary_layer_profiles()`
accumulates `face_normal_distance()` at each step, which projects onto THAT
CROSSED FACE's own normal -- correct when every crossed face's normal
already equals the march's fixed global direction (true for a structured
quad stack, where this exactly reproduces `--verify-flat-plate`'s hardcoded
column-walk result to under 1%), but wrong when a step crosses a diagonal
face whose normal is tilted away from the true wall-normal direction, where
it silently sums a distance measured along the wrong axis instead.

**Conclusion**: the "structured cell stacking" assumption this document
already disclosed is real, not hypothetical, and manifests specifically as
a `delta_99`/thickness accounting error (not a crashed or wildly wrong
march) on a genuinely unstructured near-wall triangulation. The
point-location/ray-casting alternative floated above is now a concrete
follow-up for anyone routing this module's boundary-layer-thickness output
through a mesh with unstructured near-wall triangulation (e.g. a Delaunay
mesher) -- not needed for `Cf`/`Cp`/`Cd`/`Cl`/`Cm`/`y+`, which don't use
`compute_boundary_layer_profiles()` at all and are validated to near-machine
precision by `--verify-wall-forces` on a non-orthogonal (sheared) mesh.

## Phase 3, part 1: point-location alternative (implemented 2026-07-05)

The point-location alternative motivated by the finding above is now built:
`compute_boundary_layer_profiles_point_location()`
([WallTraction.h](../include/WallTraction.h)/[.cpp](../src/WallTraction.cpp)).
Rather than choosing a neighbor cell to cross (the source of marching's
failure mode), it walks straight out from each wall face midpoint along its
fixed inward normal at `n_samples` evenly-spaced offsets up to
`max_distance`, brute-force point-locating (point-in-polygon over every
cell, no spatial acceleration structure -- same "simple first" precedent as
`WallDistance.h`'s brute-force point-to-segment search) which cell contains
each sample point and reading that cell's velocity there. This has no
near-wall-topology dependence at all, at the cost of an `O(cells)` search
per sample point (meaningfully more expensive than marching's `O(1)`-ish per
step on a large mesh).

The delta_99/displacement/momentum-thickness derivation itself (the actual
math, not the sampling) was refactored out of the marching function into a
shared `derive_boundary_layer_profile()` helper, so both methods produce a
`BoundaryLayerProfile` via identical downstream logic -- they differ only in
how the raw `(distance, u_tangential)` samples are obtained.

**Validated on the exact mesh that broke marching**: `--verify-bl-point-location`
reuses `run_verify_bl_marching_unstructured()`'s identical checkerboard-
triangulated mesh and manufactured Couette-like field, and reproduces the
exact `delta_99` to within 0.76% (vs. marching's ~20-25% error there) --
confirming this alternative doesn't share marching's failure mode, at a bar
5x tighter than what would even count as passing.

**Wiring**: a new `boundary_layer_method` case-file key (`marching` default
| `point-location`) selects which method `wall_profile_file` uses, plus
`boundary_layer_max_distance` (default: automatic, the mesh's own
bounding-box diagonal -- an over-generous distance only costs a few extra
out-of-domain samples, since point-location stops as soon as a sample point
falls outside the mesh) and `boundary_layer_n_samples` (default 200). A
matching `NavierStokesFVMSolver::compute_boundary_layer_profile_samples_point_location()`
wrapper mirrors the existing marching wrapper. See
[docs/MANUAL.md](MANUAL.md)'s Wall diagnostics section for the case-file
key table.

**Not done**: Phase 3's other half, wiring `RANSFVMSolver` into
`CaseInput`/`main.cpp`, remains untouched -- a separate, larger,
pre-existing gap (see
[docs/archive/rans-spalart-allmaras-tracker.md](archive/rans-spalart-allmaras-tracker.md))
outside this plan's own scope.
