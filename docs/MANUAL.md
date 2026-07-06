# FiniteVolume manual

A standalone, dependency-free 2D unstructured finite-volume solver. It solves
scalar diffusion, advection-diffusion of a passive scalar, the compressible
Euler equations, the compressible (laminar) Navier-Stokes equations, or the
Navier-Stokes equations closed with a Spalart-Allmaras RANS turbulence model,
on an arbitrary polygon mesh, driven by a plain-text case file, and writes
results as legacy VTK files for viewing in ParaView (or similar). See
[RANS (Spalart-Allmaras) turbulence closure](#rans-spalart-allmaras-turbulence-closure)
for what the RANS closure does and does not validate before relying on it.

This document is meant to be read start to finish by a human, and to serve
as a lookup reference (by section) for an automated coding agent working on
this codebase. Every claim below is grounded in a specific file; when in
doubt, the cited source file is the ground truth, not this document.

## Contents

1. [Overview & architecture](#overview--architecture)
2. [Building](#building)
3. [Running](#running)
4. [Case file reference](#case-file-reference)
5. [Mesh formats](#mesh-formats)
6. [Solvers](#solvers)
7. [Gradient reconstruction (shared infrastructure)](#gradient-reconstruction-shared-infrastructure)
8. [Boundary conditions](#boundary-conditions)
9. [Output format](#output-format)
10. [Monitoring & stopping criteria](#monitoring--stopping-criteria)
11. [Checkpointing & restart](#checkpointing--restart)
12. [Performance (OpenMP)](#performance-openmp)
13. [Known limitations](#known-limitations)
14. [File map](#file-map)

---

## Overview & architecture

The solver is a single C++17 executable (`FiniteVolume.exe`, built via
[../CMakeLists.txt](../CMakeLists.txt)) with no external library
dependencies. The pipeline for a normal run is:

```
case file (.case/.txt, key=value)  --CaseInput::load-->  run parameters
                                                                |
mesh file (.msh or .fvmesh)  --MeshReader::read-->  UnstructuredMesh
                                                                |
                                     boundary conditions matched by patch name
                                                                |
  UnstructuredFVMSolver / AdvectionDiffusionFVMSolver / EulerFVMSolver / NavierStokesFVMSolver / RANSFVMSolver
                                                                |
                                    step() loop (residual/checkpoint/snapshots)
                                                                |
                                        VtkWriter  -->  .vtk result file(s)
```

Everything is 2D: `Node` is `(x, y)`, a `Face` is a line segment, a `Cell` is
a closed polygon of faces. There is no 3D support and none is planned to
this data model — this also means **there is no true 3D DNS/turbulence
capability here, ever**: classical Kolmogorov scaling and the turbulent
energy cascade are 3D results (vortex stretching, the actual cascade
mechanism, doesn't exist in 2D the same way). See
[Known limitations](#known-limitations) and
[docs/navier-stokes-tracker.md](navier-stokes-tracker.md) Phase 5.

Core data structure: [`UnstructuredMesh`](../include/UnstructuredMesh.h)
(`Node`, `Face`, `Cell`, `BoundaryPatch`). A `Cell` stores its faces and
vertex indices as `std::vector<int>`, so **a cell may have any number of
faces/vertices (3, 4, 20, or more) — there is no fixed-size limit anywhere
in the mesh, geometry, or flux-assembly code.** See
[../docs/fvmesh-format.md](fvmesh-format.md) for why this matters and how to
author such a mesh directly.

Five independent solvers are wired into the case-file-driven run path and
share this one mesh representation:
[`UnstructuredFVMSolver`](../include/UnstructuredFVMSolver.h) (scalar
diffusion), [`AdvectionDiffusionFVMSolver`](../include/AdvectionDiffusionFVMSolver.h)
(advection-diffusion of a passive scalar), [`EulerFVMSolver`](../include/EulerFVMSolver.h)
(compressible Euler), [`NavierStokesFVMSolver`](../include/NavierStokesFVMSolver.h)
(compressible Navier-Stokes — Euler plus a Newtonian viscous stress tensor
and Fourier heat conduction), and [`RANSFVMSolver`](../include/RANSFVMSolver.h)
(Navier-Stokes plus a Spalart-Allmaras one-equation RANS turbulence closure —
see [RANS (Spalart-Allmaras) turbulence closure](#rans-spalart-allmaras-turbulence-closure)
for what this closure does and does not validate before relying on it).
[`src/main.cpp`](../src/main.cpp) picks one based on the case file's
`equation` key and drives its step loop directly (not inside the solver
class) so it can layer on residual tracking, periodic snapshots, stopping
criteria, and checkpointing uniformly for all five.

The latter two solvers (advection-diffusion's diffusive term, and
Navier-Stokes' viscous stress/heat conduction) share a second piece of
infrastructure: [`GradientCalculator`](../include/GradientReconstruction.h),
which reconstructs cell gradients (Green-Gauss or Least-Squares) and a
non-orthogonality-corrected face-normal derivative, needed because a naive
two-point difference is only exact on an orthogonal mesh — see
[Gradient reconstruction](#gradient-reconstruction-shared-infrastructure).

## Building

Requires CMake >= 3.20 and a C++17 compiler with OpenMP support (the
project's own toolchain is MinGW-w64 g++, built via Ninja on Windows; the
`-static` link means the resulting `.exe` has no runtime DLL dependencies).

```
cmake --build build --config Release   # or Debug
```

The single target `FiniteVolume` is produced at `build/bin/FiniteVolume.exe`
regardless of config (see the `CMAKE_RUNTIME_OUTPUT_DIRECTORY*` overrides in
[../CMakeLists.txt](../CMakeLists.txt)). Source files are picked up by a
`file(GLOB CONFIGURE_DEPENDS src/*.cpp)`, so adding a new `.cpp` under
`src/` is picked up automatically on the next build (CMake reconfigures
itself when the glob's file list changes) — no need to hand-edit
`CMakeLists.txt` for new solver/reader/writer source files.

## Running

```
FiniteVolume.exe <case_file>                     # run a simulation
FiniteVolume.exe --version                       # print build/version info
FiniteVolume.exe --validate-mesh <mesh_file>     # parse a mesh, print stats, exit
FiniteVolume.exe --verify-gradient               # self-contained GradientCalculator check, exit
FiniteVolume.exe --verify-advdiff                # self-contained advection-diffusion check, exit
FiniteVolume.exe --verify-ns-uniform             # self-contained NS uniform-flow check, exit
FiniteVolume.exe --verify-couette                # self-contained Couette flow validation, exit
FiniteVolume.exe --verify-ns-stretched-cfl       # self-contained anisotropic-mesh compute_dt() check, exit
FiniteVolume.exe --verify-wall-distance          # self-contained wall-distance module check, exit
FiniteVolume.exe --verify-sa-source              # self-contained SA source-term check, exit
FiniteVolume.exe --verify-rans-stability         # self-contained RANSFVMSolver stability check, exit
FiniteVolume.exe --verify-flat-plate             # self-contained RANS flat-plate boundary-layer check, exit
FiniteVolume.exe --verify-wall-forces            # self-contained wall force/Cf/Cp/y+/Cm check, exit
FiniteVolume.exe --verify-bl-marching-unstructured  # boundary-layer-marching unstructured-mesh stress test, exit
FiniteVolume.exe --verify-bl-point-location      # point-location BL alternative on the same stress-test mesh, exit
```

- `<case_file>` is a plain-text key=value file (see
  [Case file reference](#case-file-reference)); its `mesh_file` key points
  at a `.msh` or `.fvmesh` file, and its `output_file` key is where the
  result VTK is written.
- `--validate-mesh` parses the given mesh (dispatching on extension, same as
  a real run would) and prints node/cell/face/patch counts and the
  min/max cell volume, without running any solver — useful to sanity-check
  a hand-written or generated mesh before spending time on a full case.
- `--verify-gradient`/`--verify-advdiff`/`--verify-ns-uniform`/`--verify-couette`/
  `--verify-ns-stretched-cfl` each build a small, deterministic, self-contained
  test mesh in memory (no external file needed), run one specific correctness
  check against a known analytic result (or, for `--verify-ns-stretched-cfl`,
  a real anisotropic-mesh divergence found via Bravo integration testing),
  print a PASS/FAIL verdict with the measured error, and exit — these are this
  project's substitute for a unit-test suite (there is no `ctest`/external
  test framework; see [Building](#building)). See
  [main.cpp](../src/main.cpp)'s `run_verify_gradient()`/
  `run_verify_advection_diffusion()`/`run_verify_navier_stokes_uniform()`/
  `run_verify_couette()`/`run_verify_ns_stretched_cfl()`, and
  [docs/navier-stokes-tracker.md](navier-stokes-tracker.md) for what each one
  actually verifies and why (including two dead ends found and corrected
  while building the Couette one — worth reading before trusting a
  superficially similar test mesh in a new context).
- `--verify-wall-distance`/`--verify-sa-source`/`--verify-rans-stability`/
  `--verify-flat-plate` are the equivalent gates for the RANS
  (Spalart-Allmaras) capability — see
  [RANS (Spalart-Allmaras) turbulence closure](#rans-spalart-allmaras-turbulence-closure)
  and [docs/archive/rans-spalart-allmaras-tracker.md](archive/rans-spalart-allmaras-tracker.md)
  for what each one verifies.
- `--verify-wall-forces` checks [`WallTraction.h`](../include/WallTraction.h)'s
  `tau_wall`/`Cf`/`Cp`/`y+`/`Cm` against closed-form planar Couette flow
  values; `--verify-bl-marching-unstructured` stress-tests
  `compute_boundary_layer_profiles()`'s cell-marching heuristic against a
  genuinely unstructured (checkerboard-triangulated) near-wall mesh and
  **fails as designed**; `--verify-bl-point-location` reruns that same
  failing mesh through `compute_boundary_layer_profiles_point_location()`
  and confirms the point-location alternative doesn't share the failure —
  see [Wall diagnostics](#wall-diagnostics) and
  [docs/wall-diagnostics-plan.md](wall-diagnostics-plan.md) for what that
  failure means and why it's expected.
- Ctrl+C during a run requests a clean stop after the in-flight step
  completes (same code path as a natural convergence/step-limit stop:
  checkpoint is written if configured, final output is still written). See
  `g_interrupt_requested` in [../src/main.cpp](../src/main.cpp).

## Case file reference

Parsed by [`CaseInput::load`](../src/CaseInput.cpp). Format: one entry per
line, `key = value` or `boundary <patch> <type> <value...>`; `#` starts a
comment (to end of line); blank lines are ignored. Unrecognized `boundary`
types or an unrecognized `equation`/`euler_init` value are hard errors
(`load` returns `false`); unrecognized `key = value` keys are silently
ignored (no typo detection).

| Key | Default | Meaning |
|---|---|---|
| `mesh_file` | *(required)* | Path to a `.msh` or `.fvmesh` mesh file |
| `output_file` | *(required)* | Path to write the result `.vtk` file |
| `equation` | `diffusion` | `diffusion`, `advection_diffusion`, `euler`, `navier_stokes`, or `rans` — selects the physics |
| `nsteps` | 500 | Absolute target step count (also the resume target — see [Checkpointing](#checkpointing--restart)) |
| `residual_file` | *(disabled)* | CSV path for per-step residual history |
| `residual_interval` | 1 | Write a `residual_file` row every N steps |
| `write_interval` | 0 (disabled) | Write a numbered VTK snapshot every N steps |
| `num_threads` | 0 (OpenMP default) | OpenMP thread count |
| `checkpoint_file` | *(disabled)* | Path to save/auto-resume solver state |
| `residual_tolerance` | disabled | **Diffusion / advection-diffusion only.** Stop once residual < this |
| `residual_tolerance_rho`, `_rho_u`, `_rho_v`, `_E` | disabled | **Euler / Navier-Stokes / RANS only.** Each independently optional; converged requires every one the user set to be satisfied simultaneously |
| `residual_tolerance_nut` | disabled | **RANS only.** Checked alongside the four above |
| `scratch_dir` | *(disabled)* | Base directory for *relative* `output_file`/`checkpoint_file`/`residual_file` paths; absolute paths and `mesh_file` are unaffected |
| `output_precision` | 6 | Significant digits for every output double (VTK results/snapshots, `residual_file` CSV, `FvMeshWriter`); valid range 1-17 |
| `gradient_scheme` | `least-squares` | **Advection-diffusion / Navier-Stokes / RANS only.** `least-squares` or `green-gauss` — see [Gradient reconstruction](#gradient-reconstruction-shared-infrastructure) |

**Diffusion-only** (`equation = diffusion`, the default):

| Key | Default | Meaning |
|---|---|---|
| `alpha` | 0.01 | Diffusion coefficient |
| `dt` | 0.001 | Fixed time step (not adaptive — see [Solvers](#solvers)) |
| `initial_value` | 1.0 | phi inside the initial disc |
| `initial_radius` | 0.25 | Radius (from the origin) of the initial disc |
| `boundary <patch> dirichlet <value>` | — | Pin phi at the face to `<value>` |
| `boundary <patch> neumann <value>` | — | Prescribe flux `<value>` = `-alpha * dPhi/dn` directly |

**Advection-diffusion-only** (`equation = advection_diffusion`): same
`alpha`/`dt`/`initial_value`/`initial_radius`/`boundary dirichlet`/
`boundary neumann` keys as diffusion above (a scalar Dirichlet/Neumann BC
means the same thing for both), plus:

| Key | Default | Meaning |
|---|---|---|
| `u_adv`, `v_adv` | 0.0, 0.0 | Uniform, prescribed (not solved-for) advection velocity components |

**Euler-only** (`equation = euler`):

| Key | Default | Meaning |
|---|---|---|
| `gamma` | 1.4 | Ratio of specific heats |
| `cfl` | 0.5 | CFL number controlling the adaptive time step |
| `flux_scheme` | `rusanov` | Numerical flux at every face: `rusanov`, `hllc`, or `exact` (see [Discontinuity handling](#discontinuity-handling-shock-capturing)) |
| `exact_riemann_tol` | `1e-6` | Relative pressure tolerance for the `exact` flux scheme's Newton-Raphson star-pressure solve. Ignored unless `flux_scheme = exact`. **Not meant to be tuned in normal use** — exposed for transparency only; the default is the value this solver always used before this key existed. |
| `exact_riemann_max_iter` | `20` | Iteration cap for that same Newton-Raphson solve (stops the iteration even if `exact_riemann_tol` is never reached). Ignored unless `flux_scheme = exact`. **Not meant to be tuned in normal use**, for the same reason as `exact_riemann_tol` above. |
| `euler_init = freestream <rho> <u> <v> <p>` | — | Uniform initial state everywhere |
| `euler_init = tworegion <rho_l> <u_l> <v_l> <p_l> <rho_r> <u_r> <v_r> <p_r> <x0>` | — | Left/right split at `x_centroid = x0` (e.g. a Sod shock tube) |
| `boundary <patch> wall` | — | Slip / no-penetration |
| `boundary <patch> farfield <rho> <u> <v> <p>` | — | Fixed prescribed state |
| `boundary <patch> outflow` | — | Zero-gradient extrapolation |

**Navier-Stokes-only** (`equation = navier_stokes`): same `gamma`/`cfl`/
`flux_scheme`/`exact_riemann_tol`/`exact_riemann_max_iter` keys as Euler
above (the inviscid part reuses `EulerFVMSolver`'s exact flux machinery
unchanged), plus:

| Key | Default | Meaning |
|---|---|---|
| `mu` | 0.0 | Dynamic viscosity; `0` = no viscous stress (and no heat conduction, since conductivity scales with `mu`) |
| `prandtl` | 0.72 | Prandtl number (air's value). **Unvalidated** — a zero/negative value silently produces `Inf`/division-by-zero |
| `gas_constant` | 1.0 | Specific gas constant `R` in `p = rho*R*T`; only used for temperature/heat conduction, never pressure/sound speed (Euler-inherited physics never needed `R`). **Unvalidated**, same caveat as `prandtl` |
| `resolution_report_file` | *(disabled)* | CSV path for a resolution-adequacy diagnostic — **not** a claim about 3D DNS; see [Known limitations](#known-limitations) |
| `resolution_report_interval` | 1 | Write a `resolution_report_file` row every N steps |
| `wall_forces_file` | *(disabled)* | CSV path, one row per `(step, patch)`: `friction_drag`/`pressure_drag`/`total_drag`/`cd_friction`/`cd_pressure`/`cd_total`/`lift`/`cl`/`moment`/`cm` over every `ns_wall*`/`rans_wall*` patch, plus a `TOTAL` row — see [Wall diagnostics](#wall-diagnostics) |
| `wall_forces_interval` | 1 | Write a `wall_forces_file` row block every N steps |
| `wall_profile_file` | *(disabled)* | CSV path, one row per wall mesh node: `x`/`y`/`patch_name`/`tau_wall`/`Cf`/`Cp`/`y_plus`/`delta_99`/`displacement_thickness`/`momentum_thickness`/`shape_factor` — see [Wall diagnostics](#wall-diagnostics) |
| `wall_profile_interval` | 0 | Write a numbered `wall_profile_file` snapshot every N steps if `> 0`; the plain `wall_profile_file` path itself is always (re)written once at the run's natural end regardless |
| `reference_density`, `reference_velocity_x`, `reference_velocity_y`, `reference_pressure` | *(auto)* | Freestream values for `Cf`/`Cp`/`Cd`/`Cl`/`Cm`'s dynamic-pressure normalization and drag direction; if any is unset, defaults from the first `ns_farfield` (or, for `equation = rans`, `rans_farfield`) patch's prescribed state. Only resolved/validated (and a load-time error if no such patch exists to default from) when `wall_forces_file` or `wall_profile_file` is set |
| `reference_length` | 1.0 | Length scale for `Cd`/`Cl`/`Cm` (e.g. chord); `1.0` degrades gracefully to a per-unit-span coefficient |
| `moment_reference_x`, `moment_reference_y` | 0.0, 0.0 | Point `Cm` is taken about (e.g. an airfoil's quarter-chord) |
| `boundary_layer_method` | `marching` | `marching` (cheap, assumes a wall-normal cell stacking) or `point-location` (more expensive, no near-wall-topology assumption) — see [Wall diagnostics](#wall-diagnostics) |
| `boundary_layer_max_distance` | *(auto)* | Point-location only; furthest wall-normal distance to sample. `<= 0` (or omitted) auto-defaults to the mesh's own bounding-box diagonal |
| `boundary_layer_n_samples` | 200 | Point-location only; number of evenly-spaced sample points |
| `ns_init = freestream <rho> <u> <v> <p>` | — | Uniform initial state everywhere (own storage; same grammar as `euler_init`) |
| `ns_init = tworegion <rho_l> <u_l> <v_l> <p_l> <rho_r> <u_r> <v_r> <p_r> <x0>` | — | Left/right split at `x_centroid = x0` |
| `boundary <patch> ns_wall` | — | No-slip, stationary, adiabatic (zero heat flux) |
| `boundary <patch> ns_wall_isothermal <T>` | — | No-slip, stationary, fixed wall temperature `T` |
| `boundary <patch> ns_wall_moving <u> <v>` | — | No-slip, moving at `(u, v)`, adiabatic (e.g. Couette flow) |
| `boundary <patch> ns_wall_moving_isothermal <u> <v> <T>` | — | No-slip, moving, fixed wall temperature |
| `boundary <patch> ns_farfield <rho> <u> <v> <p>` | — | Fixed prescribed state (inviscid contribution only) |
| `boundary <patch> ns_outflow` | — | Zero-gradient extrapolation (inviscid contribution only) — **not** a periodic BC; see [Known limitations](#known-limitations) before relying on it to enforce translational invariance |

`ns_*` boundary keywords are prefixed specifically because `farfield`/
`outflow` already mean something for Euler and must land in a different
internal spec vector — see [`CaseInput.cpp`](../src/CaseInput.cpp).

**RANS-only** (`equation = rans`): same `gamma`/`cfl`/`flux_scheme`/
`exact_riemann_tol`/`exact_riemann_max_iter`/`mu`/`prandtl`/`gas_constant`/
`gradient_scheme` keys as Navier-Stokes above (`mu`/`prandtl`/`gas_constant`
govern the mean-flow molecular terms exactly as they do for Navier-Stokes),
plus:

| Key | Default | Meaning |
|---|---|---|
| `prandtl_t` | 0.9 | Turbulent Prandtl number (air's value); no Navier-Stokes analog |
| `initial_nut` | 0.0 | Uniform initial nu-tilde applied to every cell. SA is documented as sensitive to this choice — see [RANS (Spalart-Allmaras) turbulence closure](#rans-spalart-allmaras-turbulence-closure) |
| `sa_cb1`, `sa_cb2`, `sa_sigma`, `sa_kappa`, `sa_cw2`, `sa_cw3`, `sa_cv1`, `sa_cv2`, `sa_cv3` | the standard SA-noft2 set (see [SpalartAllmaras.h](../include/SpalartAllmaras.h)) | Each independently optional. **Not meant to be tuned in normal use** — exposed for experimentation, same spirit as `exact_riemann_tol` above |
| `rans_init = freestream <rho> <u> <v> <p>` | — | Uniform initial mean-flow state everywhere (own storage; same grammar as `ns_init`) |
| `rans_init = tworegion <rho_l> <u_l> <v_l> <p_l> <rho_r> <u_r> <v_r> <p_r> <x0>` | — | Left/right split at `x_centroid = x0` |
| `boundary <patch> rans_wall` | — | No-slip, stationary, adiabatic; wall nut is fixed at 0 by the SA model itself, not case-configurable |
| `boundary <patch> rans_wall_isothermal <T>` | — | No-slip, stationary, fixed wall temperature `T` |
| `boundary <patch> rans_wall_moving <u> <v>` | — | No-slip, moving at `(u, v)`, adiabatic |
| `boundary <patch> rans_wall_moving_isothermal <u> <v> <T>` | — | No-slip, moving, fixed wall temperature |
| `boundary <patch> rans_farfield <rho> <u> <v> <p> <nut>` | — | Fixed prescribed mean-flow state plus a prescribed freestream nut (inviscid contribution only) |
| `boundary <patch> rans_outflow` | — | Zero-gradient extrapolation, mean flow and nut alike (inviscid contribution only) |

`wall_forces_file`/`wall_profile_file`/`reference_*`/`boundary_layer_*` (see
Navier-Stokes above) all apply to RANS too, using `mu + rho*nu_t` as the
effective viscosity in place of `mu` — see
[Wall diagnostics](#wall-diagnostics).

A patch present in the mesh but with no matching `boundary` line defaults
to Dirichlet 0.0 (diffusion/advection-diffusion), Wall (Euler), `ns_wall` —
adiabatic, stationary (Navier-Stokes), or `rans_wall` — adiabatic, stationary
(RANS), and prints a warning — it does not fail the run.

Example (diffusion):
```
mesh_file = mesh/square.msh
output_file = output/result.vtk
alpha = 0.01
dt = 0.001
nsteps = 500
initial_value = 1.0
initial_radius = 0.25

boundary wall dirichlet 0.0
boundary inlet dirichlet 1.0
```

Example (Euler, Sod shock tube):
```
mesh_file = mesh/tube.msh
output_file = output/result.vtk
equation = euler
gamma = 1.4
cfl = 0.5
nsteps = 500

euler_init = tworegion 1.0 0.0 0.0 1.0  0.125 0.0 0.0 0.1  0.0

boundary left wall
boundary right wall
```

Example (advection-diffusion of a passive scalar):
```
mesh_file = mesh/square.msh
output_file = output/result.vtk
equation = advection_diffusion
alpha = 0.1
u_adv = 1.0
v_adv = 0.0
dt = 0.0005
nsteps = 40000

boundary left dirichlet 0.0
boundary right dirichlet 1.0
boundary bottom neumann 0.0
boundary top neumann 0.0
```

Example (Navier-Stokes, planar Couette flow — a stationary bottom wall, a
top wall moving at `U = 0.1`, no imposed pressure gradient; the exact steady
solution is `u(y) = U*y/H`, which this exact case is validated against in
`run_verify_couette()`):
```
mesh_file = mesh/channel.fvmesh
output_file = output/result.vtk
equation = navier_stokes
gamma = 1.4
cfl = 0.3
mu = 0.02
prandtl = 0.72
gas_constant = 1.0
nsteps = 20000

ns_init = freestream 1.0 0.0 0.0 1.0

boundary bottom ns_wall
boundary top ns_wall_moving 0.1 0.0
boundary left ns_outflow
boundary right ns_outflow

# Wall diagnostics (see Wall diagnostics): no ns_farfield patch exists here,
# so reference_* must be given explicitly or CaseInput::load() fails loudly.
wall_forces_file = output/wall_forces.csv
wall_profile_file = output/wall_profile.csv
reference_density = 1.0
reference_velocity_x = 0.1
reference_velocity_y = 0.0
reference_pressure = 1.0
```

Example (RANS, planar Couette flow with a freestream inlet/outlet added on
top of the Navier-Stokes example above — a physically odd combination kept
deliberately simple here just to exercise every RANS-only key; see
[RANS (Spalart-Allmaras) turbulence closure](#rans-spalart-allmaras-turbulence-closure)
before trusting this closure's turbulent behavior on a real case):
```
mesh_file = mesh/channel.fvmesh
output_file = output/result.vtk
equation = rans
gamma = 1.4
cfl = 0.3
mu = 0.02
prandtl = 0.72
prandtl_t = 0.9
gas_constant = 1.0
initial_nut = 0.06
nsteps = 20000

rans_init = freestream 1.0 0.1 0.0 1.0

boundary bottom rans_wall
boundary top rans_wall_moving 0.1 0.0
boundary left rans_farfield 1.0 0.1 0.0 1.0 0.06
boundary right rans_outflow

residual_file = output/residual.csv
checkpoint_file = output/checkpoint.bin
wall_forces_file = output/wall_forces.csv
wall_profile_file = output/wall_profile.csv
reference_length = 1.0
```

## Mesh formats

Dispatched by extension in [`MeshReader::read`](../src/MeshReader.cpp)
(`.msh` -> Gmsh reader, `.fvmesh` -> native reader). Both populate the same
[`UnstructuredMesh`](../include/UnstructuredMesh.h) and both build cell
adjacency/faces via the same shared, format-agnostic
`build_cells_faces_patches()` helper — so both formats produce byte-for-byte
equivalent geometry for the same logical mesh.

### Gmsh ASCII (`.msh`)

[`MeshReader::read_gmsh`](../src/MeshReader.cpp), supporting legacy format
2.2 and format 4.x (auto-detected from `$MeshFormat`; binary Gmsh files are
**not** supported). **Only Gmsh element types 1 (line), 2 (triangle), and 3
(quadrangle) are recognized** — this is a limitation of the *parser*, not
the solver; the underlying mesh/solver code has no cell-type restriction at
all (see below).

### FVMESH (`.fvmesh`) — this project's own format

[`MeshReader::read_fvmesh`](../src/MeshReader.cpp) /
[`FvMeshWriter::write`](../src/FvMeshWriter.cpp). Full grammar:
[docs/fvmesh-format.md](fvmesh-format.md). In short: `POINTS`, `CELLS`
(each cell is `n_verts v0 v1 ... v_{n_verts-1}`, **any `n_verts >= 3`**),
`BOUNDARY` (`node_a node_b patch_name`, patches created on demand by name —
no Gmsh-style numeric physical-group indirection). This is the format to
target if you're writing your own mesh generator, and the only format that
can express a cell with more faces than a quadrilateral.

### Arbitrary-polygon support

Neither `Cell::faces` nor `Cell::node_ids` are fixed-size — both are
`std::vector<int>` ([UnstructuredMesh.h](../include/UnstructuredMesh.h)) —
and every flux-assembly loop in all four solvers iterates
`for (int face_idx : mesh.cells[c].faces)`, so a 20-sided (or larger)
polygon cell works with zero solver changes once it's loaded. Cell
volume/centroid use the shoelace formula over however many vertices the
cell has ([MeshReader.cpp](../src/MeshReader.cpp)). The `.fvmesh` format is
how you actually get such a cell into the solver today (Gmsh's own
triangle/quad-only element types can't express one).

### No mesh-quality safeguards

**Neither reader nor the geometry code validates mesh quality.** There is no
check anywhere for: near-zero-length faces, near-zero-volume cells,
near-coincident/near-parallel adjacent face normals, or any other skewness/
degeneracy metric. `compute_geometry()`
([UnstructuredMesh.h](../include/UnstructuredMesh.h)) divides by
`face.area` unconditionally; a degenerate mesh will silently produce
`NaN`/`Inf` face normals rather than a clear error. `--validate-mesh`
reports counts and volume range, which is enough to catch a grossly wrong
mesh, but it is not a mesh-quality checker. If you're generating meshes
programmatically, validate quality upstream, not here.

## Solvers

All four solvers advance an explicit, first-order-in-time, cell-centered
finite volume scheme: integrate flux over each cell's bounding faces
(divergence theorem), then `state_new = state_old + dt/volume *
net_flux_in`. All four compute face fluxes and per-cell residuals as two
separate OpenMP-parallel passes (per-face flux, then per-cell scatter-gather)
specifically to avoid a race on cells shared by two faces. See
[docs/dns-higher-order-scheme-plan.md](dns-higher-order-scheme-plan.md) for
a deferred discussion of what a higher-order-in-time (RK3/RK4) and
higher-order-in-space (MUSCL-reconstructed face states) alternative would
require, and why "first-order-in-time, first-order-in-space" is a real
accuracy ceiling for any of these solvers, not just a detail.

### Scalar diffusion (`equation = diffusion`)

[`UnstructuredFVMSolver`](../src/UnstructuredFVMSolver.cpp) solves
`dPhi/dt = alpha * Laplacian(Phi)`. Internal faces use a two-point central
difference between neighboring cell centroids for the normal gradient;
Dirichlet boundaries use a one-point difference between the face midpoint
and the owning cell's centroid; Neumann boundaries apply the prescribed flux
directly. **`dt` is fixed** (the case file's `dt` key) — the classic
explicit-diffusion stability limit (`dt` small relative to `alpha` and cell
size) is **not checked**; an unstable choice will simply diverge (residual
goes NaN/Inf), which the run loop in [main.cpp](../src/main.cpp) detects and
stops on.

There is no shock/discontinuity handling here because the diffusion equation
is parabolic (smoothing), not hyperbolic — nothing to capture.

### Compressible Euler (`equation = euler`)

[`EulerFVMSolver`](../src/EulerFVMSolver.cpp) solves the 2D compressible
Euler equations in conserved form
(`EulerState` = `rho, rho_u, rho_v, E`, see
[EulerState.h](../include/EulerState.h)). Unlike diffusion, **`dt` is
adaptive**: recomputed every step from the CFL condition
(`EulerFVMSolver::compute_dt`), taking the most restrictive `cfl * length /
Smax` over every face, where `Smax` is the fastest wave speed (`|Vn| + speed
of sound`) touching that face. This is necessary because the stable step for
a hyperbolic system shrinks and grows as the flow develops (e.g. across a
shock).

#### Discontinuity handling (shock capturing)

The numerical flux at every face is selected per-run by the case file's
`flux_scheme` key, each scheme living in its own header:

- **`rusanov`** (default) — **Rusanov (local Lax-Friedrichs) flux**
  ([`rusanov_flux`](../include/RusanovFlux.h)):

  ```
  F* = 0.5*(F(U_L) + F(U_R)) - 0.5*S_max*(U_R - U_L)
  ```

  The second term is an upwind-biased dissipation, sized to the fastest
  local wave speed, added on top of the central average of the two physical
  fluxes. This dissipation is what keeps the scheme stable and
  non-oscillatory across a shock or contact discontinuity. It applies the
  same dissipation to every wave family, including slow-moving
  contact/shear waves that don't need it, which smears them more than
  necessary.
- **`hllc`** — **HLLC flux** ([`hllc_flux`](../include/HllcFlux.h)),
  restoring the contact wave that Rusanov/HLL collapse into the shock
  dissipation, giving sharper contact-discontinuity/shear-layer resolution
  while remaining robust and positivity-preserving. Wave speeds use Davis'
  simple estimate.
- **`exact`** — **exact Riemann (Godunov) flux**
  ([`exact_riemann_flux`](../include/ExactRiemannFlux.h)), the exact
  solution to the local Riemann problem via a per-face Newton-Raphson
  pressure solve and wave-pattern sampling, including a closed-form vacuum
  path for strong double-rarefactions. Highest accuracy of the three, but
  a per-face Newton iteration makes it meaningfully more expensive -- see
  the cost caveat in
  [docs/hllc-and-exact-riemann-plan.md](hllc-and-exact-riemann-plan.md).
  Its Newton-Raphson tolerance/iteration cap are case-file keys
  (`exact_riemann_tol`/`exact_riemann_max_iter`, see the table above) --
  exposed for transparency into how the solver behaves, not because
  they're expected to need tuning.

See [docs/hllc-and-exact-riemann-plan.md](hllc-and-exact-riemann-plan.md)
for the full algorithm/derivation references (Toro 2009 and others) for
`hllc` and `exact`.

Either way, this per-face numerical flux is the solver's only
shock-capturing mechanism. There is:
- **No reconstruction/limiter (no MUSCL, no slope limiting)** — the scheme
  is first-order in space everywhere, not just near discontinuities.
- **No separate/tunable artificial viscosity term.** Rusanov's built-in
  dissipation is blanket (applied at every face, not just near a
  discontinuity), which is simple and robust but smears shocks over more
  cells than a higher-order/limited scheme would.
- **No positivity check** on density/pressure — `pressure()`
  ([EulerState.h](../include/EulerState.h)) will silently produce a
  negative or `NaN` pressure from a sufficiently violent state without
  flagging it as the specific cause; it will eventually show up as a NaN/Inf
  residual and stop the run, but not with a diagnosis of *why*.

A pressure-sensor-based, switched artificial-viscosity enhancement (only
adding extra dissipation near large pressure jumps, rather than blanket
Rusanov dissipation everywhere) was discussed but **not implemented** — see
[docs/euler-artificial-viscosity.md](euler-artificial-viscosity.md) for the
design discussion and recommended diagnose-before-fixing approach if a run
diverges on a sharp discontinuity. (A MUSCL-reconstruction + limiter
approach, discussed in
[docs/dns-higher-order-scheme-plan.md](dns-higher-order-scheme-plan.md),
would supersede this from the opposite direction — a limiter that activates
only near real gradients, rather than blanket Rusanov dissipation — if
either is ever picked up, read them together.)

### Advection-diffusion of a passive scalar (`equation = advection_diffusion`)

[`AdvectionDiffusionFVMSolver`](../src/AdvectionDiffusionFVMSolver.cpp)
solves `dPhi/dt + div(u*Phi) = alpha*Laplacian(Phi)`, where `u = (u_adv,
v_adv)` is uniform and prescribed (not solved for — there is no momentum
equation here). The advective term is first-order upwind; the diffusive
term uses the same non-orthogonality-corrected face gradient as
Navier-Stokes' viscous term — see
[Gradient reconstruction](#gradient-reconstruction-shared-infrastructure).
`dt` is fixed and unchecked, the same footgun as scalar diffusion above.

### Compressible Navier-Stokes (`equation = navier_stokes`)

[`NavierStokesFVMSolver`](../src/NavierStokesFVMSolver.cpp) solves the
compressible Euler equations (reusing `EulerFVMSolver`'s exact inviscid flux
machinery — same `NumericalFluxScheme` choices, unchanged) plus a Newtonian
viscous stress tensor and Fourier heat conduction:

```
tau_xx = mu*(2*du/dx - (2/3)*(du/dx + dv/dy))
tau_yy = mu*(2*dv/dy - (2/3)*(du/dx + dv/dy))
tau_xy = mu*(du/dy + dv/dx)
q_x = -k*dT/dx, q_y = -k*dT/dy     (k = mu*cp/Pr, cp = gamma*R/(gamma-1))
```

built from the corrected face gradients (see
[Gradient reconstruction](#gradient-reconstruction-shared-infrastructure))
of the primitive velocity components and temperature (`T = p/(rho*R)`,
`R` = `gas_constant` — the first place in this codebase that needed a gas
constant at all, since Euler's own physics never did). `dt` is adaptive,
extending `EulerFVMSolver::compute_dt`'s inviscid CFL term with the standard
explicit-viscous-diffusion limit `2*nu/length` (`nu = mu/rho`), where
`length` is the face-normal-projected cell-centroid separation across each
face (see [Gradient reconstruction](#gradient-reconstruction-shared-infrastructure)'s
`face_normal_distance()`) — not a plain volume/area proxy, which
understates viscous stiffness on anisotropic (e.g. boundary-layer-clustered)
cells; see [Known limitations](#known-limitations).

This is a genuinely **laminar** solver — there is no turbulence closure of
any kind (no RANS Reynolds-averaging + model, no LES filtering + subgrid
model). It is *capable*, in the sense that it solves the full unmodeled
equations, of producing a Direct Numerical Simulation, but only within this
project's inherent 2D-only limitation — see
[Known limitations](#known-limitations) before treating any output as
representative of real (3D) turbulent flow. `resolution_report_file` (see
the case-file table above) reports a resolution-adequacy heuristic for
treating a run as fully-resolved 2D unsteady, computed by
`NavierStokesFVMSolver::compute_resolution_diagnostics()`.

Viscous flux is exactly zero at `ns_farfield`/`ns_outflow` boundaries (not
extrapolated) — the intended use is placing those far enough from any wall
that this doesn't matter. **`ns_outflow` is a zero-gradient/do-nothing
condition, not a true periodic one** — nothing in it forces the two ends of
a multi-cell-wide domain to match, and a slow, non-decaying drift mode can
develop and persist indefinitely if a case relies on it for translational
invariance (found and documented in detail in
[docs/navier-stokes-tracker.md](navier-stokes-tracker.md) Phase 4 — the
Couette validation's own domain is exactly one cell wide in the flow
direction specifically to avoid this).

See [docs/navier-stokes-tracker.md](navier-stokes-tracker.md) for the full
phase-by-phase development history (gradient reconstruction, non-orthogonal
face correction, advection-diffusion, the Navier-Stokes viscous terms, and
the Couette flow validation), including two genuine dead ends found and
corrected along the way — worth reading before assuming a superficially
similar test setup will behave the same way in a new context.

#### Wall diagnostics

[`WallTraction.h`](../include/WallTraction.h)/[`.cpp`](../src/WallTraction.cpp)
add per-`ns_wall*`-patch skin-friction/pressure coefficients, drag/lift/
moment integration, wall-normal `y+`, and boundary-layer thickness estimates
as a separate, on-demand **post-processing pass** over `NavierStokesFVMSolver`
— not a change to `step()`'s hot loop — following the same "recompute fresh,
retain nothing from `step()`" contract as
`compute_resolution_diagnostics()`. See
[docs/wall-diagnostics-plan.md](wall-diagnostics-plan.md) for the full
design discussion; this section summarizes the result.

**Phasing status**: Phase 1 (forces/`Cf`/`Cp`/`Cm`/`y+`), Phase 2
(boundary-layer thickness), and Phase 3 (the point-location boundary-layer
alternative, `compute_boundary_layer_profiles_point_location()` below, built
and verified in direct response to Phase 2's own stress-test finding) are all
implemented and verified. Wiring `RANSFVMSolver` into `CaseInput`/`main.cpp`
so it could adopt these same wall diagnostics was a separate, larger,
pre-existing gap outside this plan's own original scope — that gap is now
closed too (see
[RANS (Spalart-Allmaras) turbulence closure](#rans-spalart-allmaras-turbulence-closure)).

- **`compute_wall_traction()`** samples every `NoSlipWall` boundary face: the
  owning cell's own pressure (a first-order approximation — neither solver
  retains a separate boundary pressure state, the same convention as
  `AdvectionDiffusionFVMSolver`'s Neumann-boundary gradient stencil value),
  density, effective viscosity, and a signed wall-tangential shear stress
  `tau_wall` built from the same corrected-face-gradient stress-tensor
  assembly `step()` itself uses. `Cf`, `Cp`, and `y+` are derived from this
  per-face sample plus a `WallReferenceQuantities` (freestream density/
  velocity/pressure, a length scale, and a moment reference point) —
  reference-agnostic by design, which is exactly what let `RANSFVMSolver`
  adopt the same struct with `mu + rho*nu_t` in place of `mu` (see
  [RANS (Spalart-Allmaras) turbulence closure](#rans-spalart-allmaras-turbulence-closure)).
- **`compute_wall_forces()`** integrates friction + pressure force (and
  moment about `moment_reference_x/y`) per patch, plus a domain total, and
  nondimensionalizes into `Cd`/`Cl`/`Cm`. Validated to near-machine precision
  against planar Couette flow's exact `tau_wall = mu*U/H` (`--verify-wall-forces`)
  — `Cp`/moment use a looser (5%) tolerance there, not because the plumbing
  is imprecise, but because that test's deliberately small reference
  velocity (needed to keep the compressible solver close to the
  incompressible analytic solution) amplifies a known, pre-existing small
  secular-drift artifact of that exact test setup (adiabatic walls + `ns_outflow`'s
  imperfect periodicity) into a much larger relative `Cp` error — see that
  function's own comment in [main.cpp](../src/main.cpp).
- **`compute_boundary_layer_profiles()`** estimates `delta_99`,
  displacement thickness, momentum thickness, and shape factor by marching
  cell-to-cell away from a wall face along whichever neighbor's face normal
  is most aligned with the (fixed, per-march) wall-normal direction,
  sampling tangential speed at each visited cell centroid. **This assumes a
  locally wall-normal-ish structured cell stacking near the wall** (true for
  boundary-layer-clustered meshes, reproducing `--verify-flat-plate`'s
  hardcoded column-walk result to under 1%) **and this is a demonstrated
  limitation, not a hypothetical one**: `--verify-bl-marching-unstructured`
  shows a systematic ~20-25% `delta_99` error on a genuinely unstructured
  (checkerboard-triangulated) near-wall mesh, root-caused to
  `face_normal_distance()`'s per-crossed-face-normal projection not matching
  the fixed global march direction once a step crosses a diagonal face —
  see [docs/wall-diagnostics-plan.md](wall-diagnostics-plan.md)'s "Phase 2
  stress-test finding" for the full mechanism.
- **`compute_boundary_layer_profiles_point_location()`** is the point-location
  alternative built specifically to not share marching's failure mode above:
  for each wall face, it walks straight out along the fixed inward normal at
  `n_samples` evenly-spaced offsets up to `max_distance`, brute-force
  point-locating (point-in-polygon over every cell, no spatial acceleration
  structure — same "simple first" precedent as `WallDistance.h`'s brute-force
  point-to-segment search) which cell contains each sample point and reading
  that cell's velocity there — no near-wall-topology dependence at all, at
  the cost of an `O(cells)` search per sample point (meaningfully more
  expensive than marching's `O(1)`-ish per step on a large mesh).
  `--verify-bl-point-location` reruns marching's own failing stress-test mesh
  through this method and confirms it reproduces the exact `delta_99` to
  within 0.76%. Selected via the `boundary_layer_method` case-file key
  (`marching` default, or `point-location`); `boundary_layer_max_distance`
  (default: automatic, the mesh's own bounding-box diagonal) and
  `boundary_layer_n_samples` (default 200) tune it further — see the
  case-file table above.
- **Face-to-node averaging** (`average_wall_samples_to_nodes()`) reduces
  every per-face quantity onto actual wall mesh nodes (plain average for an
  interior node shared by two wall faces on the same patch; direct value at
  an open-patch endpoint; reported once per patch for a corner node shared
  by two different patches) before it reaches `wall_profile_file`, so a
  plotted surface distribution lines up with mesh-vertex coordinates rather
  than face-midpoint samples offset from them.
- `wall_forces_file` and `wall_profile_file` are Navier-Stokes/RANS-only
  case-file keys (see the tables above), wired in `run_navier_stokes()`/
  `run_rans()` following `resolution_report_file`'s exact existing precedent
  (`ensure_parent_directory`, `scratch_dir` rebasing, opened once/appended
  every N steps for `wall_forces_file`; `wall_profile_file` instead follows
  `output_file`/`write_interval`'s convention — a numbered snapshot every N
  steps if `wall_profile_interval > 0`, and the plain path always
  (re)written once at the run's natural end, including on divergence, for
  the same "still written for inspection" reason `output_file` is).
  `RANSFVMSolver::compute_wall_traction_samples()`/
  `compute_boundary_layer_profile_samples[_point_location]()` mirror
  `NavierStokesFVMSolver`'s identically-named methods exactly, substituting
  `mu + rho*nu_t` for `mu` as the effective viscosity fed to
  `compute_wall_traction()`. Note `RANSFVMSolver` has no
  `resolution_report_file`-equivalent diagnostic — that key stays
  Navier-Stokes-only.
- `--verify-flat-plate`'s Blasius/Pohlhausen extension predates this wiring
  and still calls `compute_wall_traction()`/`compute_boundary_layer_profiles()`
  directly from the verification harness as a one-off, rather than through
  `RANSFVMSolver`'s member methods — both paths reuse the same underlying
  `WallTraction.h` functions, so this is a difference in call site, not in
  what's actually computed.

### RANS (Spalart-Allmaras) turbulence closure

[`RANSFVMSolver`](../include/RANSFVMSolver.h) is a fifth, independent solver
class — not a toggle on `NavierStokesFVMSolver`, per this project's existing
"one class per equation set" pattern (see
[docs/archive/rans-spalart-allmaras-tracker.md](archive/rans-spalart-allmaras-tracker.md)'s
architecture decision). It duplicates `NavierStokesFVMSolver`'s inviscid
flux, ghost-state, and `compute_dt()` structure unchanged, and adds:

- One extra transported scalar, `nut` ("nu-tilde"), advected by the coupled
  mean-flow velocity and diffused via the same non-orthogonality-corrected
  `face_gradient()` every other diffusive term in this codebase uses, with
  production/destruction/cross-diffusion source terms from the
  Spalart-Allmaras-noft2 model (`sa_fv1`/`sa_eddy_viscosity`/
  `compute_sa_source_terms()` in
  [SpalartAllmaras.h](../include/SpalartAllmaras.h)/[.cpp](../src/SpalartAllmaras.cpp)),
  including the standard negative-`S~` robustness fix (`cv2`/`cv3`). The
  model's nine constants (`cb1`/`cb2`/`sigma`/`kappa`/`cw2`/`cw3`/`cv1`/`cv2`/`cv3`,
  `cw1` derived from the other three) live in a `SAModelConstants` struct
  rather than fixed globals, so the `sa_*` case-file keys (see the case-file
  table above) can override them — the defaults are the standard SA-noft2
  set, and overriding them is not meant to be normal practice.
- A wall-distance field (`compute_wall_distance()` in
  [WallDistance.h](../include/WallDistance.h)/[.cpp](../src/WallDistance.cpp)),
  precomputed once at construction from whichever faces' patch is
  `NSBoundaryType::NoSlipWall` — a brute-force point-to-segment distance,
  `O(cells * wall_faces)`, no spatial acceleration structure.
- Molecular viscosity/conductivity in the mean-flow viscous terms replaced by
  turbulence-inclusive effective values, `mu_eff = mu + rho*nu_t` and
  `k_eff = cp*(mu/Pr + rho*nu_t/Pr_t)` (`Pr_t` = turbulent Prandtl number,
  `nu_t = nut * fv1`).
- `RANSBoundaryCondition` wraps `NSBoundaryCondition` unchanged (a no-slip
  wall/farfield/outflow means the same thing for the mean-flow equations
  here) and adds only `farfield_nut` — the freestream `nut` a `Farfield`
  boundary prescribes; SA is documented as sensitive to this choice.
  `NoSlipWall`'s `nut` boundary value is exactly `0` by the model's own
  definition (not case-configurable); `Outflow`'s is a zero-order
  extrapolation.
- `compute_dt()`'s viscous stability term uses `nu + nu_t` (the larger of the
  two adjoining cells' `nu_t`) rather than just molecular `nu` — `nu_t` can
  be 10-1000x molecular `nu` in a real turbulent boundary layer.

**Reachable from a case file via `equation = rans`.** `run_rans()` in
[main.cpp](../src/main.cpp) drives it the same way `run_navier_stokes()`
drives `NavierStokesFVMSolver` — same checkpointing (`CheckpointEquation::RANS`
round-trips both the mean-flow state and `nut` in one file; see
[Checkpointing & restart](#checkpointing--restart)), residual tracking (a
`residual_nut` column alongside `rho`/`rho_u`/`rho_v`/`E`, plus an optional
`residual_tolerance_nut` stopping criterion), and wall diagnostics (see
[Wall diagnostics](#wall-diagnostics)) as the other solvers get. VTK output
(`write_rans_fields()`) adds `nut` and the derived `nu_t` to the usual
`rho`/`u`/`v`/`p`/`mach`/`T` fields. See the RANS-only case-file table above
for every `rans_*`/`sa_*`/`prandtl_t`/`initial_nut` key this adds over
Navier-Stokes. `--verify-rans-stability`/`--verify-flat-plate` (see
[Running](#running)) predate this wiring and still construct
`RANSFVMSolver` directly rather than through a case file — they remain this
project's correctness gate for the solver itself, independent of the
case-file plumbing around it.

**Verification target was relaxed, not fully met.** The original goal (a
flat-plate turbulent boundary layer matching the log-law velocity profile)
was not reached at this project's tractable, explicit, 2D compressible
time-stepping scale — two real attempts (`Re_L = 1e4, 1e5`) both showed SA
correctly predicting no sustained turbulence at those Reynolds numbers, not
a numerical failure. `--verify-flat-plate`'s actual check compares against
the laminar Pohlhausen profile instead. See
[docs/archive/rans-spalart-allmaras-tracker.md](archive/rans-spalart-allmaras-tracker.md)
Phase 4 for the full cost/benefit reasoning behind that decision before
assuming this solver validates real turbulent-flow predictions.

## Gradient reconstruction (shared infrastructure)

[`GradientCalculator`](../include/GradientReconstruction.h) reconstructs
per-cell gradients of an arbitrary scalar field, used by both
advection-diffusion's diffusive term and Navier-Stokes' viscous terms (once
per velocity component and temperature, sharing one `GradientCalculator`
instance since its precomputed geometry-dependent setup doesn't depend on
which field is being differentiated). Selected per-run by the case file's
`gradient_scheme` key:

- **`least-squares`** (default) — weighted least-squares fit against
  neighboring cell/boundary-face values, weighted `1/distance^2`. **Exact
  for a linear field on any non-degenerate mesh**, by construction of its
  normal equations — this holds regardless of mesh skewness/orthogonality.
- **`green-gauss`** — divergence-theorem sum of face-interpolated values
  over each cell's bounding faces. Exact for a linear field **only** when
  every face midpoint lies on the line joining its two cell centroids (true
  on a Cartesian/orthogonal mesh) — has a real, sometimes large,
  non-orthogonality error otherwise (empirically ~7% on a moderately
  skewed test mesh; see `--verify-gradient`).

A second function, `face_gradient()`, computes a non-orthogonality-corrected
face-normal derivative (the standard Jasak 1996 over-relaxed decomposition:
a direct two-point term rescaled by `1/cosTheta`, plus a correction term
using the interpolated cell-gradient vector) — this is what actually feeds
the diffusive/viscous flux, not the raw cell gradient. A naive average of
two cell gradients (or a naive two-point difference ignoring
non-orthogonality) risks decoupling the flux from `phi_R - phi_L` entirely
(checkerboard oscillation risk) or has a real, potentially large error on a
skewed mesh (found empirically to be *larger* than the exact analytic
normal derivative's own magnitude on one deliberately-adversarial test
mesh) — see `--verify-gradient`'s "corrected vs naive" comparison and
[docs/navier-stokes-tracker.md](navier-stokes-tracker.md) Phase 1 for the
full derivation and verification.

For a vector field's off-diagonal/tangential terms (needed by the viscous
stress tensor, not by a scalar diffusion flux), `face_gradient()`'s
correction to the NORMAL component only is combined with the plain
interpolated average gradient to get the full 2D gradient vector — see
`corrected_face_gradient_vector()` in
[NavierStokesFVMSolver.cpp](../src/NavierStokesFVMSolver.cpp).

A third function, `face_normal_distance()`, returns the face-normal-
projected distance between a face's two cell centroids (or a cell centroid
and the face midpoint, for a boundary face) — the same `dist * cosTheta`
quantity `face_gradient()` divides by internally, exposed standalone so
`NavierStokesFVMSolver`/`RANSFVMSolver`'s `compute_dt()` can use the same
length scale for their viscous stability estimate that the flux itself
uses (see [docs/navier-stokes-tracker.md](navier-stokes-tracker.md) Phase 6).

**No near-singular/degenerate-neighbor guard** — consistent with this
project's "no mesh-quality safeguards" stance (see
[Mesh formats](#mesh-formats)), a cell whose neighbor points are nearly
collinear could produce a near-singular least-squares normal-equation
matrix with no diagnostic.

## Boundary conditions

Applied by patch name (matched against `UnstructuredMesh::patches`, which
comes from the mesh file — see [Mesh formats](#mesh-formats)) in
[`run_diffusion`/`run_advection_diffusion`/`run_euler`/`run_navier_stokes`/`run_rans`](../src/main.cpp),
not by the mesh reader.

- **Diffusion / advection-diffusion**: `BoundaryType::Dirichlet` (pin phi)
  or `::Neumann` (prescribe flux) — see
  [UnstructuredMesh.h](../include/UnstructuredMesh.h). Identical semantics
  for both equation sets, so this one enum is shared rather than duplicated.
- **Euler**: `EulerBoundaryType::Wall` (mirror normal velocity, zero net
  normal flux), `::Farfield` (fixed prescribed conserved state, converted
  from primitive once at setup), `::Outflow` (zero-gradient extrapolation)
  — see [EulerFVMSolver.h](../include/EulerFVMSolver.h) /
  [.cpp](../src/EulerFVMSolver.cpp)'s `ghost_state`.
- **Navier-Stokes**: `NSBoundaryType::NoSlipWall` (mirrors BOTH velocity
  components about the wall's own `(wall_u, wall_v)` — not just the normal
  component mirrored about zero, unlike Euler's slip wall — plus an
  associated thermal condition: isothermal via `wall_temperature`, or
  adiabatic = zero heat flux prescribed directly), `::Farfield`,
  `::Outflow` (same meaning as Euler's, but a distinct enum/spec vector
  since a no-slip wall needs data an inviscid wall never does) — see
  [NavierStokesFVMSolver.h](../include/NavierStokesFVMSolver.h) /
  [.cpp](../src/NavierStokesFVMSolver.cpp)'s `ghost_state`. Deliberately its
  own type, not a Euler `EulerBoundaryType` extension.
- **RANS**: `RANSBoundaryCondition` wraps an `NSBoundaryCondition` unchanged
  (same `NoSlipWall`/`Farfield`/`Outflow` meanings as Navier-Stokes above)
  and adds only `farfield_nut` (the transported `nut` scalar's prescribed
  freestream value at a `Farfield` boundary) — see
  [RANSFVMSolver.h](../include/RANSFVMSolver.h) and the `rans_*` boundary
  keywords in the case-file reference above.

## Output format

[`VtkWriter`](../src/VtkWriter.cpp) writes legacy ASCII VTK
`DATASET UNSTRUCTURED_GRID`: `POINTS` (node coordinates, `z=0`), `CELLS` +
`CELL_TYPES` (every cell emitted as generic `VTK_POLYGON`, so an arbitrary
vertex-count cell round-trips correctly), then `CELL_DATA` scalar field(s).
Openable directly in ParaView. The file's title line (line 2 of the legacy
VTK header) is `FiniteVolume <version>` (e.g. `FiniteVolume 0.2.42`),
stamping every result with the exact solver build that produced it. Every
written double (coordinates and field values) uses the case file's
`output_precision` (default 6, valid range 1-17 significant digits).

- **Diffusion / advection-diffusion** writes one scalar field, `phi`.
- **Euler** writes five derived primitive fields per cell: `rho`, `u`, `v`,
  `p`, `mach` (see `write_euler_fields` in [main.cpp](../src/main.cpp)) —
  note these are *derived* from the conserved state at write time, not the
  solver's native representation.
- **Navier-Stokes** writes the same five fields plus `T` (temperature, see
  `write_navier_stokes_fields` in [main.cpp](../src/main.cpp) and
  `temperature()` in [EulerState.h](../include/EulerState.h)).

`write_interval > 0` in the case file additionally writes numbered snapshots
(`numbered_filename` in [main.cpp](../src/main.cpp) inserts a zero-padded
step number before `output_file`'s extension, e.g.
`output/result_000100.vtk`) throughout the run, independent of the final
`output_file` write.

`output_file`'s parent directory (and `checkpoint_file`'s, `residual_file`'s,
and `resolution_report_file`'s) is created automatically if it doesn't already exist
(`ensure_parent_directory` in [main.cpp](../src/main.cpp)), so a long run
doesn't fail at its very last step over a missing folder. If `scratch_dir`
is set, any of those four paths given as *relative* rebases under it first
(`resolve_output_path` in [main.cpp](../src/main.cpp)); `mesh_file` and
absolute paths are never affected.

## Monitoring & stopping criteria

A run stops (writing final output either way) on the first of:
1. **Divergence** — any residual component is NaN/Inf. Immediate stop, a
   final VTK is still written for post-mortem inspection, but **no
   checkpoint is written**, and the process exits non-zero.
2. **Convergence** — `residual_tolerance` (diffusion / advection-diffusion)
   or every `residual_tolerance_*` the user set (Euler / Navier-Stokes /
   RANS, the latter adding `residual_tolerance_nut`) is satisfied.
3. **Ctrl+C** — a clean stop after the in-flight step.
4. **`nsteps` reached** — the "natural" end.

`residual_file`/`residual_interval` write a CSV of residual norms (one
column for diffusion/advection-diffusion, four — `rho`, `rho_u`, `rho_v`,
`E` — for Euler/Navier-Stokes, or five with `nut` appended for RANS)
throughout the run, independent of whether any stopping tolerance is set;
useful for judging convergence behavior even on a fixed-`nsteps` run.

`resolution_report_file`/`resolution_report_interval` (Navier-Stokes only)
similarly write a CSV, but it is **purely a diagnostic — it never affects
stopping** (there is no `resolution_tolerance` key); see
[Gradient reconstruction](#gradient-reconstruction-shared-infrastructure)
and [Known limitations](#known-limitations) for what it does and does not
mean.

## Checkpointing & restart

[`Checkpoint`](../include/Checkpoint.h) writes a simple fixed-layout binary
file (magic + format version + equation tag [`0` = Diffusion, `1` = Euler,
`2` = AdvectionDiffusion, `3` = NavierStokes, `4` = RANS] + cell count + step
index + build number + raw field doubles — **not portable across
architectures**, only meant to round-trip on the same machine). RANS's
payload is the mean-flow `EulerState` field immediately followed by `nut`,
both resumed together via `RANSFVMSolver::set_field()`; this didn't require
a format-version bump, since the fixed header is unchanged and the payload
shape has always been a function of the equation tag alone. If `checkpoint_file` is set
and that file already exists when a run starts, the run **auto-resumes from
it** instead of using the case file's initial condition — this is how a
stopped run (any of the four reasons above except divergence, which never
checkpoints) is continued: raise `nsteps` and/or loosen a tolerance, then
rerun the same case file unchanged. The resume message reports the build
number the checkpoint was written by (e.g. `Resuming from checkpoint at
step 500 (written by build 42)`), purely as traceability metadata.

Resuming validates that the checkpoint's **format version**, equation tag,
and cell count match the current run — the mesh you resume with is not
saved in the checkpoint and must be supplied consistently. The format
version (currently 2) gates compatibility; the embedded build number does
not — it's informational only. **A checkpoint written by format version 1
(before the build-number field existed) cannot be read by this version**;
there is no migration path, since checkpoints are meant to round-trip a
single in-progress run, not archive results long-term.

## Performance (OpenMP)

All five solvers parallelize their per-face and per-cell passes with
`#pragma omp parallel for` (split into two passes specifically to avoid a
scatter-add race across cells shared by two faces — see
[Solvers](#solvers)).
`num_threads` in the case file calls `omp_set_num_threads`; `0` (default)
leaves OpenMP's own default in effect, which respects `OMP_NUM_THREADS`.
There is a caching-friendliness note in
[UnstructuredMesh.h](../include/UnstructuredMesh.h) about cell/face memory
layout (RCM reordering) that has **not** been implemented — the mesh is
used in whatever order the reader produced it.

## Known limitations

- **2D only.** Not a partial limitation to be lifted later within this data
  model — `Node`, `Face` (a line segment), and the whole geometry pipeline
  are inherently 2D. **This means there is no true 3D DNS/turbulence
  capability, ever** — classical Kolmogorov scaling and the turbulent
  energy cascade are 3D results (the cascade mechanism, vortex stretching,
  doesn't exist the same way in 2D). `NavierStokesFVMSolver` can be pushed
  toward a fully-resolved 2D unsteady simulation, and
  `resolution_report_file` gives a resolution-adequacy heuristic for that,
  but "2D DNS" and "3D DNS" are not the same claim — don't let one imply
  the other in anything built on top of this. See
  [docs/navier-stokes-tracker.md](navier-stokes-tracker.md) Phase 5.
- **`NavierStokesFVMSolver` itself has no turbulence closure** — a
  turbulent-regime flow run on a mesh that doesn't resolve every scale gives
  a real result that is *mesh-dependent* (changes character as you refine),
  not a converged approximation. A Spalart-Allmaras RANS closure exists
  (`RANSFVMSolver`, no LES option) and **is reachable from a case file**
  (`equation = rans`) — see
  [RANS (Spalart-Allmaras) turbulence closure](#rans-spalart-allmaras-turbulence-closure).
  Its own verification target (a turbulent log-law boundary layer) was not
  reached at this project's tractable Reynolds numbers either — **don't treat
  it as a validated turbulent-flow predictor**; what's actually verified is
  that it stays numerically stable and reduces correctly to a laminar
  boundary layer when SA predicts no sustained turbulence, not that its
  turbulence-sustaining behavior matches a real turbulent flow.
- **No mesh-quality checks** — see [Mesh formats](#mesh-formats). A
  degenerate mesh fails silently (NaN/Inf), not with a diagnostic. The same
  applies to `GradientCalculator`'s Least-Squares matrix inversion (no
  near-singular/degenerate-neighbor guard) and `NSBoundaryType`'s viscous
  face-normal computation (no near-zero-`cosTheta` guard).
- **Gmsh reader is triangle/quad/line only** — a parser limitation, not a
  solver one; use `.fvmesh` for anything with more faces per cell.
- **Euler solver is first-order, no limiter, no tunable artificial
  viscosity** — see [Discontinuity handling](#discontinuity-handling-shock-capturing).
  Navier-Stokes' inviscid part inherits the same limitation (it reuses this
  machinery unchanged); see
  [docs/dns-higher-order-scheme-plan.md](dns-higher-order-scheme-plan.md)
  for a deferred discussion of what fixing this would take, motivated by
  DNS accuracy but relevant to shock-capturing too.
- **`dt` for diffusion / advection-diffusion is fixed and unchecked** — an
  unstable `dt` diverges rather than being rejected up front.
- **`NSBoundaryType::Outflow` is zero-gradient/do-nothing, not truly
  periodic** — do not rely on it to enforce translational invariance across
  a domain more than one cell wide in that direction; a slow, non-decaying
  drift mode can develop and persist indefinitely. Found and worked around
  (not fixed at the BC level) in
  [docs/navier-stokes-tracker.md](navier-stokes-tracker.md) Phase 4.
- **Viscous flux is forced to exactly zero at `ns_farfield`/`ns_outflow`
  boundaries**, not extrapolated — correct when those are placed away from
  any wall, silently wrong if a case ever puts one close to a real
  velocity/temperature gradient.
- **`gas_constant`/`prandtl` are unvalidated** — a zero or negative value
  silently produces `Inf`/division-by-zero rather than a diagnostic.
- **A strongly anisotropic mesh (e.g. boundary-layer-clustered) takes a
  smaller adaptive `dt`** in `NavierStokesFVMSolver`/`RANSFVMSolver` than a
  mesh with similar-sized cells, since `compute_dt()`'s viscous stability
  term is now direction-aware rather than a plain volume/area average (see
  [docs/navier-stokes-tracker.md](navier-stokes-tracker.md) Phase 6) — more
  steps for the same physical time is the expected cost of the more
  conservative (and more correct) estimate, not a regression.
- **`Farfield`'s ghost state has no characteristic branching by local Mach
  number** — it returns the prescribed state unconditionally, a fixed-state
  Dirichlet condition rather than a proper subsonic/supersonic
  characteristic-based farfield BC. Not confirmed to cause any specific
  divergence on its own, but a real over-specification risk near Mach 1
  (one characteristic should be outgoing, not fixed); see
  [docs/ns-cfl-margin-and-farfield-bc-findings.md](ns-cfl-margin-and-farfield-bc-findings.md).
  No decision has been made on whether to fix this.
- **Checkpoints are not portable** across machines/architectures.
- **Wall pressure for `wall_forces_file`/`wall_profile_file`/`Cp` is a
  first-order approximation** — the owning cell's own pressure, not a
  boundary-extrapolated value (neither solver retains a separate boundary
  pressure state). See [Wall diagnostics](#wall-diagnostics).
- **`compute_boundary_layer_profiles()`'s cell-marching heuristic assumes a
  locally wall-normal-ish structured cell stacking near the wall, and this
  is demonstrated, not hypothetical** — `--verify-bl-marching-unstructured`
  shows a systematic ~20-25% `delta_99` error on a genuinely unstructured
  near-wall triangulation. Does not affect `Cf`/`Cp`/`Cd`/`Cl`/`Cm`/`y+`,
  which don't use this function at all. A point-location alternative,
  `compute_boundary_layer_profiles_point_location()`, does not share this
  limitation (verified on the same failing mesh by
  `--verify-bl-point-location`) at the cost of being meaningfully more
  expensive on a large mesh — select it via the `boundary_layer_method`
  case-file key when the default (`marching`) is a concern. See
  [Wall diagnostics](#wall-diagnostics) and
  [docs/wall-diagnostics-plan.md](wall-diagnostics-plan.md).

## File map

| Concern | Header | Implementation |
|---|---|---|
| Mesh data structures | [UnstructuredMesh.h](../include/UnstructuredMesh.h) | (header-only geometry) |
| Gmsh + FVMESH readers | [MeshReader.h](../include/MeshReader.h) | [MeshReader.cpp](../src/MeshReader.cpp) |
| FVMESH writer | [FvMeshWriter.h](../include/FvMeshWriter.h) | [FvMeshWriter.cpp](../src/FvMeshWriter.cpp) |
| Case file parsing | [CaseInput.h](../include/CaseInput.h) | [CaseInput.cpp](../src/CaseInput.cpp) |
| Diffusion solver | [UnstructuredFVMSolver.h](../include/UnstructuredFVMSolver.h) | [UnstructuredFVMSolver.cpp](../src/UnstructuredFVMSolver.cpp) |
| Advection-diffusion solver | [AdvectionDiffusionFVMSolver.h](../include/AdvectionDiffusionFVMSolver.h) | [AdvectionDiffusionFVMSolver.cpp](../src/AdvectionDiffusionFVMSolver.cpp) |
| Euler solver + physics | [EulerFVMSolver.h](../include/EulerFVMSolver.h), [EulerState.h](../include/EulerState.h) | [EulerFVMSolver.cpp](../src/EulerFVMSolver.cpp) |
| Navier-Stokes solver | [NavierStokesFVMSolver.h](../include/NavierStokesFVMSolver.h) | [NavierStokesFVMSolver.cpp](../src/NavierStokesFVMSolver.cpp) |
| RANS (Spalart-Allmaras) solver | [RANSFVMSolver.h](../include/RANSFVMSolver.h) | [RANSFVMSolver.cpp](../src/RANSFVMSolver.cpp) |
| SA turbulence model source terms | [SpalartAllmaras.h](../include/SpalartAllmaras.h) | [SpalartAllmaras.cpp](../src/SpalartAllmaras.cpp) |
| Wall-distance module (RANS) | [WallDistance.h](../include/WallDistance.h) | [WallDistance.cpp](../src/WallDistance.cpp) |
| Wall diagnostics (forces/Cf/Cp/y+/boundary-layer thickness) | [WallTraction.h](../include/WallTraction.h) | [WallTraction.cpp](../src/WallTraction.cpp) |
| Gradient reconstruction (shared) | [GradientReconstruction.h](../include/GradientReconstruction.h) | [GradientReconstruction.cpp](../src/GradientReconstruction.cpp) |
| Checkpoint/restart | [Checkpoint.h](../include/Checkpoint.h) | [Checkpoint.cpp](../src/Checkpoint.cpp) |
| VTK output | [VtkWriter.h](../include/VtkWriter.h) | [VtkWriter.cpp](../src/VtkWriter.cpp) |
| CLI / run loop / entry point | — | [main.cpp](../src/main.cpp) |
| FVMESH grammar spec | [docs/fvmesh-format.md](fvmesh-format.md) | — |
| Deferred Euler viscosity discussion | [docs/euler-artificial-viscosity.md](euler-artificial-viscosity.md) | — |
| Navier-Stokes development history (Phases 0-6) | [docs/navier-stokes-tracker.md](navier-stokes-tracker.md) | — |
| RANS (Spalart-Allmaras) implementation history (archived; all phases closed — `CaseInput`/`main.cpp` wiring landed afterward, see [Solvers](#solvers)) | [docs/archive/rans-spalart-allmaras-tracker.md](archive/rans-spalart-allmaras-tracker.md) | — |
| Deferred higher-order time/space scheme discussion | [docs/dns-higher-order-scheme-plan.md](dns-higher-order-scheme-plan.md) | — |
| Wall diagnostics design + Phase 2 stress-test finding | [docs/wall-diagnostics-plan.md](wall-diagnostics-plan.md) | — |
