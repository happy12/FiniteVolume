# Navier-Stokes implementation tracker

Tracks progress toward extending this solver to the compressible
Navier-Stokes equations, per the staged plan agreed with Mathieu
(2026-07-05). Each phase gets its own verification gate before the next one
starts; this file records what was added, what was actually tested, and
what's still a known gap or open question when a phase closes.

## Phase status

| # | Phase | Status |
|---|---|---|
| 0 | Cell-gradient reconstruction (Green-Gauss / Least-Squares) + verification | Done |
| 1 | Non-orthogonal corrected face gradient | Done |
| 2 | `AdvectionDiffusionFVMSolver` (new `EquationSet::AdvectionDiffusion`) | Done |
| 3 | `NavierStokesFVMSolver` (new `EquationSet::NavierStokes`, viscous stress + heat conduction) | Done |
| 4 | Couette flow validation case | Done |
| 5 | DNS exposure: resolution diagnostic (Tier 1) + higher-order scheme plan (Tier 2) | Done |
| 6 | `compute_dt()` viscous length-scale hardening (real-mesh CFL divergence found via Bravo integration) | Done |

Architecture decisions locked in for phases 2-3 (confirmed with Mathieu,
2026-07-05): both get their own solver class rather than extending
`UnstructuredFVMSolver`/`EulerFVMSolver`, matching this project's existing
"one class per equation set" pattern. `AdvectionDiffusionFVMSolver` reuses
`BoundaryType`/`BoundaryPatch` as-is (Dirichlet/Neumann on a scalar is
identical semantics to the diffusion solver). `NavierStokesFVMSolver` reuses
`EulerState.h`'s conserved-state representation and inviscid `flux()` for
its inviscid part, but needs its own `NSBoundaryCondition`/`NSBoundaryType`
(a viscous no-slip wall is genuinely different from Euler's slip wall).
Gradient scheme default is Least-Squares.

---

## Phase 0 — Cell-gradient reconstruction

**Added:**
- [include/GradientReconstruction.h](../include/GradientReconstruction.h) /
  [src/GradientReconstruction.cpp](../src/GradientReconstruction.cpp):
  `GradientScheme` (`GreenGauss` / `LeastSquares`), `Gradient2`, and
  `GradientCalculator` -- precomputes each cell's Least-Squares
  normal-equation matrix (inverted) once at construction from mesh geometry
  alone; `compute()` reconstructs gradients from a field in one pass, race-free
  under OpenMP (each cell writes only its own output entry).
- `--verify-gradient` CLI mode in [src/main.cpp](../src/main.cpp)
  (`build_skewed_verification_mesh()` + `run_verify_gradient()`): builds a
  small deterministic, self-contained test mesh with no external file
  dependency.

**Test performed:** imposed `phi(x,y) = 1.7 + 2.3*x - 1.1*y` at every cell
centroid and boundary-face midpoint of the test mesh, reconstructed the
gradient both ways, and compared against the exact analytic gradient
`(2.3, -1.1)`.

**Finding (a real bug caught before it shipped):** the first version of the
test mesh used a pure shear (`x = i*h + shear*j*h, y = j*h`) rather than a
genuinely non-affine distortion. A shear is an *affine* map of a Cartesian
grid, and affine maps preserve collinearity -- since every face midpoint
lies on the line joining its two cell centroids on any Cartesian grid, it
still does after a shear, so Green-Gauss came out exact "by accident"
(error ~1e-14, same as Least-Squares) and the test proved nothing about
non-orthogonality. Replaced it with a smooth "wavy" perturbation
(`x = X + w*sin(pi*X)*sin(pi*Y)`, `y = Y + w*sin(pi*X)*sin(pi*Y)`, boundary
untouched since the perturbation vanishes at X/Y = 0 or 1), which is
genuinely non-affine and does break centroid/face-midpoint collinearity.

**Result on the corrected (wavy) mesh:**
- Least-Squares: max gradient error ~9.3e-15 (machine precision), matching
  the theoretical expectation that it's exact for a linear field on any
  non-degenerate point distribution.
- Green-Gauss: max gradient error ~0.073 (7.3%) -- a real, expected
  non-orthogonality error, not a bug. This is the reason Phase 1 exists.

**Known limitations / open items:**
- `GradientCalculator`'s Least-Squares matrix inversion has no
  near-singular/degenerate guard (consistent with this project's existing
  "no mesh-quality safeguards anywhere" stance -- see
  [docs/MANUAL.md](MANUAL.md#known-limitations) -- but worth remembering if
  a future mesh has a cell whose neighbor points are nearly collinear).
- `boundary_phi` is a caller-filled `mesh.faces`-sized array, not routed
  through `BoundaryType`/`BoundaryPatch` -- deliberate, so the module stays
  usable for velocity/temperature later without hardcoding scalar-diffusion
  BC semantics, but it means any real solver adopting this module is
  responsible for filling it in correctly per time step.
- The verification mesh/parameters (n=8, wobble=0.15, linear coefficients)
  are hardcoded in `run_verify_gradient()`, not CLI-configurable -- matches
  `--version`'s zero-argument simplicity; nothing in the current design
  needs it configurable yet.

---

## Phase 1 — Non-orthogonal corrected face gradient

**Added:**
- `FaceGradient` struct + `face_gradient()` free function in
  [include/GradientReconstruction.h](../include/GradientReconstruction.h) /
  [src/GradientReconstruction.cpp](../src/GradientReconstruction.cpp):
  implements the standard over-relaxed decomposition (Jasak 1996) of the
  face-normal derivative into a direct two-point term (scaled by
  `1/cosTheta` to correct for the face normal not being parallel to the
  line joining the two cell centroids) plus a correction term (the
  interpolated cell-gradient vector dotted with the non-orthogonal
  remainder direction `k = n - d_hat/cosTheta`). Also returns the plain
  interpolated full gradient vector at the face (distance-weighted average
  of the two cell gradients, `grad_L` on a boundary face since there's no
  second cell) for later reuse by off-diagonal/tangential viscous-stress
  terms, which this decomposition does not itself correct.
- Extended `run_verify_gradient()` to also check every face's corrected
  normal derivative against the exact analytic value, and to report what
  today's `UnstructuredFVMSolver`-style naive two-point difference
  (`(phi_R - phi_L)/dist`, no `cosTheta` correction at all) gives on the
  same faces, for comparison.

**Test performed:** on the same wavy mesh and linear field as Phase 0, using
the Least-Squares cell gradients as input to `face_gradient()`, checked
`dphidn` at every internal and boundary face against the exact
`2.3*nx - 1.1*ny`, alongside the naive uncorrected two-point value.

**Result:** corrected `dphidn` max error ~1.8e-14 (machine precision), naive
two-point max error ~1.74 (against an exact normal derivative whose own
magnitude is ~2.55) -- confirming `UnstructuredFVMSolver::step()`'s current
diffusion flux (which uses the naive form) has a real, large, uncorrected
non-orthogonality error on this mesh, separate from and in addition to the
cell-gradient error already shown in Phase 0. Sanity-checked that this
isn't a degenerate/inverted-cell artifact: the wavy map's Jacobian
determinant is `1 + wobble*pi*sin(pi*(X+Y))`, which for `wobble = 0.15`
stays in `[0.53, 1.47]` everywhere (always positive, no inversions) --
the large naive error is genuine local non-orthogonality, not a broken
test mesh.

**Finding/limitation the static linear-field test cannot demonstrate:** the
plan's original justification for needing the two-part decomposition ("a
naive average of the two cell gradients decouples/checkerboards") is a
statement about *stability under time-marching*, not *accuracy on a static
smooth field*. For a truly linear field the gradient is spatially constant,
so a plain average of two already-accurate (Least-Squares) cell gradients
dotted with the face normal is *also* exact -- there's no way for a
single-snapshot accuracy test to expose the decoupling problem, because
decoupling is about the flux losing direct sensitivity to `(phi_R - phi_L)`
for oscillatory content, not about its accuracy on smooth content. The
corrected scheme is safe from this by construction (its dominant term is
still a direct two-point difference, just correctly scaled), but this
specific claim is asserted from the numerical-analysis argument, not
verified by the test added here. Confirming it directly would need an
oscillatory/checkerboard-mode test or an actual unsteady solve -- worth
keeping in mind if Phase 2/3's solvers ever show odd-even decoupling.
- No near-zero-`cosTheta` (face nearly tangent to the centroid-connecting
  line) guard, same "no mesh-quality safeguards" stance as Phase 0.

---

## Phase 2 — `AdvectionDiffusionFVMSolver`

**Added:**
- [include/AdvectionDiffusionFVMSolver.h](../include/AdvectionDiffusionFVMSolver.h) /
  [src/AdvectionDiffusionFVMSolver.cpp](../src/AdvectionDiffusionFVMSolver.cpp):
  new solver class for `dPhi/dt + div(u*Phi) = alpha*Laplacian(Phi)` with a
  uniform, prescribed advection velocity `(u_adv, v_adv)`. Same two-OpenMP-pass
  flux/residual pattern as `UnstructuredFVMSolver`/`EulerFVMSolver`, plus an
  extra up-front pass reconstructing cell gradients of phi each step (needed
  by every face's diffusive term). Advective flux is first-order upwind;
  diffusive flux uses Phase 1's `face_gradient()`. Reuses
  `BoundaryType`/`BoundaryPatch` as-is. A `GradientCalculator` member has to
  be constructed with an already-geometry-computed mesh, which doesn't fit
  the existing solvers' "compute_geometry() in the constructor body" pattern
  (member initialization order requires it before that point) -- worked
  around with a small `mesh_with_geometry()` helper used in the initializer
  list instead.
- New `EquationSet::AdvectionDiffusion` wired into
  [include/CaseInput.h](../include/CaseInput.h) /
  [src/CaseInput.cpp](../src/CaseInput.cpp) (`equation = advection_diffusion`,
  `u_adv`/`v_adv`/`gradient_scheme` keys; reuses the existing
  alpha/dt/initial_value/initial_radius/boundary keys) and
  `run_advection_diffusion()` in [src/main.cpp](../src/main.cpp) (same
  checkpoint/residual/stopping-criteria contract as `run_diffusion()`). Added
  `CheckpointEquation::AdvectionDiffusion = 2` to
  [include/Checkpoint.h](../include/Checkpoint.h)/[.cpp](../src/Checkpoint.cpp)
  -- no `FORMAT_VERSION` bump needed, since this only adds a new valid tag
  value to an existing 4-byte field, not a layout change.
- `--verify-advdiff` CLI mode (`run_verify_advection_diffusion()` in
  [src/main.cpp](../src/main.cpp)): extended `build_skewed_verification_mesh()`
  (shared with Phase 0/1) to tag its 4 sides as named boundary patches
  ("left"/"right"/"bottom"/"top") so different BCs can be assigned per side --
  each is exactly straight since the wavy perturbation vanishes on the outer
  boundary, so this is additive and doesn't change Phase 0/1's mesh at all.

**Test performed:** steady 1D advection-diffusion of a passive scalar
(`u_adv = (1, 0)`, `alpha = 0.1`, so `Pe = U/alpha = 10`) between Dirichlet
`phi=0` (left, x=0) and `phi=1` (right, x=1), with zero-flux Neumann on
top/bottom, on the same wavy 16x16 mesh construction as Phase 0/1 -- so the
mean flow direction is not aligned with the (warped) interior faces, even
though it's aligned with the mesh's outer boundary. Ran 40,000 explicit
steps at `dt = 0.0005` to reach steady state, then compared cell values
against the exact profile `phi(x) = (exp(Pe*x) - 1) / (exp(Pe) - 1)`.

**Result:** residual norm ~1.8e-15 at the end (confirms steady state was
actually reached, not just "ran out of steps"), max error ~8.2%, L2 error
~3.4% against the exact profile -- a reasonable result for a first-order
upwind advection scheme at this resolution and Peclet number, not a
near-machine-precision match like Phases 0-1 (this is a real PDE
discretization with real truncation error, not a linear-field
reconstruction identity). Also ran a full case-file round trip
(`equation = advection_diffusion` on a small hand-written 2x2-cell
`.fvmesh`) to confirm `CaseInput` parsing, boundary matching, and VTK output
work end-to-end, not just through the standalone verification harness.

**Known limitations / open items:**
- Neumann boundary faces are extrapolated as equal to their owning cell's
  current value for the *gradient reconstruction stencil only* (a
  zero-order approximation, since a Neumann BC prescribes a flux, not a
  value); the boundary flux itself still uses the prescribed flux directly,
  unaffected by this approximation. This happens to be exact for the test
  above (the exact profile doesn't vary with y, so the zero-flux y-direction
  boundary's true value genuinely does equal the interior cell's value), but
  is worth remembering as an approximation, not a general identity, before
  reusing it for a case where the Neumann-boundary field actually varies
  along that boundary.
- `dt` is fixed and unchecked against either the advective CFL limit or the
  explicit-diffusion stability limit, same footgun `UnstructuredFVMSolver`
  already has (documented in `docs/MANUAL.md`'s "Known limitations") --
  now shared by a second solver.
- Advection velocity is uniform and prescribed, not solved for -- there is
  no momentum equation here, by design (that's what Phase 3 adds).
- The 10% L2-error pass threshold in `--verify-advdiff` is a "does this
  discretization behave like first-order upwind should" sanity check, not a
  tight correctness bound the way Phases 0-1's 1e-9 thresholds are; a real
  regression (e.g. an accidentally-swapped sign) would very likely still be
  caught since it would push the error well past 10%, but a subtler bug
  producing, say, 12% error instead of 8% would not necessarily be caught by
  this test alone.

---

## Phase 3 — `NavierStokesFVMSolver`

**Added:**
- [include/NavierStokesFVMSolver.h](../include/NavierStokesFVMSolver.h) /
  [src/NavierStokesFVMSolver.cpp](../src/NavierStokesFVMSolver.cpp): new
  solver class (per the 2026-07-05 architecture decision -- a second class,
  not a `mu > 0` toggle inside `EulerFVMSolver`) for the compressible
  Navier-Stokes equations: `EulerFVMSolver`'s exact inviscid flux machinery
  (reuses `NumericalFluxScheme`/`rusanov_flux()`/`hllc_flux()`/
  `exact_riemann_flux()` unchanged) plus a Newtonian viscous stress tensor
  and Fourier heat conduction, assembled from corrected face gradients of
  u, v, T (three separate `GradientCalculator::compute()` calls per step,
  sharing one `GradientCalculator` instance since its precomputed
  Least-Squares matrix depends only on mesh geometry, not on which field is
  being differentiated).
  - New `NSBoundaryType`/`NSBoundaryCondition` (`NoSlipWall`/`Farfield`/
    `Outflow`) -- a no-slip wall mirrors BOTH velocity components in its
    ghost state (not just the normal one, unlike Euler's slip wall) and
    carries a thermal condition (isothermal via `wall_temperature`, or
    adiabatic = zero heat flux prescribed directly, bypassing gradient
    reconstruction for that term the same way a Neumann BC does elsewhere).
  - `corrected_face_gradient_vector()` (anonymous-namespace helper in the
    `.cpp`): since `face_gradient()` only corrects a face's NORMAL
    derivative (Phase 1), the stress tensor's tangential/off-diagonal terms
    need the full 2D gradient vector -- built by taking the interpolated
    average gradient and replacing just its normal component with the
    corrected one (`grad_f + (dphidn - grad_f.n)*n`). This is the concrete
    per-component realization of the "two-part face gradient" the plan
    called for, now applied to a vector field instead of one scalar.
  - `compute_dt()` extends `EulerFVMSolver`'s inviscid-CFL formula with the
    standard (Blazek-style) explicit-viscous-diffusion term
    `2*nu/length` added to each face's wave-speed estimate, `nu = mu/rho`.
  - Viscous flux is exactly zero at Farfield/Outflow boundaries by
    construction (not approximated/extrapolated) -- those boundaries are
    assumed placed far enough from any wall that this doesn't matter.
- `EulerState.h` gained `temperature()` (ideal-gas `p = rho*R*T`), the first
  place in this codebase to need a gas constant -- Euler itself never did
  (pressure/sound speed are R-independent).
- New `EquationSet::NavierStokes` wired fully into
  [include/CaseInput.h](../include/CaseInput.h)/[.cpp](../src/CaseInput.cpp)
  (`equation = navier_stokes`, `mu`/`prandtl`/`gas_constant` keys, `ns_init`
  mirroring `euler_init`'s grammar in its own storage, and
  `ns_wall`/`ns_wall_isothermal <T>`/`ns_farfield <rho u v p>`/`ns_outflow`
  boundary keywords -- prefixed with `ns_` specifically because "farfield"/
  "outflow" already mean something for Euler and need to land in a
  different spec vector) and `run_navier_stokes()`/`write_navier_stokes_fields()`
  in [src/main.cpp](../src/main.cpp) (same checkpoint/residual/stopping-criteria
  contract as `run_euler()`; VTK output adds a `T` field). Added
  `CheckpointEquation::NavierStokes = 3` to
  [include/Checkpoint.h](../include/Checkpoint.h) -- again no `FORMAT_VERSION`
  bump, same reasoning as Phase 2's tag addition.
- `--verify-ns-uniform` CLI mode (`run_verify_navier_stokes_uniform()` in
  [src/main.cpp](../src/main.cpp)).

**Test performed:** deliberately NOT the Couette validation (that's Phase
4's job). Phase 3's own gate: a uniform freestream state prescribed
identically as the `Farfield` BC on every side of the same wavy
non-orthogonal mesh used throughout this tracker, with `mu = 0.5` (viscosity
genuinely on), run for 20 steps. A spatially uniform state has zero gradient
everywhere, so both the (already-established) inviscid flux and the new
viscous stress/heat-flux terms must come out to exactly zero at every face
regardless of mesh non-orthogonality -- this is the direct NS analogue of
Phases 0-1's "reconstruct a known field exactly" tests, now aimed at the new
viscous flux assembly instead of a standalone scalar gradient. Also ran a
full case-file round trip (`equation = navier_stokes` on the same small
hand-written 2x2-cell `.fvmesh` from Phase 2, with a farfield inflow, a
no-slip adiabatic top/bottom, and an outflow) to confirm `CaseInput` parsing
(including the new `ns_*` boundary keywords), boundary matching, and VTK
output work end-to-end.

**Result:** max deviation from the initial uniform state after 20 steps was
exactly `0` (not just "small" -- bit-for-bit identical), confirming the
viscous flux assembly injects no spurious flux from mesh skewness when there
is genuinely nothing to diffuse. The case-file smoke test produced sane,
symmetric (top/bottom, as expected from the symmetric mesh/BCs), finite
`rho`/`u`/`v`/`p`/`T` fields with no NaNs.

**Known limitations / open items:**
- **Not yet validated against real viscous physics.** The uniform-flow test
  proves the viscous terms are silent when there's nothing to diffuse; it
  says nothing about whether the shear stress or heat conduction magnitudes
  are quantitatively correct for a real velocity/temperature gradient. That
  is explicitly Phase 4's job (Couette flow, with a known analytic linear
  velocity profile) -- do not treat Phase 3 as validating the physics, only
  the plumbing.
- Viscous flux is forced to exactly zero at Farfield/Outflow boundaries
  rather than extrapolated -- correct for the intended use (open boundaries
  placed away from walls) but would silently under-predict viscous effects
  if a case ever puts a Farfield/Outflow boundary close to a region with a
  real velocity/temperature gradient.
- `compute_dt()`'s viscous stability term (`2*nu/length`) is the standard
  literature formula but hasn't been stress-tested against an actual
  stiff/high-viscosity case yet -- only exercised at `mu = 0.5` /
  `Re ~ O(1)` in the uniform-flow test (where it doesn't matter anyway,
  since nothing changes) and `mu = 0.01` in the smoke test.
- No mesh-quality safeguards, same "by design gap" stance as every prior
  phase.
- `gas_constant`/`prandtl` default to 1.0/0.72 but there is no validation
  that `gas_constant > 0` or `prandtl > 0` -- a zero or negative value
  would silently produce `Inf`/division-by-zero rather than a diagnostic.

---

## Phase 4 — Couette flow validation

**Added:**
- `NSBoundaryCondition`/`NSBoundaryConditionSpec` gained `wall_u`/`wall_v`
  (default 0, 0 -- fully backward compatible with Phase 3's stationary-only
  walls). `ghost_state()`'s `NoSlipWall` branch now mirrors velocity about
  `(wall_u, wall_v)` instead of about zero
  (`u_ghost = 2*wall_u - u_L`), and `build_boundary_fields()` feeds
  `wall_u`/`wall_v` into the gradient stencil instead of a hardcoded 0.
  New case-file boundary keywords `ns_wall_moving <u> <v>` and
  `ns_wall_moving_isothermal <u> <v> <T>`.
- `build_structured_verification_mesh()`: the shared topology-building code
  behind the test-mesh generators was refactored to take a `displace`
  functor and separate `n_x`/`n_y` cell counts (previously hardcoded
  square), so `build_skewed_verification_mesh()` (Phases 0-2, unchanged
  behavior) and a new `build_sheared_verification_mesh()` (this phase) can
  share it instead of duplicating ~90 lines of face/cell-topology code.
- `build_sheared_verification_mesh()`: a pure shear (`x = X + shear*Y,
  y = Y`) rather than the wavy perturbation -- see the "what went wrong"
  section below for why this phase needed a different mesh than Phases 0-2.
- `--verify-couette` CLI mode (`run_verify_couette()`), using its own
  step-by-step loop with an explicit NaN/Inf divergence check (the earlier
  verification harnesses all called `solver.run(n)` blindly; this phase's
  debugging surfaced that gap -- worth carrying the same divergence guard
  back into `run_verify_navier_stokes_uniform()`/`run_verify_advection_diffusion()`
  if either is ever pushed to a less forgiving parameter regime).

**What went wrong (kept in detail because both failures are genuine,
transferable findings, not throwaway debugging noise):**

1. First attempt reused Phases 0-2's wavy mesh. Result: steady-state L2
   error ~24-50% (worse with more distortion), and it did NOT shrink with
   5x more time steps -- a genuine (wrong) steady state, not a slow one.
   Plotting `u` vs `y` at fixed `x` showed the deviation from the exact
   linear profile was shaped exactly like a parasitic Poiseuille (parabolic,
   zero at both walls, peaked at mid-height) component added on top of the
   real profile. Root cause: the wavy mesh's own distortion amplitude
   varies with `x` (peaks at `x = 0.5`), so the *discretization* is not
   x-translation-invariant even though the continuous Couette solution is
   -- that mismatch injects a spurious x-varying truncation-error "forcing."
   This hadn't mattered for Phases 0-2 because neither a static gradient
   check nor a flow-aligned-with-x advection-diffusion check depends on
   x-invariance the way Couette flow's validity as a test case does.
2. Switched to a pure shear (x-translation-invariant by construction --
   every row is an identical parallelogram just translated in x). Error
   dropped a lot (~6-9% depending on shear magnitude) but still didn't pass
   5%, residual plateaued around `1e-5`-`1e-6` instead of dropping to
   machine zero, and -- critically -- running 200,000 steps (10x) at
   `cfl = 0.3` actually **diverged** (NaN) partway through, a real stability
   finding caught only because the harness was rewritten to check for
   NaN/Inf every step instead of calling `run()` blindly. Dropping to
   `cfl = 0.1` fixed the divergence but not the accuracy: plotting `u`/`v`
   across a fixed-height row (varying x this time) showed both varying
   *smoothly and monotonically across the entire domain width* -- not
   localized to the corners. Raising the Mach number (to test a low-Mach
   acoustic-stiffness hypothesis) made no difference, which ruled that
   hypothesis out and pointed at the real cause: `NSBoundaryType::Outflow`
   is a zero-gradient/do-nothing condition, not a true periodic one. Nothing
   in it forces the two ends of an x-run to match, so a slow, smoothly
   linear-in-x drift mode has no restoring force pushing it back to the
   (also valid) x-invariant solution -- it's a real, persistent artifact of
   approximating "periodic" with "outflow," not a bug in the viscous terms
   themselves.
3. Fix: make the domain exactly one cell wide in x (`n_x = 1`,
   `n_y = 16`). This doesn't shrink the drift mode, it removes it
   structurally -- there is no room for `u`/`v` to vary with x at all when
   there's only one cell per row, while the shear still makes every
   internal (horizontal) face non-orthogonal in Jasak's sense. Result: `rho_u`
   residual `~1e-16`, max/L2 error against the exact profile `~1e-15`
   relative to `U` -- machine precision. Confirmed genuinely stable (not
   luck) by rerunning at 25x the step count (500,000 steps): residual and
   error were statistically unchanged (`9.3e-17`, unaffected).

**Test performed (final form):** a 1x16 sheared mesh (`shear = 0.3`),
stationary bottom wall, top wall moving at `U = 0.1` (`ns_wall_moving`),
zero-gradient outflow left/right, `mu = 0.02`, run 20,000 explicit steps,
compared against `u(y) = U*y/H`. Also ran a full case-file round trip
exercising the new `ns_wall_moving` keyword on the Phase 2/3 hand-written
2x2-cell `.fvmesh`, producing a sane transient (top row visibly approaching
the prescribed wall speed faster than the bottom row after only 500 steps,
left/right symmetric as expected from the symmetric outflow BCs).

**Result:** machine-precision match (`rho_u` residual `~1e-16`, L2 error
`~4e-15` relative to `U`) to the analytic linear Couette profile, confirmed
stable over 25x more steps than needed. This is the real validation Phase 3
deferred: the corrected viscous shear stress (Phase 1's face gradient,
generalized to a full tensor in Phase 3) reproduces physically correct,
quantitatively exact shear transport across a genuinely non-orthogonal
mesh, not just "stays silent when there's nothing to diffuse."

**Known limitations / open items:**
- **`NSBoundaryType::Outflow` should not be used on both ends of a
  domain that's more than one cell wide in the flow direction, for any
  problem relying on x-invariance for correctness** -- it will not enforce
  that invariance, and a slow, non-decaying drift mode can develop and
  persist indefinitely (not just "take a long time to settle"). A real
  periodic BC is the correct fix if this solver ever needs a
  multi-cell-wide x-invariant channel case; short of that, keep such cases
  one cell wide in the invariant direction, as this test now does.
- This validation is narrow by construction: no pressure gradient, no
  compressibility effects worth mentioning (`Mach ~ 0.08`), adiabatic walls
  (no viscous-heating/temperature check, as originally agreed), and a
  single cell of x-extent. It confirms the viscous shear term is
  quantitatively correct; it says nothing about pressure-gradient-driven
  (Poiseuille) flow, compressible/shock-viscous interactions, or the
  heat-conduction term's quantitative accuracy (only its plumbing was
  checked, in Phase 3's adiabatic-wall path).
- **Follow-up (same day) on the `cfl = 0.3` divergence noted above:** initially
  suspected (and documented here, briefly) as a general viscous-CFL margin
  deficiency on skewed meshes, requiring `compute_dt()`'s
  `length = volume/face.area` proxy to be hardened. Two targeted diagnostics
  disproved that before any such fix was made: (1) the SAME 16x16 sheared
  mesh/mu/step-count at `cfl = 0.3` but with `shear = 0` (orthogonal) ran the
  full 200,000 steps with no divergence (residual settled at `~1e-7`); (2)
  the already-shipped `n_x = 1` geometry (`shear = 0.3`, same mu) ran the
  full 200,000 steps at `cfl = 0.3` -- the exact value that diverged before
  -- with **no divergence and machine-precision accuracy** (`rho_u` residual
  `3.8e-17`). Together these show the divergence was entirely a symptom of
  the multi-column `Outflow`-drift-mode instability described above (item 2
  under "What went wrong"), not a general weakness in the viscous CFL
  formula. `compute_dt()` was left unchanged; `cfl = 0.3` is confirmed safe
  and restored as this test's value. Lesson for next time: attribute an
  instability to the mechanism you can actually isolate with a targeted
  diagnostic, not the mechanism you happen to be looking at when it appears
  -- the first plausible-sounding explanation (skewed-mesh CFL margin) was
  wrong, and cheap A/B tests (change one variable, rerun) caught it before
  an unneeded code change shipped.

---

## Phase 5 — DNS exposure (Tier 1 diagnostic + Tier 2 deferred plan)

Following up the Couette validation, Mathieu asked whether the code could
expose a DNS option with relevant user parameters. Key clarification before
any work started: DNS is not a solver *mode* to switch on -- it's a property
of resolution relative to the flow. `NavierStokesFVMSolver` already solves
the raw, unmodeled compressible NS equations (no RANS/LES closure), which is
the DNS-appropriate equation set by definition. The real question was
whether the code helps a user run one *correctly*, and a critical scoping
point that changed the shape of the answer: this solver is 2D-only, by
deliberate, non-negotiable architectural design (see MANUAL.md's "Known
limitations"), and classical Kolmogorov scaling / the turbulent energy
cascade are 3D results (vortex stretching doesn't exist the same way in
2D). So "DNS" here can only ever mean fully-resolved 2D unsteady, never 3D
turbulence -- that framing is threaded through every piece of Phase 5.

Split into two tiers:

**Tier 1 (implemented) -- resolution diagnostic, no new physics:**
- `NavierStokesFVMSolver::compute_resolution_diagnostics()`
  ([NavierStokesFVMSolver.h](../include/NavierStokesFVMSolver.h)/[.cpp](../src/NavierStokesFVMSolver.cpp)):
  a posteriori estimate from the CURRENT flow field -- per cell, strain rate
  `S_ij` from the existing `grad_u`/`grad_v` (no new gradient machinery),
  dissipation `epsilon = 2*nu*(S_11^2+S_22^2+2*S_12^2)`, Kolmogorov length
  `eta = (nu^3/epsilon)^0.25`, ratio `h/eta` (`h = sqrt(volume)`). Cells with
  `epsilon` below a small floor (no local shear -- e.g. `mu = 0`, or before
  any flow develops) are excluded from the aggregate rather than producing
  `eta -> infinity`; returns domain min/max/mean plus `n_active` (how many
  cells had something to report).
- New `resolution_report_file`/`resolution_report_interval` case-file keys
  (Navier-Stokes only), wired into `run_navier_stokes()` exactly like
  `residual_file`, including the `scratch_dir` rebasing the other three
  output paths already got (a real gap caught and fixed while wiring this
  in -- `resolution_report_file` had been missing from `main()`'s
  `resolve_output_path` calls).
- **Test performed:** smoke test on the existing moving-wall case
  (Phase 3/4's hand-written 2x2-cell `.fvmesh`, `mu = 0.05`, moving top
  wall) with `resolution_report_interval = 50`. Result: all 4 cells
  reported as active by step 50 (shear has diffused through the whole
  domain), ratios finite and smoothly evolving (~1.40-1.41 across 500
  steps), no NaN. Full regression pass on all four `--verify-*` flags
  confirmed unaffected.
- **Known limitation:** the classical Kolmogorov scaling underlying `eta`'s
  formula is a 3D turbulence result; using it here is a pragmatic resolution
  heuristic for 2D unsteady motion, not a validated 2D theory. Don't cite
  the `h/eta` ratio as if it carries the same meaning it would in a real 3D
  DNS quality check.

**Tier 2 (planning only, not implemented) -- higher-order time/space
discretization:** fully written up in
[docs/dns-higher-order-scheme-plan.md](dns-higher-order-scheme-plan.md)
(RK3/RK4 time integration vs. MUSCL spatial reconstruction, with a
recommendation) -- not restated here, to avoid the two docs drifting apart;
read that doc for the detail. Neither option is implemented.

**Consolidated documentation pass** (also part of Phase 5, since MANUAL.md
and CLAUDE.md had never been touched across Phases 0-4 despite four new
solvers, a new shared gradient module, and a dozen new case-file keys
landing): MANUAL.md gained a new "Gradient reconstruction" section, full
case-file/boundary-condition/output-format coverage for
`advection_diffusion` and `navier_stokes`, and expanded "Known limitations"
(2D-vs-3D DNS framing, no turbulence closure, `Outflow` non-periodicity,
unvalidated `gas_constant`/`prandtl`, viscous flux forced to zero at
open boundaries). CLAUDE.md gained matching architecture bullets sized to
its existing terse style (pointers + the one non-obvious decision per
topic, not MANUAL.md's depth).

**Known limitations / open items:**
- Tier 1's diagnostic has no stopping-criterion integration (there is no
  `resolution_tolerance` key) -- it is purely informational, by design; see
  [Monitoring & stopping criteria](../docs/MANUAL.md#monitoring--stopping-criteria).
- Tier 2 remains unimplemented; `NavierStokesFVMSolver` (and
  `EulerFVMSolver`, which it reuses unchanged) stay first-order in time and
  space until/unless that plan is picked up.
- The 2D-vs-3D DNS distinction is documented in three places now (this
  tracker, MANUAL.md, CLAUDE.md) specifically so it can't be quietly lost
  the way `resolution_report_file`'s `scratch_dir` rebasing almost was --
  but it still depends on whoever extends this code next actually reading
  one of them before implying broader capability than the 2D architecture
  supports.

---

## Phase 6 — `compute_dt()` viscous length-scale hardening

Triggered by Bravo integration testing hitting the first real (non-synthetic)
Navier-Stokes divergence: a laminar run on a real boundary-layer-clustered
airfoil mesh from Bravo's `AirfoilMesherRunner` (`blFirstLayerHeight =
0.001`, `blGrowthRatio = 1.05`, Mach 0.98, `mu = 0.02`, `cfl = 0.5`) diverged
to NaN by step 25. Full write-up of the original finding, bisection, and
Bravo's own client-side mitigation (lowering its default `cfl` 0.5 -> 0.3) is
in [docs/ns-cfl-margin-and-farfield-bc-findings.md](ns-cfl-margin-and-farfield-bc-findings.md)
-- not restated here, to avoid the two docs drifting apart.

**Not a contradiction of Phase 4's own prior finding.** Phase 4's "Known
limitations" section above already investigated `compute_dt()`'s `length =
volume/face.area` proxy once, on a *different* mesh and mechanism: a
`cfl = 0.3` divergence on a 16x16 sheared (skewed but NOT wall-normal-
stretched) mesh, which two targeted A/B tests traced entirely to the
`Outflow`-BC drift-mode instability, not the viscous length proxy --
`compute_dt()` was explicitly left unchanged at the time, and that
conclusion still stands for that case. This phase's finding is a real,
different deficiency in the same proxy, on a genuinely anisotropic
boundary-layer mesh (a short wall-normal cell dimension against a much
longer streamwise one), not exercised by any prior phase's mesh -- caught by
a real Bravo case, not this tracker's own validation suite.

**Root cause:** `volume/face.area` is a single scalar averaged over both of
a cell's dimensions. On an anisotropic boundary-layer cell this
systematically biases toward the longer (streamwise) dimension rather than
the short wall-normal one where the viscous term's real stiffness lives, so
`compute_dt()`'s viscous stability estimate (`2*nu/length`) understates how
restrictive the true limit should be. Separately (found by direct
inspection, not by this divergence): at a domain boundary, the old formula
used the wall cell's *full* height (`volume/area`) rather than the
*centroid-to-wall* distance -- exactly half that, for a rectangular cell --
a real 2x error at every boundary face independent of any skew.

**Fix:** replaced `length` in both `NavierStokesFVMSolver::compute_dt()`
and `RANSFVMSolver::compute_dt()` (which had copied the identical pattern
for its combined laminar+turbulent viscous term -- RANS's whole use case is
wall-resolved anisotropic meshes, so it inherited the same exposure) with a
new shared free function, `face_normal_distance()`
([GradientReconstruction.h](../include/GradientReconstruction.h)/[.cpp](../src/GradientReconstruction.cpp)).
It returns the face-normal-projected distance between `face.cell_left`'s
centroid and the opposite point (`face.cell_right`'s centroid, or the face
midpoint for a boundary face) -- exactly the `dist * cosTheta` quantity
`face_gradient()`'s own over-relaxed decomposition already computed inline
for the viscous/diffusive flux itself. This makes the stability estimate
and the flux calculation consistent (previously they used two unrelated
length scales), which is also the physically correct choice: the explicit
diffusion stability limit `dt <= length^2/(2*nu)` needs `length` to be the
actual spacing the gradient is taken over. `face_gradient()` was refactored
to call the same new helper for its own `dist`/`cosTheta` rather than
duplicating the formula. `EulerFVMSolver::compute_dt()` was left unchanged
-- no viscous term, not implicated.

A single `length` value is used for both the numerator (advective CFL) and
the viscous denominator, not two different length scales for each --
mixing the old, larger `volume/area` estimate into one term while using the
new, smaller face-normal distance in the other would under-conservatively
combine the two and could leave a residual version of the same bug.

**Test performed:** all eight pre-existing `--verify-*` gates re-run and
confirmed unaffected (`--verify-gradient`, `--verify-advdiff`,
`--verify-ns-uniform`, `--verify-couette`, `--verify-wall-distance`,
`--verify-sa-source`, `--verify-rans-stability`, `--verify-flat-plate`) --
none of their meshes are anisotropic enough for the fix to change their
result. A new `--verify-ns-stretched-cfl` gate (`run_verify_ns_stretched_cfl()`
in [src/main.cpp](../src/main.cpp)) was added specifically to exercise the
fix: a laminar `NavierStokesFVMSolver` run on `build_flat_plate_mesh()`
(RANS Phase 4's wall-normal-stretched mesh generator -- the only existing
mesh builder with a genuinely thin first cell), `first_cell_height = 1e-4`,
`growth_ratio = 1.05`, `mu = 0.02`. Confirmed via direct A/B (not just "add
the fix and hope," per Phase 4's own lesson above): at `cfl = 2.0` on this
mesh, the pre-fix formula diverges to NaN by step 4; the post-fix formula
completes the full 2000 steps with no divergence. (This synthetic mesh is
purely orthogonal, unlike Bravo's curved-airfoil mesh, so its cliff sits at
a much higher `cfl` than Bravo's real case -- consistent with the boundary-
face 2x error alone being enough to matter here, without needing skew on
top of it.)

**Result:** PASS on all nine gates. The A/B test above is the one that
actually demonstrates the fix; the other eight are non-regression evidence
that isotropic/mildly-skewed meshes are unaffected.

**Known limitations / open items:**
- Not validated against Bravo's actual airfoil mesh/case file directly (that
  mesh isn't part of this repo) -- only against a synthetic mesh built to
  exercise the same mechanism (thin, anisotropic boundary-layer cells). The
  underlying formula change is the same either way, but the specific
  `cfl = 0.4`/`0.5` boundary Bravo observed on its real mesh was not
  independently re-measured after this fix.
- This fix makes `dt` smaller (more steps for the same physical time) on
  any strongly anisotropic mesh, not just the one that triggered this
  investigation -- an intentional cost of correctness, not a regression, but
  worth knowing before assuming a boundary-layer-resolving Navier-Stokes run
  will run at the same wall-clock cost as before.
- The `Farfield` ghost-state limitation flagged alongside the original
  finding (no characteristic branching by local Mach number -- see
  [docs/ns-cfl-margin-and-farfield-bc-findings.md](ns-cfl-margin-and-farfield-bc-findings.md))
  is unrelated to this fix and remains open; no decision has been made on
  it.
- No mesh-quality safeguards, same "by design gap" stance as every prior
  phase in both trackers -- in particular, `face_normal_distance()`
  inherits `face_gradient()`'s existing unguarded assumption that
  `cosTheta` isn't (near-)zero.
