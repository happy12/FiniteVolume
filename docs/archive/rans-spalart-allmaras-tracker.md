# RANS (Spalart-Allmaras) implementation tracker

**Archived (2026-07-05): all 4 planned phases closed** (Phase 4 landed on a
relaxed laminar-boundary-layer validation target rather than the originally
planned turbulent log-law match -- see that phase's notes below for why).
Kept for its debugging history and architecture rationale; see
[CLAUDE.md](../../CLAUDE.md) and [docs/MANUAL.md](../MANUAL.md) for the
current, load-bearing summary of what `RANSFVMSolver` actually does.

Tracks progress toward a RANS turbulence-modeling capability built on the
existing `NavierStokesFVMSolver` infrastructure, per the plan agreed with
Mathieu (2026-07-05, immediately following
[docs/navier-stokes-tracker.md](../navier-stokes-tracker.md)'s Phase 5). Same
discipline as that tracker: each phase gets its own verification gate
before the next one starts, and this file records what was added, what was
actually tested, and what's still a known gap or open question when a
phase closes.

## Architecture decisions locked in (confirmed with Mathieu, 2026-07-05)

- **New `RANSFVMSolver` class**, not a `turbulence_model` parameter added
  to `NavierStokesFVMSolver`. Consistent with this project's existing "one
  class per equation set" pattern (Euler vs. Navier-Stokes are separate
  classes even though Navier-Stokes reuses Euler's inviscid flux
  unchanged) -- `RANSFVMSolver` duplicates `NavierStokesFVMSolver`'s
  inviscid/viscous flux, ghost-state, and `compute_dt()` structure, and
  adds the turbulence transport equation, wall-distance field, and
  `mu_eff` substitution on top.
- **SA-noft2 variant**: no trip term (`ft2`), fully-turbulent-flow
  assumption everywhere. Standard modern default; one fewer function
  (and one fewer set of constants, `Ct1`-`Ct4`) to implement and verify.
  Model constants: `Cb1 = 0.1355`, `Cb2 = 0.622`, `sigma = 2/3`,
  `kappa = 0.41`, `Cw2 = 0.3`, `Cw3 = 2.0`, `Cv1 = 7.1`,
  `Cw1 = Cb1/kappa^2 + (1+Cb2)/sigma` (derived, not independent). Later
  expanded in Phase 3, re-confirmed with Mathieu at the time rather than
  assumed, to add the standard negative-`S~` robustness fix's `Cv2 = 0.7`,
  `Cv3 = 0.9` -- see that phase's notes for why.
- **Full validation committed**: Phases 1-3 (infrastructure + stable
  coupling, sanity-checked but not yet quantitatively accurate) are
  explicitly *not* the finish line -- Phase 4 (flat-plate turbulent
  boundary layer vs. the log-law velocity profile) is in scope for this
  initiative from the start, not deferred as a follow-up. Expect this to
  need real iteration, comparable to (plausibly harder than) the Couette
  validation's own debugging journey in
  [docs/navier-stokes-tracker.md](../navier-stokes-tracker.md) Phase 4 --
  SA is known in the literature to be sensitive to the freestream `nut`
  value, and a wall-resolved (`y+ ~ 1`) mesh plus a proper turbulent inflow
  condition are both harder setup problems than Couette's fully-developed,
  x-invariant flow was.

## What SA adds on top of `NavierStokesFVMSolver`

One extra transported scalar, `nut` ("nu-tilde"), plus a derived eddy
viscosity `nu_t = nut * fv1`. It couples into the existing viscous flux in
exactly one place: molecular `mu` is replaced by
`mu_eff = mu + rho*nu_t` in the stress tensor, and thermal conductivity
gets a turbulent contribution via a turbulent Prandtl number
`Pr_t ~ 0.9`. Everything else about the Navier-Stokes flux assembly
(inviscid flux, ghost states, `face_gradient()`-based viscous terms) is
unchanged.

The `nut` transport equation is structurally an advection-diffusion
equation (like [docs/navier-stokes-tracker.md](../navier-stokes-tracker.md)
Phase 2's `AdvectionDiffusionFVMSolver`, but advected by the REAL flow
velocity from the coupled Navier-Stokes solution, not a prescribed one)
plus two new terms:
- **Production** `P = Cb1*S_tilde*nut`, where `S_tilde` is a modified
  vorticity magnitude -- computable directly from the already-existing
  `grad_u`/`grad_v` cell gradients (`Omega = |dv/dx - du/dy|` in 2D, plus
  SA's near-wall correction using `fv2`). No new gradient infrastructure
  needed.
- **Destruction** `D = Cw1*fw*(nut/d)^2`, where `d` is distance to the
  nearest wall -- a genuinely new piece of infrastructure; nothing in this
  codebase currently has a wall-distance concept at all.

## Phase status

| # | Phase | Status |
|---|---|---|
| 1 | Wall-distance module + verification | Done |
| 2 | SA source terms in isolation + sanity check | Done |
| 3 | Coupling into `RANSFVMSolver`'s viscous flux + turbulence-aware `compute_dt()` | Done |
| 4 | Flat-plate turbulent boundary layer vs. log-law validation | Attempted -- target relaxed, see below |

## Phase 1 — Wall-distance module

**Added:**
- [include/WallDistance.h](../../include/WallDistance.h) /
  [src/WallDistance.cpp](../../src/WallDistance.cpp): `compute_wall_distance()`,
  a free function computing, per cell, the minimum Euclidean point-to-segment
  distance from that cell's centroid to the nearest face in a caller-supplied
  list of wall face indices. Brute-force `O(cells * wall_faces)`, no
  acceleration structure -- consistent with the project's existing "simple
  first, document scaling limits" stance (see CLAUDE.md's
  mesh-quality-safeguards precedent). Deliberately has no dependency on
  `RANSFVMSolver`/`NSBoundaryType` (neither exists yet) -- it depends only on
  mesh geometry and a list of face indices, the same reasoning as
  `GradientCalculator`'s precomputed Least-Squares matrix
  ([GradientReconstruction.h](../../include/GradientReconstruction.h)). A future
  `RANSFVMSolver` just needs to hand it the indices of whichever faces have
  `NSBoundaryType::NoSlipWall`, once wired in at Phase 3. Returns
  `+infinity` per cell if the wall face list is empty, rather than `0.0` --
  there being no wall to be far from is not the same as being distance zero
  from one.
- `--verify-wall-distance` CLI mode (`run_verify_wall_distance()` in
  [src/main.cpp](../../src/main.cpp)): reuses
  [docs/navier-stokes-tracker.md](../navier-stokes-tracker.md) Phase 4's
  `build_sheared_verification_mesh()` (1x16 sheared mesh, `shear = 0.3`)
  purely for its geometry -- no solver involved, since the module doesn't
  need one.

**Test performed:** exactly the planned verification -- the "bottom" patch
of the 1x16 sheared mesh is a straight horizontal line at `y = 0` for every
`x` (the shear `x = X + shear*Y, y = Y` leaves `Y = 0`'s row exactly at
`y = 0` regardless of shear), so the exact distance from any cell centroid
to it is simply that centroid's `y_centroid`. Computed `compute_wall_distance()`
against only that one wall face and compared every cell's result to its
exact `y_centroid`.

**Result:** max error `0` (bit-for-bit exact, not just "small") -- the
perpendicular-distance case has no floating-point subtraction of nearly
equal quantities to lose precision on, unlike e.g. Phase 0/1 of the NS
tracker's gradient checks (`~1e-14`/`~1e-15`). Full regression pass on all
four pre-existing `--verify-*` flags (`--verify-gradient`, `--verify-advdiff`,
`--verify-ns-uniform`, `--verify-couette`) confirmed unaffected.

**Known limitations / open items:**
- Only ever tested against a single flat wall face, per the tracker's own
  planned scope -- not yet exercised with multiple wall patches (e.g. both
  Couette walls) or a curved/piecewise wall boundary. The point-to-segment
  distance math has no reason to behave differently in that case, but it
  hasn't been checked with a test that would catch a wrong-wall-picked bug
  if one existed.
- No spatial acceleration structure, so this does not scale past
  `O(cells * wall_faces)` -- fine at this project's current mesh scale, would
  need revisiting (e.g. a tree over wall face midpoints) if a much larger
  wall-face count is ever used.
- Not yet wired into any solver or case file -- `RANSFVMSolver` doesn't exist
  yet (Phase 3's job); this phase is geometry-only infrastructure, same
  relationship Phase 0 of the NS tracker had to Phase 2/3's solvers.

## Phase 2 — SA source terms in isolation

**Added:**
- [include/SpalartAllmaras.h](../../include/SpalartAllmaras.h) /
  [src/SpalartAllmaras.cpp](../../src/SpalartAllmaras.cpp): the SA-noft2 model
  constants (`SA_CB1`, `SA_CB2`, `SA_SIGMA`, `SA_KAPPA`, `SA_CW2`, `SA_CW3`,
  `SA_CV1`, derived `SA_CW1`, exactly as fixed in this tracker's architecture
  decision), `sa_fv1()`/`sa_eddy_viscosity()`, and `compute_sa_source_terms()`
  -- a pure per-cell function returning `SASourceTerms{production,
  destruction, cross_diffusion}` from local `nut`/`nu`/vorticity/wall-distance/
  `grad(nut)` values. Deliberately does NOT assemble the transport equation's
  conservative diffusion term `div((nu+nu_t)/sigma * grad(nut))` -- that needs
  mesh face geometry and reuses `GradientReconstruction.h`'s already-verified
  `face_gradient()`/`GradientCalculator` unchanged, deferred to Phase 3's
  actual solver coupling, the same way `NavierStokesFVMSolver` reuses them for
  the viscous stress tensor. `destruction` is forced to exactly `0.0` when
  `nut <= 0` rather than evaluated through the `r`/`g`/`fw` chain -- a real bug
  avoided before it shipped: `r`'s formula divides by `S~ * kappa^2 * d^2`,
  which is a literal `0/0` when `nut == 0` AND `omega == 0` (`S~ == 0` too in
  that state), and `nut * NaN` is `NaN`, not the physically-correct `0`, so
  the `nut <= 0` guard is load-bearing, not defensive boilerplate.
- `--verify-sa-source` CLI mode (`run_verify_sa_source_terms()` in
  [src/main.cpp](../../src/main.cpp)): reuses `build_sheared_verification_mesh()`
  (Phase 4 of the NS tracker) and `compute_wall_distance()` (Phase 1 of this
  tracker) for geometry, and `GradientCalculator` (Phase 0 of the NS tracker)
  to reconstruct `grad(nut)` for Test 2 -- no new mesh infrastructure needed.

**Test performed:** exactly the two planned checks.
1. Zero vorticity, `nut` initialized to 0 everywhere on the 1x16 sheared
   mesh: confirmed every cell's production/destruction/cross_diffusion is
   exactly `0.0` (not just small), then ran a 1000-step explicit-Euler loop
   of `nut` using those source terms and confirmed `nut` stayed at exactly
   `0.0` in every cell throughout -- a genuine fixed-point check, not a
   single-evaluation snapshot.
2. A manufactured, smooth, always-positive `omega(y) = 1 + y` and
   `nut(y) = 0.01*(1+y)` field (not physically self-consistent, per the
   tracker's own wording) on the same mesh, chosen so `S~` stays comfortably
   positive everywhere (deliberately avoiding the known negative-`S~` gap
   noted below): checked every cell's three source terms for finiteness and
   non-negativity.

**Result:** both PASS. Test 1: all source terms exactly `0.0` every step,
`nut` exactly `0.0` after 1000 steps. Test 2: max production `0.0439688`, max
destruction `0.707294`, max cross-diffusion `9.33e-05` across the mesh, all
finite and non-negative -- confirming the `fv1`/`fv2`/`r`/`g`/`fw` nonlinear
closure chain doesn't blow up or invert sign on a smooth field, including
near the wall where the `1/d^2` term in `S~` gets large (`d` as small as
`~0.031` on this mesh) but the standard `r` clip (`r <= 10`) and `fw`'s own
boundedness kept every value sane without needing any extra guard beyond
the `nut <= 0` one already in `compute_sa_source_terms()`. Full regression
pass on all five pre-existing `--verify-*` flags confirmed unaffected.

**Known limitations / open items:**
- **No negative-`S~` robustness fix at the time this phase closed** (the
  `Cv2`/`Cv3`-based modification some SA implementations add when
  `fv2*nut/(kappa^2*d^2)` is large and negative) -- not part of the constant
  set this tracker's architecture decision fixed, and not exercised here:
  Test 2's manufactured field was deliberately chosen to keep `S~` positive.
  **Update (Phase 3):** this was in fact hit on the very first real coupled
  `RANSFVMSolver` run, exactly as anticipated below -- the fix was added in
  Phase 3, see that phase's "What went wrong" section for the full story.
- `compute_sa_source_terms()` has no near-zero-`wall_distance` guard (`d == 0`
  would divide by zero in both the `S~` and `r` terms) -- consistent with this
  project's "no mesh-quality safeguards" stance, and not reachable by either
  test here since no cell centroid coincides with a wall.
- Only ever tested against a single wall-distance field derived from one flat
  wall face (Phase 1's own known limitation) -- Test 2's manufactured `omega`/
  `nut` fields are synthetic, but the wall-distance input they're combined
  with is not yet exercised with multiple wall patches.
- Not yet wired into any transport equation, solver, or case file --
  `RANSFVMSolver` doesn't exist yet (Phase 3's job); this phase is
  pointwise-formula infrastructure only, the same relationship Phase 0 of the
  NS tracker had to Phase 2/3's solvers.

## Phase 3 — Coupling into the viscous flux

**Added:**
- [include/RANSFVMSolver.h](../../include/RANSFVMSolver.h) /
  [src/RANSFVMSolver.cpp](../../src/RANSFVMSolver.cpp): the third solver class
  (per this tracker's architecture decision), duplicating
  `NavierStokesFVMSolver`'s inviscid flux, ghost-state, and `compute_dt()`
  structure and adding: the wall-distance field (Phase 1, precomputed once
  at construction from whichever faces' patches are `NSBoundaryType::NoSlipWall`),
  the full `nut` transport equation (advected by the coupled mean-flow
  velocity, diffused via a face-interpolated `(nu+nu_t)/sigma` and Phase 1's
  `face_gradient()`, sourced by Phase 2's `compute_sa_source_terms()` as a
  per-cell volumetric term added directly into the `nut` residual -- NOT a
  face flux), and `mu_eff`/`k_eff` substituted into the mean-flow stress
  tensor/heat conduction (`mu_eff = mu + rho*nu_t`, `k_eff = cp*(mu/Pr +
  rho*nu_t/Pr_t)`, both face-interpolated as an arithmetic mean of the two
  adjoining cells' local values, with a `NoSlipWall` boundary face's "R" side
  evaluated from `nut = 0` exactly -- so the eddy viscosity face-blends
  toward molecular-only right at a wall, not the interior cell's full value).
  `RANSBoundaryCondition` wraps `NSBoundaryCondition` unchanged, adding only
  `farfield_nut` (the one new setting the transport equation needs);
  `NoSlipWall`'s `nut = 0` and `Outflow`'s zero-order extrapolation are not
  case-configurable, per the SA model's own definition.
  `compute_dt()`'s viscous term now uses `nu + nu_t` (the larger of the two
  adjoining cells' `nu_t`, mirroring how `rho_min` already picks the more
  restrictive side for the molecular term).
- [include/SpalartAllmaras.h](../../include/SpalartAllmaras.h) /
  [src/SpalartAllmaras.cpp](../../src/SpalartAllmaras.cpp) gained the standard
  negative-`S~` robustness fix (Spalart & Allmaras 1994: `SA_CV2 = 0.7`,
  `SA_CV3 = 0.9`) -- added mid-phase after a real coupled run hit exactly the
  failure mode this fix exists for (see "What went wrong" below). This
  expands the constant set beyond what Phase 2's architecture decision
  originally fixed; re-confirmed with Mathieu before adding it, rather than
  assumed.
- `--verify-rans-stability` CLI mode (`run_verify_rans_stability()` in
  [src/main.cpp](../../src/main.cpp)): reuses `run_verify_couette()`'s exact
  1x16 sheared-mesh planar Couette setup (stationary bottom wall, top wall
  moving at `U`, zero-gradient outflow left/right), with SA genuinely turned
  on via a literature-typical freestream `nut = 3*nu`.

**What went wrong (kept in detail, matching this project's own precedent in
[docs/navier-stokes-tracker.md](../navier-stokes-tracker.md) Phase 4, for why
the finding is genuine and transferable, not throwaway debugging noise):**

The first attempt (freestream `nut = 3*nu`, matching the literature's usual
recommendation for a fully-turbulent inflow) diverged in under 120 steps --
not a slow drift, an exponential blow-up (`nut` roughly 10x-ing every 2
steps once it started). Instrumenting every step's residual and `nut`
min/max showed the divergence was NOT the `compute_dt()` risk the phase's
own plan called out in advance (an under-conservative viscous time-step
limit under a large `nu_t`) -- `nut` was still small (`O(nu)`) when the
blow-up began. Hand-evaluating the model's own formulas at the initial
condition found the actual cause: at `t=0` the flow is uniform (`omega = 0`
almost everywhere) and `chi = nut/nu = 3`, which lands `fv2` at
approximately `-1.478` -- `fv2` is genuinely negative for a wide range of
`chi` (roughly 1 to 18), a real feature of the SA model, not a bug. Near a
thin near-wall cell (`d ~ 0.03` on this mesh), the `S~ = omega +
(nut/(kappa^2*d^2))*fv2` term's near-wall singular factor amplifies that
negative `fv2` into a hugely negative `S~` (`~ -540` at the very first
step), which both flips `production`'s intended sign and, once destruction's
own `r/g/fw` chain reacts to the resulting sign/magnitude swings, produces
runaway explicit-time-stepping oscillation. This is exactly the
negative-`S~` gap Phase 2 had already flagged as a known, unimplemented
limitation ("a future case that drives `S~` negative... revisit if that's
ever actually encountered in a real coupled run") -- it was encountered on
the very first real coupled run, as anticipated.

Confirmed via targeted A/B testing (not just "add the fix and hope"): a much
smaller freestream `nut/nu` (0.1, then 1.0) that stays outside the
negative-`fv2` zone ran the full 20,000 steps with no divergence at all
(residuals settling to machine precision), isolating the cause to the
negative-`S~` mechanism specifically, not e.g. a sign error elsewhere in the
new flux assembly or an under-sized `dt`. That safe run also surfaced a
second, smaller finding: at this shear rate (`U = 0.1`, `mu = 0.02`, giving
a shear Reynolds number `~5`), `nut` monotonically decays toward zero over
the whole run regardless of the (safe) starting value -- destruction
dominates production everywhere, and the flow relaxes back to effectively
laminar. This is physically correct SA behavior at this Reynolds number
(the model has no business sustaining turbulence this weak), not a bug, but
it meant that test never actually exercised a large, sustained `nu_t` --
exactly the scenario the phase's own plan was most worried about for
`compute_dt()`. Asked Mathieu whether to (a) ship Phase 3 with the smaller,
safe freestream value and defer the negative-`S~` fix to Phase 4, or (b)
implement the standard fix now; chose (b). Adding
`SA_CV2`/`SA_CV3` and the piecewise `S~` formula (see SpalartAllmaras.h) let
the original `nut = 3*nu` setup run the full 20,000 steps with no
divergence, residuals settling to machine precision -- confirming the fix
resolves the specific mechanism found, on the specific case that triggered it.

**Test performed:** `run_verify_rans_stability()`'s final form: 1x16 sheared
mesh (`shear = 0.3`), stationary bottom wall, top wall moving at `U = 0.1`,
zero-gradient outflow left/right, `mu = 0.02`, freestream/initial
`nut = 3*nu = 0.06`, `cfl = 0.3`, 20,000 explicit steps, checking every
step for NaN/Inf and the final `rho_u`/`nut` residuals against a fixed
"settled" threshold (`1e-6`) rather than a strict step-over-step
comparison, since by the end both residuals are down at floating-point
noise level (`~1e-17`-`~1e-29`) where a strict `<=` comparison is itself
noise, not signal.

**Result:** PASS. `rho_u` residual settles to `~4.6e-17`, `nut` residual to
`~5.9e-29`, no NaN/Inf at any step. As with the safe-value diagnostic runs,
`nut` still decays toward zero over the full run (min/max both `~1e-27`-`~1e-28`
by the end) -- this specific Couette shear rate is simply too weak to
sustain turbulence under SA, confirmed physically sane rather than
suppressed by this phase's plan. The negative-`S~` mechanism is directly
exercised and survived (unlike before the fix), but a genuinely large,
sustained `nu_t` was still not reached in this test -- see "Known
limitations" below.

**Known limitations / open items:**
- **The negative-`S~` fix is not exhaustively validated.** It resolves the
  specific mechanism found here (large near-wall singular term dominating a
  negative `fv2` at `chi` in the unstable range, with `omega` small), but it
  does not make destruction's own `r/g/fw` chain well-behaved at every
  combination of `nut`/`omega`/`d` -- see SpalartAllmaras.h's own "Known
  limitation" note. Revisit if a future case (plausibly Phase 4's
  wall-resolved mesh, with much smaller `d` than this tracker's test meshes
  use) hits a different stiffness symptom.
- **No test here reaches a genuinely large, sustained `nu_t`.** This Couette
  shear rate relaxes to laminar regardless of the (safe) starting `nut`, so
  `compute_dt()`'s `nu + nu_t` viscous term -- the specific risk this
  phase's plan called out in advance -- was exercised only at small `nu_t`
  values in practice, not the "10-1000x molecular `nu`" regime the plan was
  actually worried about. That risk is deferred, not resolved; Phase 4's
  flat-plate boundary layer (which needs sustained turbulence by
  construction) is the first case likely to exercise it for real.
- `RANSBoundaryCondition::farfield_nut` is untested by this phase's
  verification (no `Farfield` boundary appears in the Couette setup, only
  `NoSlipWall`/`Outflow`) -- Phase 4's flat-plate case, which needs a real
  turbulent inflow, will be the first to exercise it.
- Not yet wired into `CaseInput`/any case-file `equation` key or `main.cpp`
  run function -- this tracker's phases never scoped that (unlike the NS
  tracker's phases, which each added case-file wiring alongside their
  solver); `RANSFVMSolver` is exercised only through
  `--verify-rans-stability` today.
- No mesh-quality safeguards, same "by design gap" stance as every prior
  phase in both trackers.

## Phase 4 (planned) — Flat-plate turbulent boundary layer validation

The standard SA verification target: a flat-plate turbulent boundary layer
compared against the log-law velocity profile,
`u+ = ln(y+)/kappa + B` (`kappa = 0.41`, `B ~ 5.0`), in the log-law region
away from the wall and away from the edge of the boundary layer.

**Known setup difficulty, flagged in advance rather than discovered
mid-debugging:**
- Needs a wall-resolved mesh (`y+ ~ 1` for the first cell off the wall) --
  a much more demanding resolution requirement near the wall than Couette
  flow's uniform mesh needed, and this project has no mesh-quality/
  stretching-ratio tooling to help build one (see MANUAL.md's "No
  mesh-quality checks").
- Needs a genuine turbulent inflow/freestream `nut` value -- SA is
  documented in the literature as sensitive to this choice; an
  inconsistent value can produce spuriously laminar or spuriously
  turbulent behavior downstream of the inlet before it "figures out" the
  right eddy viscosity.
- Needs enough streamwise extent for the boundary layer to actually
  develop before the comparison station -- an "entrance length" concern
  analogous to (but likely more demanding than) the `n_x = 1`
  translational-invariance fix the Couette validation needed for a
  different reason.

**Planned verification:** velocity profile in the log-law region matches
`u+ = ln(y+)/kappa + B` within an appropriate tolerance (to be set once a
first attempt establishes what's realistically achievable at this
project's mesh scale) -- not yet run, so no result to report.

---

## Phase 4 — Flat-plate boundary layer: attempted, target relaxed

**Added:**
- `build_flat_plate_mesh()` in [src/main.cpp](../../src/main.cpp): a flat-plate
  mesh generator (uniform streamwise cells over `[0, L]`, wall-normal cells
  geometrically stretched from the wall over `[0, H]`) -- reuses
  `build_structured_verification_mesh()` (the same generic topology builder
  every mesh in this file uses) with a `displace()` that looks up each
  reference node's precomputed stretched `y`, rather than duplicating the
  face/cell-topology code again.
- `--verify-flat-plate` CLI mode (`run_verify_flat_plate_boundary_layer()`
  in [src/main.cpp](../../src/main.cpp)).

**What went wrong (kept in detail, matching this project's own precedent in
[docs/navier-stokes-tracker.md](../navier-stokes-tracker.md) Phase 4, for why
each finding here is genuine and transferable, not throwaway debugging
noise):**

1. First attempt: `Re_L = 1e4` (chosen to keep this project's explicit,
   uniform-mesh 2D compressible time-stepping tractable -- a wall-resolved
   `y+ ~ 1` mesh, sized via a Schlichting turbulent flat-plate skin-friction
   estimate applied only to the mesh design, not the actual comparison).
   Ran stably (720 cells, ~0.83 ms/step) and converged within ~3000 of
   200,000 steps -- but `nu_t` decayed to a negligible fraction of `nu`
   (`~3e-6`), exactly like the Couette finding in Phase 3: the flow relaxes
   to an effectively laminar state rather than sustaining turbulence. The
   resulting velocity profile, compared against the turbulent log-law,
   showed a large, growing deviation (`u+` far above `ln(y+)/kappa + B`
   everywhere in the nominal log-law region, `y+` in `[30, 90]`) -- not a
   resolution problem (the exit-station profile spanned that range with 7
   cells), a genuine physical mismatch: comparing a laminar profile against
   a turbulent reference.
2. Diagnosed via the local friction Reynolds number, not by guessing: at
   `Re_L = 1e4`, `Re_tau ~ 7-8` from the simulation's own computed wall
   shear -- far below the `Re_tau ~ 150-200` minimum where a real log-law
   region can physically exist, regardless of numerics. Confirmed with
   Mathieu before spending more compute: pushed `Re_L` to `1e5` (10x, ~2x
   the mesh, ~15x the cost) as a second, real data point rather than
   theorizing further. Result: the same qualitative outcome (`nu_t/nu` still
   `~3e-6`, no positive trend at all toward sustained turbulence). Two
   consistent data points, no trend, pointed at the threshold being much
   higher than either attempt -- plausibly near real flat-plate transition
   Reynolds numbers (`~5e5-1e6`), independently estimated at roughly
   100x+ this phase's original cost (translating to hours, with no
   guarantee of success on the first try). Presented this to Mathieu with
   the concrete cost estimate; decided to relax the verification target
   rather than spend that budget chasing an uncertain outcome.
3. Relaxed target: at the (sub-transition) `Re_L` this project can run
   tractably, SA correctly predicts no sustained turbulence, so the
   mean-flow equations should produce an ordinary LAMINAR flat-plate
   boundary layer -- comparable against the classical Pohlhausen quartic
   approximation `u/U = 2*eta - 2*eta^3 + eta^4` (`eta = y/delta_99`)
   instead of the turbulent log-law originally targeted. `delta_99` is
   extracted from the simulation's own exit-station profile (not an
   independent analytic estimate), isolating the profile's self-similar
   SHAPE from any mismatch in the boundary-layer growth rate itself. Reused
   the `Re_L = 1e4` setup (already validated, cheap: ~16s for 20,000 steps).
   Result: **PASS** -- `nu_t/nu ~ 3e-6` (confirming the laminar regime), 18
   cells inside the boundary layer, `L2` error against the Pohlhausen
   profile `7.5%` (under the 15% tolerance).
4. A further finding while choosing the domain height `H`: `delta_99`
   tracks `H` itself (`H = 0.1, 0.2, 0.6` gave `delta_99/H ~ 0.92, 0.93,
   0.73` respectively, checked directly, not assumed), not a fixed value
   independent of `H` the way a genuinely isolated flat-plate boundary
   layer under an untouched freestream should behave. At this tractable Re,
   molecular viscosity is large enough that the viscous-affected region
   reaches the domain's top `Farfield` boundary well before `x = L` -- this
   setup is closer to a developing, shallow channel/duct entrance flow than
   a classical isolated flat plate with a clean freestream above it. The
   Pohlhausen comparison is still meaningful (a laminar viscous shear
   layer's self-similar shape, not specifically tied to being externally
   unbounded), but this is a real caveat on what was actually validated,
   disclosed rather than papered over. `H = 0.2` (`delta_99/H ~ 0.93`,
   `L2 = 7.5%`) was kept as the final setting since it already passes
   comfortably; making `H` large enough to avoid this entirely would need a
   proportionally larger `Re_L`/mesh, running back into finding 2's
   compute-cost wall.

**Test performed (final form):** `40x18` stretched flat-plate mesh
(`L = 1`, `H = 0.2`, first cell height `~2.6e-3` after the mesh generator's
rescale, `growth_ratio = 1.15`), `Re_L = 1e4`, freestream `nut = 3*nu`,
stationary adiabatic no-slip wall (bottom), `Farfield` inlet (left) and
undisturbed-freestream top boundary (both the same prescribed state),
zero-gradient `Outflow` (right), `20,000` steps. Exit-station (`x = L`)
velocity profile compared against the Pohlhausen quartic laminar
approximation using a self-consistently-extracted `delta_99`.

**Result:** PASS. `L2` error `7.5%`, max error `10.4%`, both under the 15%
tolerance; `nu_t/nu ~ 3e-6` confirms the comparison is legitimately against
a laminar (not turbulent) reference, consistent with what SA itself predicts
at this Reynolds number.

**Known limitations / open items:**
- **The turbulent log-law target originally set for this phase was not
  reached, and is not expected to be reachable at this project's tractable
  2D explicit-solver scale.** Two real attempts (`Re_L = 1e4, 1e5`) both
  showed no positive trend toward sustained turbulence; reaching a regime
  where SA plausibly sustains itself would need `Re_L` near real transition
  values (`~5e5-1e6`), independently estimated at 100x+ this phase's cost
  (hours, not minutes, with no guarantee of success on the first attempt).
  This is a genuine scope boundary of this codebase's explicit,
  uniform-CFL, 2D time-stepping approach, not a bug in `RANSFVMSolver` or
  the SA implementation -- see `NavierStokesFVMSolver`'s own documented 2D
  DNS scope boundary (docs/navier-stokes-tracker.md Phase 5) for a similar
  precedent of a real, disclosed architectural limit rather than a solved
  problem.
- **The flat-plate setup that DID pass is closer to a shallow channel/duct
  entrance flow than an isolated flat plate with a clean external
  freestream** (finding 4 above) -- `delta_99` tracks the domain height `H`
  rather than converging to a value independent of it. The Pohlhausen
  comparison's SHAPE match is still meaningful, but "flat-plate boundary
  layer" is a looser description of what was actually run than the phase's
  title implies.
- `RANSBoundaryCondition::farfield_nut` (added in Phase 3, untested there)
  is exercised here for the first time, but only in a regime where `nu_t`
  stays negligible -- its behavior when `nu_t` is actually significant
  remains untested by either tracker.
- No mesh-quality safeguards, same "by design gap" stance as every prior
  phase in both trackers.
- **This closes out the RANS (Spalart-Allmaras) plan as originally scoped**
  (all 4 phases attempted, each with a real verification result and a
  documented gate) -- but the plan's own stated finish line (Phase 4's
  log-law match) was not reached, by design after weighing the cost against
  the payoff, not by exhausting the available options. A future revisit
  wanting the original log-law target would need to either accept a
  multi-hour-plus run at a much higher `Re_L`, or reconsider whether this
  project's explicit 2D solver is the right tool for that specific
  validation target at all.

## Post-close addendum — `compute_dt()` viscous length-scale hardening

`RANSFVMSolver::compute_dt()`'s `length` (Phase 3's "Added" section above)
had the same `volume/face.area` deficiency on anisotropic meshes as
`NavierStokesFVMSolver::compute_dt()`, found via real Bravo integration
testing on a laminar (non-RANS) case. Since Phase 3 explicitly duplicated
`NavierStokesFVMSolver`'s `compute_dt()` structure for its own combined
laminar+turbulent viscous term, it inherited the identical exposure --
RANS's whole point is wall-resolved anisotropic meshes (e.g. Phase 4's
`build_flat_plate_mesh()` above), so this was a real, not hypothetical,
gap. Fixed identically in both solvers by replacing `length` with the new
shared `face_normal_distance()` helper
([GradientReconstruction.h](../../include/GradientReconstruction.h)). Full
root cause, fix, and verification detail (not restated here, to avoid the
two docs drifting apart) is in
[docs/navier-stokes-tracker.md](../navier-stokes-tracker.md) Phase 6. All of
this tracker's own `--verify-*` gates (`--verify-wall-distance`,
`--verify-sa-source`, `--verify-rans-stability`, `--verify-flat-plate`)
were re-run after the change and confirmed unaffected -- none of their
meshes are anisotropic enough for it to change their result.
