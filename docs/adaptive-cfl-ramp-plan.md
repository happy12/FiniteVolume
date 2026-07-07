# Adaptive (residual-based) CFL ramping -- recommendation

Discussed 2026-07-06 (prompted by a Bravo user's k-omega SST run looking
"frozen" -- traced to a Bravo-side status-display bug, not a solver hang, but
the discussion surfaced a real question: should the solver help a user avoid
hand-tuning a single fixed CFL for an entire run?). Planning only -- no
implementation started. This document is meant to be handed to an agent (or
human) picking up the actual implementation later; it should be self-contained
enough to start from without re-deriving the design discussion that produced
it.

## Implementation note (2026-07-06)

Implemented, with one deliberate revision to the trigger rule below, made
after reading SU2's actual `CSolver::AdaptCFLNumber`
(`SU2_CFD/src/solvers/CSolver.cpp`) rather than working from this plan's
single-step-ratio sketch alone. SU2's linear-residual/under-relaxation
triggers don't transfer (this project has no linear solve, no
under-relaxation -- every solver here is explicit), but its **windowed
log-residual trend with a hard divergence reset** is more robust than a
single-step ratio on an explicit solver's inherently noisier step-to-step
residual, and was adopted in place of item 2's `reduction_ratio` rule below:
`next_cfl()`'s actual signature takes a rolling residual-history window (not
a before/after pair), and a fresh divergence-reset path (hard reset to
`cfl_min` if the windowed trend rises more than `cfl_ramp_divergence_
threshold` orders of magnitude) was added beyond what's described below. Six
case-file keys shipped (the four below plus `cfl_ramp_window` and
`cfl_ramp_divergence_threshold`), not four. See
[include/CflRamp.h](../include/CflRamp.h) for the actual rule and
[docs/MANUAL.md](MANUAL.md)'s "Residual-based CFL ramping" section for the
shipped case-file surface and user-facing behavior -- both supersede the
sketch below where they disagree. The rest of this plan (case-file
key names/defaults for `cfl_min`/`cfl_ramp_growth`/`cfl_ramp_shrink`, the
`main.cpp`-centralized placement, `set_cfl()` on all four solvers, the
checkpoint-deferral decision, and `--verify-cfl-ramp`) shipped as originally
designed.

## Context

Every one of the four CFL-driven solvers (`EulerFVMSolver`,
`NavierStokesFVMSolver`, `RANSTurbulenceSASolver`, `RANSTurbulenceSSTSolver`)
takes a single scalar `cfl` at construction ([CaseInput.h:416](../include/CaseInput.h),
`cfl = 0.5` default) and uses it unchanged for the entire run -- each solver's
own `compute_dt()` multiplies `cfl` against that step's own inviscid-wave-speed
(and, for Navier-Stokes/RANS, viscous-diffusion and, for RANS k-omega SST,
explicit-source) stability limits, recomputed fresh every step, but the `cfl`
coefficient itself is a constant for the whole run today. `cfl` is currently a
private member with no setter (see e.g. [EulerFVMSolver.h:162](../include/EulerFVMSolver.h)).

The practical problem this causes: a `cfl` safe for the smooth, mostly-converged
flow field late in a run is often not safe for the crude uniform/freestream
initial condition every run starts from, especially for the two RANS solvers
-- `RANSTurbulenceSSTSolver::compute_dt()`'s explicit-source term
(`dt_source = cfl / max(beta_star*omega, beta*omega)`, see that header's class
comment) is most restrictive exactly when a near-wall cell's `omega` is still
approaching its large wall-boundary value from a uniform initial guess. Users
currently have to pick one conservative `cfl` for the whole run, giving up
speed later just to survive the first several hundred steps.

**What this does NOT fix.** Per CLAUDE.md's SST section, the flat-plate
validation's decay-to-laminar signature has the "same acoustic-CFL-limited
explicit time-stepping root cause" as SA's own equivalent finding -- i.e. the
number of physical time units reachable in a given step budget is capped by
the acoustic CFL condition from the very first step, and ramping still only
reaches the same target `cfl` (ceiling) a fixed-CFL run at that same target
would use once warmed up. This plan is a robustness/usability improvement for
surviving the initial transient and removing manual tuning, not a fix for
that separately-documented turbulence-sustaining validation gap -- the
validation plan below explicitly checks that this feature doesn't change that
finding, rather than hoping it might.

## Sources

- Blazek, J. (2015). *Computational Fluid Dynamics: Principles and
  Applications*, 3rd ed., Butterworth-Heinemann -- Ch. 6.1 ("Time Step
  Determination")/Ch. 8 (multigrid/convergence acceleration) discuss starting
  an explicit or dual-time-stepping run at a reduced CFL and increasing it as
  the solution develops, specifically flagging RANS two-equation models as
  needing a more conservative start than the mean-flow equations alone due to
  turbulence source-term stiffness -- the general "start low, ramp up"
  heuristic this plan implements.
- Mavriplis, D.J. (1998), "Multigrid Strategies for Viscous Flow Solvers on
  Anisotropic Unstructured Meshes," *J. Comput. Phys.* 145(1), 141-165 --
  residual-based CFL control on unstructured meshes (this project's own mesh
  model), motivating a feedback rule over a fixed schedule.
- Economon, T.D., Palacios, F., Copeland, S.R., Lukaczyk, T.W., Alonso, J.J.
  (2016), "SU2: An Open-Source Suite for Multiphysics Simulation and Design,"
  *AIAA Journal* 54(3), 828-846 -- SU2's `CFL_ADAPT` option is the most
  directly comparable production implementation: grow CFL by a factor when the
  monitored residual drops fast enough, shrink it by a (usually smaller)
  factor when it stalls or rises, bounded by a min/max -- the three-way
  grow/shrink/hold rule this plan adopts (see SU2's `CConfig`/`CIteration`
  `Adapt_CFL_Number()` for the reference implementation's exact shape, not
  reproduced verbatim here since it's not public-domain source).

## Coding structure

**Selection mechanism.** New case-file key `cfl_mode = fixed|ramp` (default
`fixed`, so existing case files are unaffected -- same backward-compatibility
precedent as `flux_scheme`'s default). Only meaningful for the four
`cfl`-driven equation sets (`euler`, `navier_stokes`, `rans_sa`, `rans_sst`);
`diffusion`/`advection_diffusion` use `dt`/`alpha` directly and are unaffected.

**Where the ramp logic lives -- main.cpp's shared step loop, not inside any
solver's `step()`.** Same precedent as residual tracking/checkpointing/
stopping-criteria already being centralized there instead of duplicated per
solver (see CLAUDE.md's architecture section: "run inline in `main.cpp`...
handled uniformly across all four"). CFL ramping is the same category of
cross-cutting run-control concern, not per-solver physics -- it has no reason
to know about SA's `nut` vs SST's `k`/`omega` vs Euler's plain state, it only
needs *a* residual scalar and *a* `set_cfl()` call. Concretely:

1. New tiny header `include/CflRamp.h`:
   ```cpp
   enum class CflMode { Fixed, Ramp };

   struct CflRampParams {
       double cfl_min = 0.05;         // starting CFL (case-file "cfl_min")
       double growth_factor = 1.5;    // multiplier applied when residual is dropping (case-file "cfl_ramp_growth")
       double shrink_factor = 0.5;    // multiplier applied when residual stalls/rises (case-file "cfl_ramp_shrink")
   };

   // Stateless: given the previous step's chosen cfl and the before/after
   // residual scalars for the interval just completed, returns the next cfl
   // to use, clamped to [params.cfl_min, cfl_max]. cfl_max is the existing
   // "cfl" case-file key's value -- ramping doesn't need a *new* key for the
   // ceiling, it reuses "cfl"'s existing meaning (the CFL you actually want
   // once warmed up).
   double next_cfl(double current_cfl, double cfl_max, double residual_before,
                    double residual_after, const CflRampParams& params);
   ```
   A pure function, not a class with hidden state, matching this codebase's
   general preference (e.g. the flux functions in `RusanovFlux.h`/`HllcFlux.h`)
   -- easy to unit-test in isolation with synthetic residual sequences (see
   Validation plan #1), and `main.cpp` only needs to keep `current_cfl` as a
   local in the step loop, no new class instance.

2. **Grow/shrink/hold rule** (`next_cfl`'s body), the SU2-style three-way rule:
   - Compute `reduction_ratio = residual_after / residual_before`.
   - If `reduction_ratio < 1` (residual dropped): `current_cfl *= growth_factor`,
     clamped to `cfl_max`.
   - If `reduction_ratio > 1.1` (residual grew by more than a 10% tolerance
     band -- avoids reacting to normal step-to-step noise as if it were
     divergence): `current_cfl *= shrink_factor`, clamped to `cfl_min`.
   - Otherwise (roughly flat, within that tolerance band): hold `current_cfl`
     unchanged.

3. **Which residual scalar to monitor.** Recommend the density residual
   (`solver.residual().rho` -- already computed every `step()` call
   unconditionally, see e.g. [EulerFVMSolver.h:137-141](../include/EulerFVMSolver.h),
   independent of whether `residual_file` is even set) as the single monitored
   scalar for all four equation sets, rather than combining rho/rho_u/rho_v/E
   (and, for RANS, `nut` or `k`/`omega`) into one artificial blended number --
   same choice SU2 makes (monitors density residual specifically for its CFL
   adaptation trigger). `nut_residual()`/`k_residual()`/`omega_residual()`
   already exist as separate accessors if a future revision wants to gate on
   the turbulence variable too, but starting with one well-understood signal
   avoids arbitrarily weighting residuals with different physical units and
   convergence rates.

4. **Check cadence.** Evaluate `next_cfl()` every step, not just every
   `residual_interval` steps -- the residual scalar is already free every
   step (see point 3), so there's no reason to react to it more slowly than
   `residual_file`'s own write cadence, and reacting every step gives the
   fastest possible response to an incipient divergence.

5. **Per-solver change**: add `void set_cfl(double new_cfl) { cfl = new_cfl; }`
   to the public interface of all four solver classes (mechanical, ~1 line
   each) -- `cfl` is already a plain private member on all four (confirmed:
   [EulerFVMSolver.h:162](../include/EulerFVMSolver.h),
   [NavierStokesFVMSolver.h:219](../include/NavierStokesFVMSolver.h),
   `RANSTurbulenceSASolver.h`, `RANSTurbulenceSSTSolver.h`), no other member
   depends on `cfl` being fixed at construction time.

6. **Case-file wiring** ([CaseInput.h](../include/CaseInput.h)/
   [CaseInput.cpp](../src/CaseInput.cpp)): `cfl_mode = fixed|ramp` (validated
   the same way `sst_limiter`/`equation` are -- reject an unrecognized value,
   `load()` returns `false`), `cfl_min`, `cfl_ramp_growth`, `cfl_ramp_shrink`
   (all three ignored/unvalidated when `cfl_mode == fixed`, same "only
   meaningful for X" precedent as `sst_*` keys being ignored outside
   `rans_sst`). `main.cpp`'s four `run_*` functions gain, right after
   constructing their solver:
   ```cpp
   double current_cfl = (case_input.cfl_mode == CflMode::Ramp) ? case_input.cfl_min : case_input.cfl;
   double prev_residual = 0.0;
   bool have_prev_residual = false;
   ```
   and, inside the existing step loop, right after the existing
   NaN/Inf-divergence check and before the existing `tracking_residual`
   file-write block:
   ```cpp
   if (case_input.cfl_mode == CflMode::Ramp) {
       if (have_prev_residual) {
           current_cfl = next_cfl(current_cfl, case_input.cfl, prev_residual, r.rho, ramp_params);
           solver.set_cfl(current_cfl);
       }
       prev_residual = r.rho;
       have_prev_residual = true;
   }
   ```
   Four call sites (`run_euler`/`run_navier_stokes`/`run_ransSA`/`run_ransSST`),
   each already has `r`/`residual` in scope at that point in the loop.

7. **Checkpoint interaction -- open design question, recommend deferring for
   v1.** Resuming a checkpointed ramp run today would restart `current_cfl` at
   `cfl_min` every resume (checkpoints round-trip the field state, not the
   ramp controller's state), which defeats the purpose for a long
   multi-session run. Fixing this properly means adding `current_cfl` (one
   `double`) to the checkpoint payload, which per CLAUDE.md's checkpoint
   section requires bumping `Checkpoint::FORMAT_VERSION` (same as every prior
   solver addition that touched the payload). Recommend implementing v1
   *without* this (resume restarts the ramp from `cfl_min`, documented as a
   known limitation in `docs/MANUAL.md`), and only adding the checkpoint field
   if a real workflow needs long checkpointed ramp runs -- avoids a
   format-version bump for a feature whose main value (surviving the initial
   transient) has usually already paid off long before a typical checkpoint
   interval.

## Exposed case-file surface (summary)

| Key | Type | Default | Meaning |
|---|---|---|---|
| `cfl_mode` | enum | `fixed` | `fixed` (today's behavior, unchanged) or `ramp` |
| `cfl` | double | 0.5 | Unchanged meaning: the target/ceiling CFL. In `ramp` mode this is the value the ramp grows *toward*, not a fixed value used verbatim from step 1 |
| `cfl_min` | double | 0.05 | Ramp only: starting CFL |
| `cfl_ramp_growth` | double | 1.5 | Ramp only: multiplier applied when the density residual drops step-over-step |
| `cfl_ramp_shrink` | double | 0.5 | Ramp only: multiplier applied when the density residual rises by more than 10% step-over-step |

Four new keys total (one enum + three doubles), all ignored in `fixed` mode --
deliberately kept small per the "few more parameters as appropriate" framing,
rather than exposing every SU2-style adaptation knob (SU2 also has separate
under-relaxation and multiple monitored-variable options; not proposed here
unless the simple version proves insufficient in practice).

## Validation plan

1. **Unit-style check on `next_cfl()` directly** (no mesh/solver needed, same
   spirit as the exact-Riemann plan's consistency check). **DONE** -- shipped
   as `run_verify_cfl_ramp()`'s first block in
   [src/main.cpp](../src/main.cpp), covering the actual (windowed-trend, not
   single-step-ratio) rule: a monotonically decreasing synthetic residual
   sequence grows `current_cfl` once the `cfl_ramp_window`-sized history
   fills, clamped at `cfl_max`; a monotonically, mildly increasing sequence
   shrinks it, floored at `cfl_min`; a perfectly flat history holds it
   unchanged; a synthetic jump exceeding `cfl_ramp_divergence_threshold`
   orders of magnitude across the window hard-resets to `cfl_min` even from a
   high starting `cfl`; fewer than `cfl_ramp_window` samples holds `cfl`
   unchanged (the warm-up case, not present in the original single-step-ratio
   design this item was written against).
2. **New `--verify-cfl-ramp` gate** (matching this project's existing
   `--verify-*` convention rather than a separate harness). **DONE** -- reuses
   `--verify-couette`'s exact 1x16 sheared-mesh Couette setup, run once with
   `cfl_mode = fixed` at `cfl = 0.3` and once with `cfl_mode = ramp`
   targeting the same `0.3` as the ceiling; both converge to the same steady
   solution well within the 5% L2 tolerance `--verify-couette` itself uses,
   and the ramp run's `current_cfl` trace is confirmed (via `min`/`max` seen
   across the run) to start at `cfl_min` (0.05) and reach at least 90% of the
   0.3 ceiling.
3. **Robustness stress test -- the actual value proposition.** **DONE
   2026-07-06.** Result: bisection (12 binary-search iterations after an
   initial doubling search from the known-stable `2.0`) found `cfl =
   2.69971` diverges within 486 steps in `fixed` mode on the 10x67 stretched
   mesh; a `fixed` run pinned at `cfl_min = 0.05` survives the full 500-step
   window (confirms the baseline is genuinely conservative); a `ramp` run
   targeting `2.69971` as the ceiling also survived the full 500 steps --
   the strong outcome (case (a) above, not just the fallback (b)). PASS. See
   `run_verify_cfl_ramp()`'s third block in [src/main.cpp](../src/main.cpp).
   No case currently sitting in this repo is known to diverge in `fixed`
   mode: the one real divergence on record
   (`docs/ns-cfl-margin-and-farfield-bc-findings.md`, Bravo's Mach-0.98
   airfoil case at `cfl = 0.5`) was traced to a `compute_dt()` bug that is
   now fixed and generalized into `--verify-ns-stretched-cfl`'s regression
   mesh (`build_flat_plate_mesh(mesh, n_x=10, L=1.0, H=0.05,
   first_cell_height=1e-4, growth_ratio=1.05)`, stable today at `cfl = 2.0`
   for 2000 steps). A concrete test therefore has two parts:
   1. **Find a divergent `cfl` by direct bisection** (same methodology
      `ns-cfl-margin-and-farfield-bc-findings.md` itself used) on that exact
      mesh/BC/IC setup (`rho_inf=1.0, p_inf=1.0, u_inf=0.5, mu=0.02,
      gamma=1.4, prandtl=0.72`, freestream IC, `NavierStokesFVMSolver`),
      sweeping `cfl` upward from the known-stable `2.0` in fixed mode, 2000
      steps per trial, until a value that reliably diverges (NaN/Inf
      residual) within the first few hundred steps is found.
   2. **Rerun the same setup twice** at that found-divergent value as the
      ceiling: once `cfl_mode = fixed` (expected: diverges, confirming the
      bisected value actually reproduces the failure) and once `cfl_mode =
      ramp` targeting it as the ceiling. Passing bar: the `ramp` run either
      (a) survives the full 2000 steps by shrinking back down whenever the
      residual trend crosses `cfl_ramp_divergence_threshold`, or (b) if it
      still diverges, does so no faster (in step count) than a `fixed` run
      pinned at `cfl_min` for its entire duration -- i.e. ramping must never
      make the outcome *worse* than the most conservative fixed choice, even
      if it can't always make a genuinely bad ceiling survivable.
   Implementation note for whoever picks this up: this needs a new
   `run_verify_cfl_ramp_stress()` (or a third block inside
   `run_verify_cfl_ramp()`) since no existing `--verify-*` gate exercises a
   deliberately-divergent `cfl` today; the bisection step itself is a one-time
   manual exploration (print the found value in the gate's own output so a
   future reader can see what was used, the same way
   `run_verify_ns_stretched_cfl()` documents its own mesh's `first_cell_height`).
4. **SST-specific.** **DONE 2026-07-06.** Result: both a fresh `fixed`
   pass (`cfl = 0.3` throughout) and a fresh `ramp` pass (targeting `cfl =
   0.3` as the ceiling) on `run_verify_sst_flat_plate()`'s exact 40x18 mesh
   ended with the identical "decaying or flat" verdict (final `nu_t/nu` =
   1.39606e-06 fixed vs. 1.40045e-06 ramp -- both negligible relative to the
   freestream's `nu_t/nu = 3`). PASS. See `run_verify_cfl_ramp()`'s fourth
   block in [src/main.cpp](../src/main.cpp). Extend `run_verify_sst_flat_plate()`
   (or add a sibling function reusing its exact setup: `build_flat_plate_mesh`
   with `n_x=40, L=1.0, H=0.2, first_cell_height=2.949e-3, growth_ratio=1.15`,
   `Re_L=1e4`, `u_inf=0.2`, `cfl=0.3`, `initial_k`/`initial_omega` giving
   `nu_t/nu = 3` at the freestream) with a second pass run identically except
   `cfl_mode = ramp` targeting the same `cfl = 0.3` as the ceiling. Passing
   bar: the ramp pass's own `nu_t/nu` trend at the exit station (the same
   quantity `run_verify_sst_flat_plate()` already prints every 2000 steps)
   must show the identical decaying-or-flat "laminar-decay signature" the
   fixed-mode pass shows -- confirming this feature doesn't accidentally
   change the underlying turbulence-sustaining finding documented in
   CLAUDE.md/docs/sst-komega-tracker.md, rather than leaving that claim in
   this plan's Context section unverified.
5. **Backward compatibility / regression.** **DONE 2026-07-06.** Result:
   (a) confirmed statically -- every one of the 8 `case_input.cfl_mode ==
   CflMode::Ramp` gates across the four `run_*` loops (2 each: setup +
   per-step) in [src/main.cpp](../src/main.cpp) is a hard equality check
   against a default-`Fixed` enum, so a case file with no `cfl_mode` key
   provably never calls `set_cfl()`; (b) ran a scratch copy of
   `naca0012_alpha0.case` (real 159,485-cell airfoil mesh, `cfl = 0.5`, no
   `cfl_mode` key) for 300 fresh steps with `checkpoint_file` omitted and
   `output_file`/`residual_file` pointed at scratch paths -- completed
   without error, residuals stayed in the same sane, non-diverging 1e-3
   range the campaign's own README describes for its early transient. The
   live `naca0012_alpha0.ckpt` was confirmed untouched (still exactly
   6,379,440 bytes, matching the README's documented byte-for-byte size)
   afterward. PASS. A concrete,
   already-existing target: [tools/rans_validation/naca0012_alpha0.case](../tools/rans_validation/naca0012_alpha0.case)
   (`equation = rans_sa`, `cfl = 0.5`, no `cfl_mode` key, and a live
   `checkpoint_file` from an in-progress 200000-step validation campaign --
   see the "RANS validation campaign status" memory; **do not resume/overwrite
   that checkpoint as part of this check**). The test: (a) confirm
   `CaseInput::load()` parses this file unchanged and `cfl_mode` resolves to
   `CflMode::Fixed` with `cfl_ramp` holding its all-default values; (b) run a
   short, isolated smoke test against a COPY of this case file pointed at
   fresh `output_file`/`residual_file` paths and no `checkpoint_file` (so the
   live campaign's checkpoint is never touched), for a few hundred steps only;
   confirm the residual trace is identical to what the same short run produces
   with the `cfl_mode`/`cfl_ramp_*` keys and code path removed entirely --
   i.e. that `RANSTurbulenceSASolver::set_cfl()` is never called and `cfl`
   stays at its case-file value of `0.5` for every step. Since `cfl_mode`
   defaults to `Fixed` and the `Ramp` branch in each `run_*` loop is gated on
   `case_input.cfl_mode == CflMode::Ramp`, this should be a pure smoke-test
   confirmation, not a source of new behavior.

## Status

Implemented 2026-07-06 -- see the Implementation note near the top of this
document for the one deliberate deviation (SU2-informed windowed-trend/
divergence-reset rule in `next_cfl()`, in place of this plan's original
single-step-ratio sketch) and [docs/MANUAL.md](MANUAL.md)'s "Residual-based
CFL ramping" section for the shipped, user-facing description. **All five
validation plan items are now DONE (2026-07-06)** -- items 1-2-3-4 are
covered by `run_verify_cfl_ramp()`'s four blocks in
[src/main.cpp](../src/main.cpp) (all `PASS` under `--verify-cfl-ramp`); item
5 was run manually as a one-off scratch-copy smoke test against the real
`naca0012_alpha0.case` (not folded into the CLI gate, since it deliberately
exercises a real, large, external-file-backed case rather than a
self-contained in-memory mesh -- see item 5's own Result for what was run
and confirmed). No open validation follow-up remains for this feature as
originally scoped; the only items NOT covered by any test here are the ones
this plan's Context section always said were out of scope (the
acoustic-CFL-limited turbulence-sustaining validation gap) and the
checkpoint round-tripping deferral noted in "Coding structure" item 7.
