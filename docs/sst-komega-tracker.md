# k-omega SST implementation tracker

Tracks progress toward a second RANS turbulence closure, Menter's k-omega
SST, alongside the existing `RANSTurbulenceSASolver` (Spalart-Allmaras, see
[docs/archive/rans-spalart-allmaras-tracker.md](archive/rans-spalart-allmaras-tracker.md)).
Same discipline as that tracker and
[docs/navier-stokes-tracker.md](navier-stokes-tracker.md): each phase gets
its own verification gate before the next one starts, and this file records
what was added, what was actually tested, and what's still a known gap or
open question when a phase closes.

## Why SST, and why it matters for later work

k-omega SST (Menter 1994; revised Menter, Kuntz & Langtry 2003) is, along
with Spalart-Allmaras, one of the two default RANS closures for external
aerodynamics in the literature and in industrial practice (NASA's
Turbulence Modeling Resource and the AIAA CFD Drag Prediction Workshops run
both as their baseline closures). It generally out-performs SA in adverse
pressure gradient / separation-prone flow, at the cost of two transport
equations instead of one and a materially harder wall boundary condition.

**SST is also the standard base for the γ-Reθ (Langtry-Menter 2009)
transition model**, which directly targets a gap already found in this
project's own SA validation work: `airfoilNACA4412.fvmesh` was found
unsuitable for a fully-turbulent closure specifically because the flow at
its intended Reynolds number is dominated by a laminar separation bubble
(see the `rans-validation-campaign` memory entry). If this SST
implementation keeps strain-rate magnitude, vorticity magnitude, wall
distance, and the F1 blending function available as clearly-separated
quantities (not fused/optimized away), adding γ-Reθ later is two more
transport equations plus a production/destruction hook into the k-equation
-- not a rewrite. Keep that in mind when structuring Phase 1-2 below.

## Naming, following the 2026-07-06 rename precedent

Per the same rename that turned `RANSFVMSolver` into `RANSTurbulenceSASolver`
(disambiguating the closure model in the class name once a second RANS
closure was on the table): this solver is **`RANSTurbulenceSSTSolver`**,
PascalCase, matching every other solver class in this codebase
(`EulerFVMSolver`, `NavierStokesFVMSolver`, `AdvectionDiffusionFVMSolver`,
`RANSTurbulenceSASolver`). That rename landed in two passes (2026-07-06),
settling on a mixed convention -- match it exactly, don't invent a
different placement for `SST`:
- Enum tags use an underscore-separated `_SA` suffix: `EquationSet::RANS_SA`,
  `CheckpointEquation::RANS_SA`, case-file value `equation = rans_sa`. The
  `SST` equivalents follow the same underscored form:
  `EquationSet::RANS_SST` (CaseInput.h), case-file value `equation = rans_sst`,
  `CheckpointEquation::RANS_SST` (Checkpoint.h) -- per Checkpoint.h's own
  documented precedent (adding the `RANS_SA` tag never required a version
  bump, since the tag is stored as a plain numeric value and the 40-byte
  header layout doesn't change), adding `RANS_SST = 5` does **not** require
  bumping `Checkpoint::FORMAT_VERSION` either. Only bump it if this work
  also changes the fixed header layout itself.
- Case-file boundary/init keywords, the per-patch BC struct, and the
  `main.cpp` driver/writer functions instead run `SA`/`SST` directly onto
  `rans`/`run_rans`/`write_rans_fields`/`RANSBoundaryCondition` with no
  underscore in between, then an underscore before whatever comes after:
  `ransSA_wall`/`ransSA_farfield`/`ransSA_outflow`/`ransSA_init` (case-file
  keywords), `RANSBoundaryConditionSA`/`RANSBoundaryConditionSpecSA`
  (solver-side/`CaseInput`-side structs), `run_ransSA()`/`write_ransSA_fields()`
  (main.cpp driver/writer). The `SST` equivalents should mirror this exactly:
  `ransSST_wall`/`ransSST_farfield`/`ransSST_outflow`/`ransSST_init`,
  `RANSBoundaryConditionSST`/`RANSBoundaryConditionSpecSST`,
  `run_ransSST()`/`write_ransSST_fields()` -- **not** `rans_sst_wall` or
  `run_rans_sst()`, which would silently reintroduce the inconsistency this
  rename was done to remove.

Confirm this exact mapping with Mathieu before implementing Phase 3/4 (where
these identifiers actually get created) rather than re-deriving it from
scratch -- it was deliberated over several rounds this session specifically
because the "obvious" `_sst` suffix was considered and rejected in favor of
the `SST`-directly-after-`rans` form above.

**Confirmed by Mathieu (2026-07-06):** this exact mapping, including
`run_ransSST()`/`write_ransSST_fields()` as genuinely separate sibling
functions to `run_ransSA()`/`write_ransSA_fields()` (not a single generic
function branching internally on `EquationSet`) -- created when Phase 3/4
actually wire the solver in, not before.

## Architecture decision to confirm before Phase 1

Following this project's "one class per equation set" precedent (Euler vs.
Navier-Stokes vs. RANS-SA are separate classes even though each later one
reuses the previous one's inviscid flux unchanged): **default assumption is
a new `RANSTurbulenceSSTSolver` class**, duplicating
`NavierStokesFVMSolver`'s inviscid/viscous flux, ghost-state, and
`compute_dt()` structure exactly as `RANSTurbulenceSASolver` does today, not
a `turbulence_model` toggle on either existing RANS solver. Confirm this
explicitly with Mathieu before writing code -- don't assume silently, the
same way every prior solver-class decision in this codebase was confirmed
first.

**Confirmed by Mathieu (2026-07-06):** new `RANSTurbulenceSSTSolver` class, as above.

## On literature constants: verify before coding, don't hand-derive from memory

The equation **forms** below (production/destruction/cross-diffusion
structure, blending function roles, eddy-viscosity limiter) are standard
and well-established. The **numeric constants** are a different matter:
published SST implementations disagree at the 3rd-4th decimal depending on
which paper/revision/code they trace to (e.g. `gamma1` is derived as
`beta1/beta_star - sigma_omega1*kappa^2/sqrt(beta_star) ~ 0.553` in Menter's
own derivation, but several widely-used codes hardcode `5/9 ~ 0.556`
instead), and there are two genuinely different eddy-viscosity-limiter
variants in circulation (original SST uses vorticity magnitude `Omega` in
the limiter; the later "SST-V" variant uses strain-rate magnitude `S`
instead -- NASA's Turbulence Modeling Resource documents both as distinct
verification cases with different reference results). **Pull the exact
constant set and limiter variant from a single primary source before
coding** -- Menter (1994), Menter/Kuntz/Langtry (2003), or NASA's
Turbulence Modeling Resource SST page (turbmodels.larc.nasa.gov) -- and
record which one this implementation follows in this file once chosen, the
same way the SA tracker recorded its exact `Cv2`/`Cv3` provenance
(Spalart & Allmaras 1994 / Allmaras, Johnson & Spalart 2012). Do not mix
constants from different sources.

**Primary references:**
- Menter, F.R. (1994), "Two-Equation Eddy-Viscosity Turbulence Models for
  Engineering Applications," *AIAA Journal* 32(8), pp. 1598-1605.
- Menter, F.R., Kuntz, M., Langtry, R. (2003), "Ten Years of Industrial
  Experience with the SST Turbulence Model," *Turbulence, Heat and Mass
  Transfer 4*.
- Wilcox, D.C. (2006), *Turbulence Modeling for CFD*, 3rd ed. -- background
  on the underlying k-omega/k-epsilon formulations SST blends between.
- Langtry, R.B., Menter, F.R. (2009), "Correlation-Based Transition
  Modeling for Unstructured Parallelized CFD Codes," *AIAA Journal* 47(12)
  -- the future γ-Reθ extension this architecture should stay compatible
  with (see above).
- NASA Turbulence Modeling Resource (turbmodels.larc.nasa.gov) -- exact
  constant tables for both SST and SST-V, plus reference verification cases
  (flat plate, bump-in-channel) usable as an external cross-check independent
  of this project's own `--verify-sst-*` gates below.

**Confirmed by Mathieu (2026-07-06): NASA Turbulence Modeling Resource**,
now at `tmbwg.github.io/turbmodels/sst.html` -- the historical
`turbmodels.larc.nasa.gov` URL 301-redirects to a generic `nasa.gov` landing
page which itself points to that new site; recording the working URL here so
a future reader doesn't hit the same dead-end. Exact constant table pulled
from that page:

| Constant | Value |
|---|---|
| `beta_star` | 0.09 |
| `kappa` | 0.41 |
| `a1` | 0.31 |
| `sigma_k1` | 0.85 |
| `sigma_k2` | 1.0 |
| `sigma_omega1` | 0.5 |
| `sigma_omega2` | 0.856 |
| `beta1` | 0.075 |
| `beta2` | 0.0828 |

`gamma1`/`gamma2` are derived, not independent (`gamma_i = beta_i/beta_star -
sigma_omega_i*kappa^2/sqrt(beta_star)`), giving `gamma1 ~ 0.553` and
`gamma2 ~ 0.44` -- confirming this source uses the derived value, not the
`5/9 ~ 0.556` some codes hardcode instead (see the "verify before coding"
note above).

**Also confirmed by Mathieu (2026-07-06): implement BOTH eddy-viscosity
limiter variants**, not one -- `SSTLimiterVariant::Vorticity` (original
Menter 1994 SST: `nu_t = a1*k / max(a1*omega, Omega*F2)`, production clipped
at `20*beta_star*rho*omega*k`) and `SSTLimiterVariant::StrainRate` (the
Menter/Kuntz/Langtry 2003 revision: `nu_t = a1*k / max(a1*omega, S*F2)`,
production clipped at `10*beta_star*rho*omega*k`), selectable by the caller
rather than hardcoded. Each variant's production-clip coefficient is paired
with its limiter choice (not independently configurable), per NASA TMR's own
pairing of the two -- this keeps each variant internally consistent with a
single revision rather than mixing a clip coefficient from one with a
limiter from the other.

## Planned phases

| # | Item | Status |
|---|------|--------|
| 1 | k/omega transport equations (production/destruction/cross-diffusion) + F1/F2 blending, in isolation | Closed -- `--verify-sst-source` PASS |
| 2 | Eddy-viscosity closure + coupling into viscous flux and `compute_dt()` | Closed -- `--verify-sst-stability` PASS |
| 3 | Wall/farfield boundary conditions + Kato-Launder production limiter | Closed -- Farfield gate + full `CaseInput`/`main.cpp`/`Checkpoint` wiring, case-file-driven smoke-tested |
| 4 | Flat-plate validation | Closed -- `--verify-sst-flat-plate` PASS (laminar Pohlhausen target, as predicted) |

## Phase 1 (planned) -- k/omega transport equations in isolation

Implement production, destruction, and cross-diffusion source terms for the
`k` and `omega` transport equations, plus the `F1`/`F2` blending functions,
as free functions/a struct mirroring `SpalartAllmaras.h`'s
`compute_sa_source_terms()` pattern (constants collected into a struct
rather than fixed globals, so a case-file override is possible later, same
rationale as `SAModelConstants`).

`F1` blends between the inner (near-wall, k-omega) and outer (far-field,
k-epsilon-equivalent) constant sets based on wall distance, `k`, `omega`,
molecular viscosity, and the cross-diffusion term `CDkomega` -- this needs
`compute_wall_distance()` (`WallDistance.h`), reused unchanged from the SA
work, since it depends only on mesh geometry and which faces are walls, not
on which solver/model is asking.

**Planned verification:** a new `--verify-sst-source` CLI gate, mirroring
`run_verify_sa_source_terms()`: manufactured `k`/`omega` fields on
`build_sheared_verification_mesh()` (reused from the NS/SA trackers, purely
for its geometry), checking each production/destruction/cross-diffusion term
and both blending functions against their analytic values at a few
hand-picked points -- not yet run, so no result to report.

**Closed 2026-07-06, `--verify-sst-source` PASSES.** Implemented as
`SSTKOmega.h`/`.cpp` (mirroring `SpalartAllmaras.h`/`.cpp`'s free-function
shape): `SSTModelConstants`, `SSTLimiterVariant`, `sst_gamma1()`/`sst_gamma2()`,
`sst_strain_rate_magnitude()`/`sst_vorticity_magnitude()` (both exposed as
their own callable functions per this file's "Why SST matters" note),
`sst_F1()`/`sst_F2()`, `sst_eddy_viscosity()`, and `compute_sst_source_terms()`
returning an `SSTSourceTerms` struct with `k_production`/`k_destruction`/
`omega_production`/`omega_destruction`/`cross_diffusion` plus the intermediates
(`S`, `Omega`, `F1`, `F2`, `nu_t`).

**Resolved scoping note:** the k-equation's own production term (`Pk =
min(mu_t*S^2, clip_coef*beta_star*rho*omega*k)`, `mu_t = rho*nu_t`)
structurally requires the eddy viscosity even though this tracker nominally
assigns "eddy-viscosity closure" to Phase 2. `compute_sst_source_terms()`
computes `nu_t` internally regardless, since `Pk` cannot be evaluated
without it -- Phase 2's remaining scope is coupling this SAME `nu_t` formula
into the mean-flow viscous flux and `compute_dt()`, not deriving `nu_t` for
the first time. The `omega`-equation's production term, by contrast, reduces
algebraically to `gamma*rho*S^2` (the `nu_t` in Menter's `(gamma/nu_t)*P`
cancels via `mu_t = rho*nu_t`, `P = mu_t*S^2`) and never needed `nu_t` at all.

**Genuine deviation from `SASourceTerms`'s precedent, not an inconsistency:**
every SST term here is `rho`-weighted (`rho` appears explicitly), unlike
`SASourceTerms`, which omits `rho` entirely. SA's `nut` is not a mass-specific
quantity in its classical (Spalart & Allmaras 1992) form; `k`/`omega` are
conventionally transported in Menter's compressible `rho`-weighted form (the
NASA TMR equations this module follows literally include `rho` in every
term). Recorded here so a future reader doesn't "fix" this as a bug by
analogy with SA.

**Verification methodology actually used:** Test 1 mirrors SA's zero-fixed-point
check (`k=omega=0`, zero mean flow, stays exactly zero under 1000 explicit-Euler
steps). Test 2 uses a manufactured LINEAR mean-flow velocity field
(`u=a*x+b*y`, `v=c*x+d*y`) and linear-in-y `k`/`omega` fields on the 1x16
sheared mesh, both reconstructed EXACTLY by `GradientCalculator`'s
Least-Squares scheme (per `GradientReconstruction.h`'s own documented
exactness for linear fields) -- giving fully independent, hand-computable
exact `S`/`Omega` values to check against (`max` error ~1e-16, machine
precision). `F1`/`F2` themselves are checked against their known asymptotic
limits (-> 1 approaching the wall, -> 0 far from it) at the nearest/farthest
cell rather than against a hand-derived value at an arbitrary point (which
would just re-derive the same formula) -- confirmed both limits hold for both
`SSTLimiterVariant` options. Production/destruction/cross-diffusion are
checked two ways: finite-and-non-negative on every cell, and via an
independent re-aggregation formula built from the RETURNED `F1`/`F2`/`nu_t`
intermediates (max error 0, i.e. exact) -- this catches an aggregation bug
(wrong blend, missing clip, sign error) even though it doesn't independently
re-derive `F1`/`F2` themselves at that point. A separate standalone check of
`sst_eddy_viscosity()` with hand-picked `k=omega=1`, `F2=1`, `S=1.0`,
`Omega=0.1` confirms the two `SSTLimiterVariant` options give genuinely
different `nu_t` (1.0 vs. 0.31, exact) -- the manufactured field above turned
out unsuitable for this specific check, since its `omega` values are large
enough that the `a1*omega` floor dominates `max(a1*omega, limiter*F2)` for
BOTH variants everywhere in that field, which would have passed vacuously.

## Phase 2 (planned) -- Eddy viscosity closure + coupling

`nu_t = a1*k / max(a1*omega, S_or_Omega*F2)` (Bradshaw shear-stress
limiter -- confirm `S` vs `Omega` per the "verify before coding" note
above before implementing). Wire `nu_t` into the viscous flux exactly as
`RANSTurbulenceSASolver::step()` does today (`mu_eff = mu + rho*nu_t`,
`k_eff` via `prandtl_t`), and extend `compute_dt()` for the new diffusion
terms' stability limit using the existing `face_normal_distance()` helper
(the same fix already applied to both `NavierStokesFVMSolver::compute_dt()`
and `RANSTurbulenceSASolver::compute_dt()` post-close, per
[docs/ns-cfl-margin-and-farfield-bc-findings.md](ns-cfl-margin-and-farfield-bc-findings.md)).

**Planned verification:** a new `--verify-sst-stability` CLI gate, mirroring
`run_verify_rans_stability()`'s planar-Couette stress test on the same 1x16
sheared mesh, watching specifically for near-wall `omega` gradient
stiffness (a different failure mode than SA's `S~`/`fw`/`fv2` chain --
`omega`'s `1/d^2` near-wall behavior is structurally stiffer, so do not
assume SA's fix or its `compute_dt()` margin transfers unmodified) -- not
yet run, so no result to report.

**Closed 2026-07-06, `--verify-sst-stability` PASSES for both `SSTLimiterVariant`
options.** Implemented `RANSTurbulenceSSTSolver`/`RANSBoundaryConditionSST`
(`include`/`src`/`RANSTurbulenceSSTSolver.h`/`.cpp`), duplicating
`NavierStokesFVMSolver`/`RANSTurbulenceSASolver`'s inviscid/viscous
flux/ghost-state structure per the confirmed architecture decision.

**Real architecture decision made during this phase, recorded here since the
tracker's own phase split didn't anticipate it:** `k`/`omega` are transported
as the CONSERVED quantities `rho*k`/`rho*omega` (matching how `EulerState`
conserves `rho_u`/`rho_v`, not `u`/`v`, directly) -- required by Phase 1's own
decision that `SSTSourceTerms` are `rho`-weighted (Menter's compressible
form), unlike `RANSTurbulenceSASolver`'s `nut` (never `rho`-weighted, so SA
never faced this question). `RANSTurbulenceSSTSolver` does NOT store
`rho*k`/`rho*omega` as persistent members, though: `k`/`omega` (PRIMITIVE,
same external shape as SA's `nut`) are the stored/exposed fields, and
`step()` does the `rho*k`/`rho*omega` bookkeeping as local scratch values
inside its per-cell integration loop -- compute `rho*k_old = U[c].rho*k[c]`
before updating `U[c]`, add the (already `rho`-weighted) flux+source
residual, then divide by the FRESHLY-updated `U[c].rho` to recover primitive
`k`/`omega`. This keeps the external interface (accessors, `set_field()`)
identical in shape to SA's, and keeps `k`/`omega`'s own diffusive flux
differentiating PRIMITIVE `k`/`omega` (via `GradientCalculator`, exactly like
NS/SA differentiate primitive `u`/`v`/`T`, not `rho_u`/`rho_v`/`E`), while the
advective flux upwinds the CONSERVED `rho*k`/`rho*omega` directly (same
upwind structure as SA's `nut_advective`, just on the conserved quantity).

**Second real decision: the diffusive-flux model form genuinely differs from
SA's, not just in whether `rho` is present.** SA's nut diffusion divides the
WHOLE `(nu_lam+nu_t)` sum by a single constant `sigma`:
`(1/sigma)*div((nu+nut)*grad(nut))`. SST's k/omega diffusion instead scales
ONLY the turbulent part by a blended (`F1`-weighted) `sigma_k`/`sigma_omega`,
leaving molecular viscosity unscaled: `div((mu+sigma_k*mu_t)*grad(k))`. Implemented
as `mu_eff_k = mu + sigma_k_face*turb_dyn_visc_face` (dynamic form, mirroring
NS's own `mu_eff` pattern), with `sigma_k_face`/`sigma_omega_face` blended
from each face's arithmetic-mean `F1` (same face-interpolation convention as
`turb_dyn_visc_face` itself).

**`compute_dt()`'s diffusion term needed no new sigma-scaled candidate,
despite the tracker's "extend compute_dt()" framing suggesting one might be
needed:** since `sigma_k`/`sigma_omega` are both `<= 1.0` everywhere (see the
NASA TMR constant table above), `nu_lam + sigma*nu_t` never exceeds
`nu_lam + nu_t` -- the SAME (unscaled) term `RANSTurbulenceSASolver::compute_dt()`
already uses. Reusing that exact term (just substituting SST's own `nu_t`)
is therefore already a safe, conservative bound for k/omega's diffusion
stability too, not merely "close enough."

**A genuinely new term WAS needed, but for the SOURCE terms, not diffusion --
the real near-wall stiffness this phase's planned verification anticipated
turned out to live there, not in the diffusion operator.** `omega`'s wall
boundary value (`60*nu/(beta1*d1^2)`) is large on a well-resolved mesh (on
this tracker's own 1x16 test mesh, `d1 = 0.03125` gives `omega_wall ~ 16400`),
and once the near-wall cell's own `omega` rises toward that value via
diffusion, `omega`'s destruction term (`beta*omega` as a per-unit-mass rate)
becomes large enough that a plain diffusion-CFL `dt` can overshoot it in a
single explicit step. `compute_dt()` now includes a `dt_source =
cfl / (max(beta_star, beta1, beta2) * omega_max)` candidate per face (using
the largest of the three constants as a safe stand-in for the true
`F1`-blended `beta`, avoiding a fresh gradient recompute inside `compute_dt()`
purely for this) -- this is what actually kept the run stable; removing it
during development reproduced exactly the divergence this phase's planned
verification predicted.

**Also needed and not previously called out: `compute_dt()`'s existing
`nu_t_max` term (used for the mean-flow viscous CFL limit) can't be computed
exactly without fresh velocity gradients** (`nu_t` depends on `S`/`Omega`,
unlike SA's `nu_t` which only depends on `nut`/`nu`). Resolved with a cheap,
exactly-derived, always-conservative bound instead of recomputing gradients
here: since `nu_t = a1*k / max(a1*omega, limiter*F2)` and the `max()`'s
denominator is always `>= a1*omega`, `nu_t <= k/omega` unconditionally -- used
directly, no `a1` needed (it cancels).

**Wall BC methodology note:** `omega_wall`'s required per-face `d1` (distance
from the wall face to the first interior cell's centroid, NOT the general
`compute_wall_distance()` field SA uses) reuses `GradientReconstruction.h`'s
existing `face_normal_distance()` rather than adding a new function -- it
already computes exactly this quantity (verified algebraically: on this
tracker's flat, un-rotated-normal wall faces, it reduces to that cell's own
`y_centroid`, matching `run_verify_wall_distance()`'s own established result
for the same mesh). This satisfies the "needs a genuinely new per-face `d1`
value" requirement functionally (a different tool than SA's wall-distance
field) without a new geometry function.

**Observed finding, consistent with (not contradicting) this tracker's own
Phase 4 risk callout:** on this phase's mild Couette setup (`U = 0.1`), `k`
decays to numerical zero (`~1e-40`) by the end of the run for BOTH limiter
variants -- the same decay-to-laminar signature the SA tracker's Phase 4 and
the `rans-validation-campaign` memory entry already found, surfacing here
even earlier (in a basic stability smoke test, not a dedicated log-law run).
`RANSTurbulenceSASolver`'s own `nut` field decays the same way on this exact
setup (confirmed by re-running `--verify-rans-stability`, `nut` settling to
`~1e-12`) -- this is a property of the mild test flow, not specific to the
SST implementation. Both mean-flow and k/omega residuals settle cleanly
(final residuals `~1e-12` to `~1e-40`), and `omega` settles to a bounded,
non-uniform near-wall profile (max `~1339` at this `mu`/mesh, well below the
`~16400` wall BC value, consistent with a diffusion-destruction equilibrium
decaying inward from the wall) rather than diverging -- satisfying this
phase's actual bar (stable, settles; not a log-law accuracy claim).

## Phase 3 (planned) -- Boundary conditions

- Wall: `k_wall = 0` (Dirichlet); `omega_wall = 60*nu/(beta1*d1^2)` where
  `d1` is the wall-face-to-first-cell-centroid distance. This is
  analytically singular and is expected to be the single hardest piece of
  this implementation -- budget real iteration here, comparable to (per
  this project's own precedent) the Couette validation's debugging journey
  in [docs/navier-stokes-tracker.md](navier-stokes-tracker.md) Phase 4.
  Needs a genuinely new per-face `d1` value, not the general wall-distance
  field already computed for SA (which reports distance from arbitrary mesh
  points, not specifically the first off-wall cell centroid's distance).
- Farfield: `k_inf`/`omega_inf` derived from freestream turbulence intensity
  `Tu` and an eddy-viscosity ratio `nu_t/nu`, following the same
  case-file-key convention SA already established for `initial_nut`
  (new keys, naming TBD -- e.g. `sst_turbulence_intensity`,
  `sst_eddy_viscosity_ratio`).
- Kato-Launder production limiter (clips `Pk` to prevent stagnation-point
  turbulence over-production) -- SST-specific, no SA analogue.
- Wire into `CaseInput`/`main.cpp` (`equation = rans_sst`), new
  `CheckpointEquation::RANS_SST` tag, VTK output fields (`k`, `omega`,
  `nu_t`), and `docs/MANUAL.md`.

**Planned verification:** extend `--verify-sst-stability` (or a new gate) to
exercise a `Farfield` boundary, since Phase 2's Couette-only setup has no
farfield patch and would leave `k_inf`/`omega_inf` completely untested --
same gap the SA tracker's Phase 3 explicitly flagged for `farfield_nut` at
the time -- not yet run, so no result to report.

**Closed 2026-07-06.** Scope actually completed, in order:

1. **Kato-Launder limiter** -- added as a `kato_launder` trailing parameter to
   `compute_sst_source_terms()` (default `false`, every Phase 1/2 call site
   unaffected), replacing production's `S^2` with `S*Omega` in both
   `k_production` (pre-clip) and `omega_production`'s `gamma*rho*S^2`
   simplification, consistently -- both derive from the same Boussinesq `P`.
2. **Farfield verification** -- extended `--verify-sst-stability` with a third
   scenario (one `SSTLimiterVariant`, since this is a plumbing check, not a
   second full sweep): the same mesh with the left `Outflow` patch replaced by
   a `Farfield` inflow carrying prescribed `farfield_k`/`farfield_omega`.
   **PASSES** -- stays finite, settles, closing the exact gap this section
   flagged.
3. **`CaseInput`/`main.cpp`/`Checkpoint` wiring** -- `EquationSet::RANS_SST` +
   `equation = rans_sst`; `RANSBoundaryConditionSpecSST`; `ransSST_init`/
   `ransSST_wall*`/`ransSST_farfield`/`ransSST_outflow` case-file keywords
   (`ransSST_farfield` takes 6 trailing values, `rho u v p farfield_k
   farfield_omega`, explicit like SA's `farfield_nut` -- NOT derived from
   `sst_turbulence_intensity`/`sst_eddy_viscosity_ratio`, see below);
   `sst_limiter`/`sst_kato_launder`/`sst_beta_star`/etc. keys; `CheckpointEquation::RANS_SST = 5`
   with a new triple-payload `Checkpoint::write()`/`read()` overload (`U`
   immediately followed by `k`, immediately followed by `omega` -- no format-version
   bump, per Checkpoint.h's own precedent); `run_ransSST()`/`write_ransSST_fields()`
   in `main.cpp` (the exact names this file's own Naming section specified,
   genuine siblings to `run_ransSA()`/`write_ransSA_fields()`, per Mathieu's
   confirmed direction), adding `k`/`omega`/`nu_t` VTK fields and
   `residual_k`/`residual_omega` columns; wall diagnostics wired identically to
   RANS (Spalart-Allmaras) (`mu + rho*nu_t` effective viscosity, `nu_t`
   recomputed from fresh gradients since it needs `S`/`Omega`/`F2`, unlike
   SA's simpler `fv1(nut, nu)`).
4. **`initial_k`/`initial_omega` derivation** -- `sst_turbulence_intensity`
   (`Tu`)/`sst_eddy_viscosity_ratio` case-file keys derive
   `initial_k = 1.5*(Tu*|U|)^2`, `initial_omega = initial_k/(ratio*nu)` from
   the `ransSST_init` freestream velocity and `mu`/`rho`, ONLY when
   `initial_k`/`initial_omega` weren't set explicitly AND the initial
   condition is `freestream` (a `tworegion` IC has no single velocity
   magnitude to derive `Tu` from) -- a friendlier input than SA's
   directly-set-only `initial_nut`, per the tracker's original Phase 3 note.
5. **`docs/MANUAL.md`** -- new "RANS (k-omega SST) turbulence closure"
   section mirroring the SA one, plus case-file table, boundary-conditions,
   output-format, checkpointing, and file-map updates.

**Verified end-to-end via a hand-built case file** (not just the CLI verify
gates above, which construct the solver directly): a 4x4-cell square mesh
with `NoSlipWall`/`Farfield`/`Outflow` patches, `sst_limiter = strain_rate`,
`sst_kato_launder = true`, and both `initial_k`/`initial_omega` explicit AND
(in a separate run) `Tu`/ratio-derived. Confirmed: the run completes and
writes `k`/`omega`/`nu_t` VTK fields; `residual_file` reports
`residual_k`/`residual_omega` columns; `wall_forces_file`/`wall_profile_file`
produce sane per-patch/per-node rows; checkpoint-resume (200 steps, then
resumed to 400) correctly reports "Resuming from checkpoint at step 200" and
continues; the `Tu`/ratio-derived run's near-initial `k`/`omega` snapshot
matched the hand-computed expected values (`k ~ 9.2e-4` vs. `1.5*(0.05*0.5)^2
= 9.375e-4`; interior `omega ~ 4.7e-3` vs. `k/(10*nu) = 4.6875e-3`) at cells
away from the wall, with near-wall cells already visibly pulled toward the
much larger wall-BC value after just one step -- consistent with, not
contradicting, Phase 2's own finding.

**Found and fixed during this phase, worth recording:** a case-file path
starting with `/c/...` (Git Bash's MSYS-mangled form of a Windows path) fails
to open when read directly by this native (non-MSYS) executable via
`std::ifstream` -- MSYS/Git Bash only rewrites such paths when they appear as
direct command-line ARGUMENTS to a non-MSYS binary (e.g. `--validate-mesh
/c/foo`), never when they're bytes a program reads out of a file it opens
itself. Not a bug in this codebase; recorded so a future manual case-file
test doesn't lose time on the same red herring -- use native `C:/...` paths
inside case files.

## Phase 4 (planned) -- Validation

Reuse `build_flat_plate_mesh()` and the log-law target from the SA
tracker's Phase 4 as the first attempt.

**Known risk, flagged in advance rather than discovered the hard way:** the
SA flat-plate campaign's actual finding was that this project's explicit,
density-based, acoustic-CFL-limited time-stepping (no low-Mach
preconditioning or implicit stepping) could not sustain turbulence at any
`Re_L` that was tractable to run (`Re_tau ~ 7-8` achieved at both
`Re_L = 1e4` and `1e5`, versus the `Re_tau ~ 150-200` minimum where a real
log-law region can physically exist) -- `nu_t` decayed toward the laminar
limit in both attempts rather than sustaining itself (see that tracker's
Phase 4 and the `rans-validation-campaign` memory entry). `omega`'s
`1/d^2` near-wall stiffness is expected to make this **at least as hard,
plausibly harder** -- a stiffer near-wall equation does not relax the
acoustic-CFL ceiling that caused the original failure, and may force an
even smaller stable `dt` for the same wall-clock budget.

**Planned methodology:** before committing to a multi-hour run, do a short,
cheap check of whether `nu_t/nu` is growing or holding steady over the
first several thousand steps. If the same decay-to-laminar signature from
the SA campaign reappears early, stop and fall back to the same relaxed
laminar (Pohlhausen quartic profile) validation target SA landed on, rather
than re-spending a multi-hour campaign to rediscover the same ceiling. If
it does sustain turbulence, compare against the log-law
`u+ = ln(y+)/kappa + B` (`kappa = 0.41`, `B ~ 5.0`) as originally intended.

**Closed 2026-07-06, `--verify-sst-flat-plate` PASSES against the laminar
target, exactly the outcome this section anticipated.** Implemented as a
single continuous run (this project's own verify-mesh scale, 720 cells,
finished in ~20 seconds -- not the "multi-hour" case the planned methodology
was guarding against, so no early truncation was actually needed once
observed), monitoring exit-station wall-cell `nu_t/nu` every 2000 steps
throughout, mirroring `--verify-flat-plate`'s exact setup (same `L`/`H`/`n_x`/
`first_cell_height`/`growth_ratio`/`Re_L = 1e4`/`rho_inf`/`p_inf`/`u_inf`) for
a direct apples-to-apples comparison, with freestream `k`/`omega` chosen so
`nu_t/nu = 3` at the inflow (matching SA's own `initial_nut = 3*nu`
convention exactly).

**Observed `nu_t/nu` trend (one `SSTLimiterVariant`, Vorticity -- see below
for why only one):** monotonically decaying at every check, `1.62e-5` (step
2000) down to `1.40e-6` (step 20000) -- never once increasing. This selected
the laminar Pohlhausen comparison automatically (the code branches on the
observed trend, not a hardcoded assumption), which **passed**: `L2 = 0.0752`
against the quartic profile (threshold `0.15`), and the general
`compute_boundary_layer_profiles()`/Blasius-Cf/Pohlhausen-thickness-ratio
regression checks all passed too, with the identical tolerances
`--verify-flat-plate` uses.

**Cross-model consistency, not independently re-derived:** SST's converged
solution lands on essentially the SAME laminar profile SA's own
`--verify-flat-plate` converges to on the identical mesh/BCs (`L2 = 0.0752`
for SST vs. `0.0752` for SA, re-run in the same session for direct
comparison) -- exactly what should happen once both models' turbulent
viscosity has decayed to negligible: the mean-flow equations reduce to the
same laminar Navier-Stokes problem regardless of which closure decayed away,
so both should (and do) converge to the same answer. This is corroborating
evidence the decay finding is a genuine property of this project's tractable
Reynolds-number/time-stepping regime, not an SST-specific implementation bug.

**Only one `SSTLimiterVariant` was run, not both:** Phase 2's own Couette
stability check already showed `k` decaying to numerically negligible
(`~1e-40`) regardless of which limiter variant was used -- once `k`/`nu_t`
are negligible, `nu_t = a1*k/max(a1*omega, limiter*F2)`'s choice of limiter
(`S` vs. `Omega`) cannot matter to the outcome, since both terms in that
`max()` are themselves negligible. Re-running the second variant here would
not have tested anything the first didn't already establish.

**What this does and does not validate, same caveat as SA's Phase 4:** this
confirms the plumbing (BCs, source terms, `compute_dt()`'s new stability
term, wall diagnostics) holds together end-to-end on a real (if small)
external-flow-shaped mesh and produces a physically sane laminar result --
it does NOT validate SST's turbulence-sustaining behavior against any real
turbulent reference, for the identical reason SA's own Phase 4 couldn't:
this project's explicit, density-based, acoustic-CFL-limited time-stepping
cannot reach the Reynolds number a real turbulent boundary layer needs to
exist, regardless of which one-equation or two-equation closure is asked to
model it. Don't let this PASS be read as "k-omega SST validated against
turbulent flow" -- see docs/MANUAL.md's own caveat in the new SST section for
the same point stated for a reader who hasn't read this tracker.
