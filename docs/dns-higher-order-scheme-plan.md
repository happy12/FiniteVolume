# Higher-order time/space discretization for accurate 2D DNS (deferred)

Discussed 2026-07-05. Exploratory only — no implementation started or
agreed to. Phase 5's "Tier 2" of the DNS-exposure discussion in
[docs/navier-stokes-tracker.md](navier-stokes-tracker.md) — see that
document for the Tier 1 resolution diagnostic this pairs with, and for why
"DNS" here always means *fully-resolved 2D unsteady*, never 3D turbulence
(this solver is 2D-only, by deliberate, non-negotiable architectural
design — see [MANUAL.md](MANUAL.md)'s "Known limitations"; classical
Kolmogorov scaling and the energy cascade are 3D results, so nothing below
changes that ceiling, it only improves how accurately the code resolves
whatever 2D unsteady motion the flow actually has).

## Why this matters for DNS specifically (and not for the verification
work done so far)

Every existing verification test (Phases 0-4) targets a *steady* or
*quasi-1D analytic* solution — exactness at the final state is what
mattered, and the path taken to get there was irrelevant. DNS is different:
the entire *time history* of an unsteady flow's small-scale motion is the
result being computed, not just where it settles. Two properties of the
current numerics that were invisible to every steady-state test become
directly limiting there:

1. **Forward-Euler time integration is 1st-order accurate in time.**
   `EulerFVMSolver::step()`/`NavierStokesFVMSolver::step()` both advance
   `U_new = U_old + dt * RHS(U_old)` once per call. This damps and
   phase-shifts high-frequency content faster than a higher-order integrator
   would at the same `dt` — for a *steady* target this just means "converges
   a bit differently," but for genuinely unsteady turbulent-frequency content
   it means the resolved small scales decay faster than the physics says
   they should, purely from time-integration error stacked on top of the
   (separately real) spatial dissipation below.
2. **The inviscid flux is 1st-order in space everywhere** (piecewise-constant
   cell values fed directly into Rusanov/HLLC/exact — no MUSCL, no slope
   reconstruction; see [MANUAL.md](MANUAL.md)'s "Discontinuity handling").
   Rusanov in particular adds dissipation proportional to the fastest wave
   speed at *every* face, not just where a real discontinuity needs it. For
   a smooth turbulent velocity field with no shocks, this blanket dissipation
   would quietly act like an uncontrolled, un-tunable implicit LES filter —
   killing small-scale structures the mesh could otherwise resolve, in a way
   that gets *worse* as resolution increases relative to the flow's smallest
   scales (the opposite of what should happen when refining toward a DNS
   target).

Neither of these is a bug — they're exactly the right tradeoff for the
low-order, robust, shock-capturing solver this project already is. They
just aren't compatible with treating a run as a *quantitatively trustworthy*
DNS of anything beyond a smooth, low-Reynolds-number flow like the Couette
validation.

## Option A: higher-order explicit time integration (RK3/RK4)

Replace the single `U += dt * RHS(U)` forward-Euler update with a multi-stage
explicit Runge-Kutta scheme (e.g. SSP-RK3, or classical RK4), evaluating the
same `RHS()` (the existing per-face flux assembly, unchanged) 3-4 times per
step at intermediate states, then combining with the scheme's weights.

- **Scope:** touches `step()` in `EulerFVMSolver`/`NavierStokesFVMSolver` —
  restructure the current single-pass flux-assembly-plus-update into a
  reusable `compute_rhs(U) -> residual` call, then a small stage-combination
  loop around it. The flux assembly itself (inviscid + viscous, gradient
  reconstruction, boundary conditions) is untouched.
- **Cost:** 3-4x more flux evaluations per step. Some of this is offset by
  RK3/4's better stability region allowing a somewhat larger `cfl`, but
  the net cost per unit simulated time is still meaningfully higher — this
  buys time-accuracy (correct phase/amplitude of resolved unsteady content),
  not raw throughput.
- **Risk:** low. This is a mechanical, well-understood restructuring with no
  new physics and no change to spatial discretization; existing verification
  tests (all steady-state or quasi-steady) should be unaffected once the
  scheme reduces to consistent behavior at `dt -> 0`.

## Option B: higher-order, low-dissipation spatial reconstruction

Replace piecewise-constant cell values with a linear (MUSCL-style)
reconstruction of the LEFT/RIGHT face states before they enter the Riemann
solver, using the SAME per-cell gradients `GradientCalculator` already
computes (`grad_u`, `grad_v`, and for the inviscid Euler variables, gradients
of `rho`/`rho_u`/`rho_v`/`E` or their primitive equivalents) — i.e., this
reuses Phase 0's infrastructure rather than adding a new one:

```
U_face_L = U_L + grad_L . (x_face - x_L)     // extrapolate owner cell to the face
U_face_R = U_R + grad_R . (x_face - x_R)     // extrapolate neighbor cell to the face
```

then feed `U_face_L`/`U_face_R` (instead of the raw cell-centered `U_L`/`U_R`)
into the existing `rusanov_flux()`/`hllc_flux()`/`exact_riemann_flux()`
unchanged. A limiter (minmod, or an unstructured-mesh variant like
Barth-Jespersen) would still be needed near any real discontinuity to avoid
introducing oscillations — but a DNS-targeted *smooth* turbulent region
could run this reconstruction unlimited (or very lightly limited), getting
genuine 2nd-order spatial accuracy with far less blanket dissipation than
today's piecewise-constant + Rusanov combination.

- **Scope:** the larger of the two options — touches the per-face flux-assembly
  loop in both `EulerFVMSolver::step()` and `NavierStokesFVMSolver::step()`,
  adds a limiter (a genuinely new piece of numerics, not a reuse of existing
  code), and needs its own verification (e.g. re-running something like
  Phase 0's linear-field test, but checking that the *reconstructed face
  value* — not just the gradient — is exact for a linear field, plus a real
  discontinuity test confirming the limiter still captures a shock cleanly).
- **Side benefit:** this is close to what
  [docs/euler-artificial-viscosity.md](euler-artificial-viscosity.md)'s
  deferred "targeted shock sensor" discussion was reaching for, from the
  opposite direction — a limiter that only activates near real gradients,
  rather than a blanket viscosity term layered on top of Rusanov. Revisit
  that doc together with this one if either is picked up, rather than
  designing two independent dissipation-control mechanisms.
- **Risk:** medium. Reconstruction + limiting is standard, well-documented
  numerics, but is more code and more failure modes (a bad limiter can
  either oscillate near discontinuities or over-clip smooth extrema) than
  Option A.

## Recommendation (not yet agreed to; for discussion when this is picked up)

1. Don't build either speculatively — wait for a concrete case that actually
   needs it (e.g. an unsteady vortex-shedding or transitional-flow validation
   target that the current first-order-in-time, first-order-in-space scheme
   demonstrably gets wrong).
2. If/when picked up, **Option B (spatial reconstruction) first**: it
   benefits the Euler solver's existing shock-capturing use case too (not
   just DNS), directly follows up the deferred artificial-viscosity
   discussion, and reuses `GradientCalculator` infrastructure that already
   exists and is already verified (Phase 0/1). Option A (RK3/RK4) is smaller,
   lower-risk, and largely independent — could be done first instead if
   time-accuracy is the more pressing need, or in parallel since the two
   don't interact much (one changes the RHS's spatial accuracy, the other
   changes how it's integrated in time).
3. Either way, this stays scoped to "improve the accuracy of the 2D unsteady
   solution this code already computes" — it does not and cannot make this
   a 3D DNS tool. Keep that framing explicit in whatever documentation
   eventually describes this work, so a future reader doesn't infer more
   physical scope than the 2D architecture can ever support.
