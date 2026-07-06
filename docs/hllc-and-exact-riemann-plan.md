# HLLC and exact Riemann (Godunov) flux schemes -- implementation plan

Discussed 2026-07-05. Planning only -- no implementation started. This
document is meant to be handed to an agent (or human) picking up the actual
implementation later; it should be self-contained enough to start from
without re-deriving the design discussion that produced it.

## Context

The Euler solver ([src/EulerFVMSolver.cpp](../src/EulerFVMSolver.cpp)) uses
a single, hardcoded numerical flux: Rusanov (local Lax-Friedrichs), via
`rusanov_flux` in [include/EulerState.h](../include/EulerState.h). Rusanov
is the most robust and cheapest option but the most dissipative -- it smears
contact discontinuities and shear layers because it applies the same amount
of dissipation (scaled to the fastest wave present) to every wave family,
including slow-moving ones that don't need it.

This plan adds two more selectable flux schemes:
- **HLLC** -- restores the contact wave that plain HLL/Rusanov discard,
  giving much sharper contact/shear resolution while remaining robust and
  positivity-preserving, at moderate extra cost/complexity.
- **Exact Riemann solver (Godunov's original method)** -- the exact solution
  to the local Riemann problem; highest possible accuracy for a first-order
  scheme, but requires a per-face Newton-Raphson iteration and is
  meaningfully more expensive (see the cost caveat below).

Goal: make the flux scheme a case-file-selectable option
(`flux_scheme = rusanov|hllc|exact`, default `rusanov` so existing case
files are unaffected), without changing anything else about the solver
(`compute_dt()`'s wave-speed estimates and `ghost_state()`'s boundary
construction are both independent of which interior flux is chosen).

## Sources

- Toro, E.F. (2009). *Riemann Solvers and Numerical Methods for Fluid
  Dynamics*, 3rd ed., Springer -- the primary reference for both schemes:
  **Ch. 4** (exact Riemann solver), **Ch. 9** (HLL family overview,
  wave-speed estimates), **Ch. 10** (HLLC).
- Toro, Spruce, Speares (1994), "Restoration of the contact surface in the
  HLL-Riemann solver," *Shock Waves* 4(1), 25-34 -- the original HLLC paper.
- Harten, Lax, van Leer (1983), "On upstream differencing and Godunov-type
  schemes for hyperbolic conservation laws," *SIAM Review* 25(1), 35-61 --
  the HLL framework HLLC extends.
- Davis, S.F. (1988), "Simplified second-order Godunov-type methods,"
  *SIAM J. Sci. Stat. Comput.* 9(3), 445-473 -- the simple `SL`/`SR`
  wave-speed estimate.
- Einfeldt, B. (1988), "On Godunov-type methods for gas dynamics,"
  *SIAM J. Numer. Anal.* 25(2), 294-318 -- refined wave-speed estimates
  (HLLE), a common robustness upgrade for HLLC.
- Sod, G.A. (1978), "A survey of several finite difference methods...,"
  *J. Comput. Phys.* 27(1), 1-31 -- the shock-tube test itself, already this
  project's existing Euler example (see `docs/MANUAL.md`'s Euler example
  case file).

## Coding structure

**Selection mechanism.** Add
`enum class NumericalFluxScheme { Rusanov, HLLC, Exact };` in
[include/EulerFVMSolver.h](../include/EulerFVMSolver.h), alongside the
existing `EulerICMode`/`EulerBoundaryType` (same category: solver
configuration, not state). `EulerFVMSolver`'s constructor gains a
`flux_scheme` parameter, stored as a member; `step()`'s per-face parallel
loop switches on it once per face -- a plain `switch`, not a
function-pointer/virtual dispatch, since the value is constant for the whole
run and this is a hot OpenMP loop (matches the existing code's avoidance of
unnecessary abstraction).

**Where the new flux functions live -- decided: option (b), one file per
scheme.** `rusanov_flux`/`hllc_flux`/`exact_riemann_flux` (moved out of
`EulerState.h`) each live in their own header --
[include/RusanovFlux.h](../include/RusanovFlux.h),
[include/HllcFlux.h](../include/HllcFlux.h),
[include/ExactRiemannFlux.h](../include/ExactRiemannFlux.h) -- leaving
`EulerState.h` focused on the state type + basic physics (`pressure`,
`sound_speed`, `flux`). All three are header-only (`inline` functions, no
`.cpp`), matching `EulerState.h`'s existing convention and preserving
inlining into `EulerFVMSolver.cpp`'s hot per-face OpenMP loop (this
project's build has no LTO, so a separate `.cpp` would prevent that).
`EulerFVMSolver.cpp` includes all three headers instead of `EulerState.h`
directly.

**Case file wiring.** New key `flux_scheme = rusanov|hllc|exact` in
[CaseInput](../include/CaseInput.h) (default `rusanov`), parsed and
validated in [CaseInput.cpp](../src/CaseInput.cpp) the same way `equation`
is (reject unrecognized values with a descriptive stderr message, `load()`
returns `false`). `main.cpp` passes `case_input.flux_scheme` straight into
the `EulerFVMSolver` constructor. No other call site changes.

## HLLC -- algorithm plan

Both new schemes need the state rotated into the face-normal/tangential
frame first (`un = u*nx + v*ny`, `ut = -u*ny + v*nx`), since the standard 1D
HLLC/exact derivations are along a single direction; the resulting
normal/tangential momentum flux is rotated back afterward. This is the
standard technique (how SU2, OpenFOAM's compressible solvers, etc.
structure it) -- Rusanov's current implementation avoids needing this only
because its dissipation term doesn't require an intermediate star-state.

1. Wave speed estimates `SL`, `SR` -- start with **Davis' simple estimate**
   (`SL = min(Vn_L - c_L, Vn_R - c_R)`, `SR = max(Vn_L + c_L, Vn_R + c_R)`),
   the simplest choice and adequate for most cases; **Einfeldt's estimate**
   is the standard robustness upgrade if Davis' proves too permissive on
   strong shocks during validation.
2. Contact speed `S*` (Toro eq. 10.37), then the star-state conserved
   variables `U*_L`/`U*_R` (Toro eq. 10.39/10.43) -- tangential momentum is
   carried through unchanged from whichever side (L/R) applies, since a
   contact wave doesn't affect tangential velocity.
3. Flux selection by which wave-speed interval contains 0 (`F_L` / `F*_L` /
   `F*_R` / `F_R`), then rotate the resulting momentum components back to
   `(x, y)`.

## Exact Riemann (Godunov) -- algorithm plan

1. Newton-Raphson solve for star-region pressure `p*` (Toro eq. 4.85/4.86,
   the pressure function and its derivative built from shock/rarefaction
   relations), seeded with a two-shock or PVRS initial guess (Toro eq. 4.47)
   for faster/more robust convergence. **Updated after initial
   implementation**: the tolerance/iteration cap (`1e-6` relative, 20
   iterations) were originally hardcoded internal constants, but are now
   user-facing case-file keys (`exact_riemann_tol`/`exact_riemann_max_iter`,
   see [docs/MANUAL.md](MANUAL.md)) at Mathieu's request, for transparency
   into the solver's behavior -- not because they're expected to need
   tuning in normal use. Same default values as before.
2. Star velocity `u*` (Toro eq. 4.9), then **sample the solution at the
   face** (always at `x/t = 0` for a stationary face, Toro Ch. 4.5) --
   determine which of the five regions (left state / left fan interior /
   star-left / star-right / right fan interior / right state) the face sits
   in, evaluate density/velocity/pressure there, and compute the physical
   flux directly from that sampled state (this *is* Godunov's flux, by
   definition).
3. **Vacuum handling**: if `2*(c_L+c_R)/(gamma-1) <= u_R - u_L`, the two
   states are separating fast enough that a vacuum forms between them and
   the standard iteration doesn't apply -- Toro Ch. 4.7 gives a closed-form
   vacuum-front solution to implement explicitly, rather than silently
   falling back to a different scheme for that face (which would make
   results a hard-to-reason-about hybrid).
4. **Cost caveat, stated plainly**: this is a Newton iteration with several
   `sqrt`/`pow` evaluations per iteration, **per face, per step**, inside
   the existing OpenMP-parallelized flux loop -- expect on the order of
   10-50x the per-face cost of Rusanov. Given the solver is first-order with
   no limiter everywhere (see `docs/MANUAL.md`'s discontinuity-handling
   section), the exact solver's accuracy benefit over HLLC is usually small
   relative to spatial truncation error at realistic mesh resolutions; its
   main value here is as a *reference*, not necessarily a routine production
   option. Implement it as a real third `flux_scheme` choice as asked, but
   the validation phase should measure and report the actual overhead on
   this codebase so that's an informed choice, not a guess.

## Validation plan

1. **Consistency check**: `U_L == U_R` (uniform flow) must reproduce the
   exact physical flux for all three schemes (`F*(U,U) = F(U)`) -- a basic
   unit-style check any correct numerical flux must pass.
2. **Reflection symmetry**: `flux(U_L, U_R, n) == -flux(U_R, U_L, -n)` for
   HLLC and Exact -- catches sign/rotation bugs.
3. **Sod shock tube against the exact solution** -- since the exact Riemann
   solver *is* the reference-solution generator (same math, used once as an
   initial-condition sample instead of per-step), this is nearly free once
   the exact solver exists: run Rusanov/HLLC/Exact on the standard Sod case
   (`docs/MANUAL.md`'s Euler example), compare each against the analytic
   profile at the same physical time.
4. **Contact-sharpness comparison**: measure how many cells span the
   density jump across the contact discontinuity for each scheme -- expect
   Exact ~= HLLC (both resolve it essentially exactly at first order) and
   both visibly sharper than Rusanov.
5. **Robustness stress tests**: a strong-shock two-region case (large
   pressure ratio) confirming no negative density/pressure for Rusanov/HLLC
   and that Exact's Newton iteration converges (or hits the documented
   vacuum path) rather than diverging or hanging; a strong double-rarefaction
   case specifically to exercise the vacuum-formation path.
6. **Regression smoke test**: rerun the existing quad and 20-gon test meshes
   (see the round-trip/validation work already done for the `.fvmesh`
   format) under all three `flux_scheme` values, confirming no crash/NaN --
   these have no discontinuity, so this mostly checks the plumbing
   (mesh/solver wiring), not the physics.
7. **Performance measurement**: time a representative Sod-tube run under
   each scheme, reporting the actual per-step overhead of HLLC and Exact
   relative to Rusanov on this codebase, to make the cost trade-off from the
   algorithm section concrete rather than a textbook estimate.

## Status

**HLLC and exact Riemann both implemented** (2026-07-05): `rusanov_flux`,
`hllc_flux`, `exact_riemann_flux` each live in their own header (see file
layout above); `NumericalFluxScheme` (`Rusanov`/`HLLC`/`Exact`) lives in
[EulerFVMSolver.h](../include/EulerFVMSolver.h); case-file key
`flux_scheme = rusanov|hllc|exact` wired through `CaseInput`/`main.cpp`.
**Full validation plan run (2026-07-05)**, via a standalone header-only
harness (compiled directly with g++, not part of the CMake target) plus
real `FiniteVolume.exe` runs for the mesh/plumbing and performance items:

1. **Consistency** (`F*(U,U)=F(U)`): PASS for all three schemes, all test
   states/normals, machine-precision errors (<=1e-16).
2. **Reflection symmetry** (`flux(L,R,n)=-flux(R,L,-n)`): PASS for all
   three, machine-precision errors.
3. **Sod tube vs. analytic** (400 cells, t~0.059): L1 density error vs. the
   exact solution -- Rusanov 6.59e-3 > HLLC 4.66e-3 > Exact 4.47e-3, the
   expected ordering. The analytic reference (generated by reusing
   `ExactRiemannDetail`'s own sampling code at general x/t, not just the
   x/t=0 used for the flux) reproduced the textbook Sod values (p*=0.30313,
   u*=0.92745, rho*_L=0.42632, rho*_R=0.26557) to 5 decimal places.
4. **Contact sharpness**: Rusanov smears the contact over 16 cells vs. 15
   for both HLLC and Exact, within a window bounded by the rarefaction
   fan's tail and the shock (to isolate the contact specifically). Real but
   modest, smaller than textbook "near-exact" resolution, because this
   solver is first-order with no limiter, so truncation error accumulates
   as the contact is advected over many steps regardless of scheme.
5. **Robustness**: strong shock (p_l/p_r = 100000) -- no negative
   density/pressure, any scheme. Strong double rarefaction (u_l=-8,
   u_r=+8, vacuum-forming) -- all three stayed positive; Exact's minimum
   density (5.3e-5) tracked the true near-zero vacuum far more closely
   than Rusanov (0.021) or HLLC (0.030), consistent with only Exact having
   genuine vacuum handling.
6. **Regression smoke test**: the existing quad Sod mesh, plus a new
   arbitrary-polygon mesh (a 20-gon and a 10-gon sharing a 7-segment
   interface, built to stress "no fixed face count") both ran clean under
   all three `flux_scheme` values -- no crash, no NaN/Inf.
7. **Performance** (Release build, 5000-cell 1D-strip mesh, 5000 steps, 2
   reps each): Rusanov ~2.37s avg, HLLC ~2.49s avg (+9%), Exact ~2.93s avg
   (+23%) -- far below the "10-50x" textbook estimate above, likely because
   the two-shock/PVRS seed converges in very few Newton iterations on this
   mild case, and per-step fixed costs (`compute_dt()`, residual reduction)
   dilute the flux-only difference. Worth re-measuring on a harder case if
   the overhead becomes a real concern later.

All validation-plan items are now complete.
