# Euler solver: artificial viscosity / shock-divergence (deferred)

Discussed 2026-07-04. Exploratory only — no implementation started or agreed to.

Mathieu asked how to manage numerical dissipation in the Euler solver
([src/EulerFVMSolver.cpp](src/EulerFVMSolver.cpp)) since it can diverge on
sharp discontinuities, and whether adding artificial viscosity or diffusion
averaging (with user parameters: on/off + tunable factors) would help.

## Key findings

- The solver already uses the Rusanov flux (in `RusanovFlux.h`), which is
  itself a dissipative/artificial-viscosity-like scheme (adds a
  `Smax`-scaled diffusive term to the central flux). It's first-order,
  cell-centered, no MUSCL/reconstruction, explicit forward-Euler in time,
  CFL-limited via `compute_dt()`.
- Given that, divergence on sharp discontinuities is more likely caused by
  instability (CFL sizing, or negative pressure/density from very strong
  jumps) than by insufficient dissipation. Stacking a second, blanket
  artificial-viscosity term on top would likely just smear the whole
  solution rather than address the actual failure mode.
- Proposed direction (not yet agreed/implemented): a targeted,
  user-toggleable pressure-based shock sensor (Jameson-style) that only
  blends in extra viscosity near cells with large pressure jumps/gradients.
  This avoids dissipating well-resolved, smooth regions of the flow.
- Alternative next step offered: first reproduce/diagnose the actual
  divergence (e.g. run a strong Sod shock-tube case, check for negative
  pressure or dt blow-up) before picking a fix, rather than assuming the fix
  is "add more viscosity."

## Recommendations

1. Diagnose before fixing: reproduce the divergence on a strong Sod-tube (or
   whatever case triggers it), and check specifically whether it's a
   dt/CFL blow-up or a negative pressure/density event, rather than assuming
   lack of dissipation is the cause.
2. Don't add a second blanket artificial-viscosity term on top of Rusanov —
   Rusanov already supplies dissipation everywhere; another blanket term
   mostly just smears smooth/well-resolved regions without necessarily
   fixing the instability.
3. If a shock-capturing knob is still wanted, prefer a targeted, switched
   approach over blanket viscosity: a pressure-jump/gradient-based shock
   sensor (Jameson-style) that only activates extra dissipation near cells
   with large pressure jumps, exposed as case-file parameters
   `artificial_viscosity = on/off`, `av_kappa` (strength), and
   `av_pressure_sensor_threshold` (activation threshold) in
   [include/CaseInput.h](include/CaseInput.h), with the actual blending
   applied per-face in [src/EulerFVMSolver.cpp](src/EulerFVMSolver.cpp)'s
   flux pass.
4. Only pursue #3 after #1 confirms the failure mode is genuinely
   under-dissipated shocks (as opposed to a CFL/robustness bug), since
   fixing the wrong root cause would leave the real instability in place.
