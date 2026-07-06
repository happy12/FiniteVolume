# Navier-Stokes findings from Bravo integration testing: `compute_dt()` CFL margin + `Farfield` BC limitation

**Status: the `compute_dt()` CFL margin issue below is RESOLVED** -- see
"Resolution" under that section. The separate `Farfield` ghost-state
limitation is still open (no decision made).

## Symptom

The first real (non-synthetic-validation-mesh) Navier-Stokes run driven from
Bravo diverged almost immediately:

```
Error: solver diverged (NaN/Inf residual) at step 25
```

out of a requested 20000 steps. Case: Mach 0.98 (`ns_init = freestream 1 0.98
0 0.714285714`), `mu = 0.02`, `flux_scheme = rusanov`, `cfl = 0.5`, on a real
airfoil mesh from Bravo's `AirfoilMesherRunner` with boundary-layer cell
clustering (`blFirstLayerHeight = 0.001`, `blGrowthRatio = 1.05`,
`blDivisions = 15`, against a `tunnelHeight = 3.5`-chord farfield) —
`boundary airfoil ns_wall`, `boundary inlet ns_farfield ...`, `boundary
outlet ns_outflow`.

## Root cause (confirmed by direct bisection, not just theory)

Reran `FiniteVolume.exe` directly on the identical case file, changing only
`cfl`, 2000 steps each:

| `cfl` | Result |
|---|---|
| 0.02, 0.1, 0.25, 0.4 | Stable — residuals decrease cleanly on all 4 conserved variables |
| 0.5 | Diverges to NaN by step 25, reproduced twice |

The stability boundary is sharp and sits between 0.4 and 0.5 — not a
wrong-by-10x margin, just thin exactly where a conventionally "safe" CFL of
0.5 lands once the mesh has genuine boundary-layer-style anisotropic cells.

Suspected mechanism: `NavierStokesFVMSolver::compute_dt()` uses `length =
volume/face.area` as a single scalar cell-size proxy for both the inviscid
CFL term and the viscous `2*nu/length` term (see
[docs/navier-stokes-tracker.md](navier-stokes-tracker.md) Phase 3, which
documents this formula as validated only against a shear-skewed
(non-stretched) mesh and a small `mu=0.01` hand-built smoke mesh — never
against a mesh with real boundary-layer clustering). This proxy has no
directional awareness: on a highly anisotropic cell (very short wall-normal
dimension, much longer streamwise dimension — exactly what boundary-layer
clustering produces), `volume/face.area` likely biases toward the longer
dimension rather than the short one where the real viscous stiffness lives,
understating how restrictive the true viscous stability limit should be.

## Recommended fix

A direction-aware (or at minimum more conservative) length estimate for the
viscous stability term on anisotropic cells — e.g. the actual face-normal
cell thickness rather than `volume/face.area` — so a nominally "safe" CFL
doesn't sit this close to the real stability boundary on real (not just
synthetic validation) meshes. Worth re-running the existing
`--verify-couette`/uniform-flow checks plus a new stretched-mesh case after
any change, since Phase 4's own history (see the "what went wrong" section of
the tracker) shows this solver has already had more than one CFL-related
divergence traced to a different root cause than first suspected — don't
assume this fix is right without a targeted A/B test the way that phase
eventually did.

## Resolution

Fixed: `compute_dt()` in both `NavierStokesFVMSolver` and `RANSFVMSolver`
(which had copied the identical pattern for its own combined
laminar+turbulent viscous term) now uses a direction-aware length -- the
face-normal-projected cell-centroid separation, via a new shared
`face_normal_distance()` helper in
[GradientReconstruction.h](../include/GradientReconstruction.h)/[.cpp](../src/GradientReconstruction.cpp)
-- instead of `volume/face.area`. This is the same length scale
`face_gradient()` already uses for the viscous/diffusive flux itself, so
the stability estimate and the flux calculation are now consistent, and it
also corrects a real 2x error at every boundary face (the old formula used
a wall cell's full height instead of its centroid-to-wall distance).

Confirmed via direct A/B on a synthetic anisotropic mesh (not just "add the
fix and hope"): at a `cfl` high enough to expose the mechanism, the pre-fix
formula diverges within a handful of steps; the post-fix formula does not.
Full root cause, fix detail, and verification (including why this doesn't
contradict an earlier, unrelated `compute_dt()` finding from the Couette
validation work) is in
[docs/navier-stokes-tracker.md](navier-stokes-tracker.md) Phase 6 -- not
restated here.

Not independently re-measured against Bravo's actual airfoil mesh/case file
(not part of this repo) -- the fix was verified against a synthetic mesh
built to exercise the same anisotropic-cell mechanism, not the original
case itself. The underlying formula change is the same either way, but
Bravo's own client-side `cfl` default (lowered 0.5 -> 0.3 as a workaround
before this fix landed) does not need to be raised back based on anything
tested here.

## Related, separately-flagged limitation: `Farfield`'s ghost state has no characteristic branching

While reading `NavierStokesFVMSolver.cpp` for the above, noticed
`ghost_state()`'s `Farfield` branch returns the prescribed `farfield_state`
unconditionally, regardless of local Mach number — a fixed-state Dirichlet
condition, not a proper subsonic/supersonic characteristic-based farfield BC.

This is **not implicated** in the divergence above (lowering `cfl` alone
fixed it with Mach 0.98 and everything else unchanged), so it's not a
confirmed bug — but it's flagged separately because near-Mach-1 farfields are
exactly the regime where a fixed-state condition is a real over-specification
(one characteristic should be outgoing, not fixed), and unlike
`NSBoundaryType::Outflow`'s documented non-periodicity, this gap isn't called
out anywhere in `docs/MANUAL.md`'s "Known limitations" — it may be an
undocumented gap rather than a deliberately accepted one, worth a decision
either way.

## Not in scope for this note

Bravo has already lowered its own default `cfl` from 0.5 to 0.3 (for both
Euler and Navier-Stokes) as an immediate mitigation on its side — that's a
workaround in the client, not a fix to `compute_dt()` itself. This note is
about whether/how the solver's own stability estimate should be hardened so
future callers (Bravo or otherwise) don't have to discover this margin
empirically on their own real meshes.
