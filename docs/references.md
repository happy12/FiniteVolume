# Literature references and validation provenance

Every non-trivial algorithm in this codebase that comes from published
literature, organized by the module implementing it: the full citation, what
was implemented from it, and — per this document's whole purpose — exactly
how each was validated and what it was compared against (cross-referencing
the `--verify-*` gates in [MANUAL.md](MANUAL.md)'s "Running" section, since
this project has no separate unit-test framework). Sections with no cited
source are labeled as such explicitly, rather than omitted, so this document
is a complete map of "where did this come from," not just a bibliography of
the parts that happen to cite something.

This document was compiled by searching every existing code comment and doc
for citation-style text, then independently verifying each citation's exact
bibliographic details and, in one case, its accuracy (see
[Spalart-Allmaras turbulence model](#spalart-allmaras-turbulence-model)
below) — it isn't a from-scratch literature survey, it's a consolidation and
fact-check of what this project already drew on.

## Spalart-Allmaras turbulence model

Implemented in [include/SpalartAllmaras.h](../include/SpalartAllmaras.h) /
[src/SpalartAllmaras.cpp](../src/SpalartAllmaras.cpp), coupled into a real
solver in [RANSTurbulenceSASolver](../include/RANSTurbulenceSASolver.h). Two distinct sources,
previously conflated in this codebase's own comments (corrected as part of
compiling this document — see the note below):

- **Base model** (`sa_fv1`/`sa_eddy_viscosity`, the `production`/`cross_diffusion`
  terms, the `r`/`g`/`fw` destruction chain, `cw1`'s derivation):
  Spalart, P.R. & Allmaras, S.R. (1992), "A One-Equation Turbulence Model for
  Aerodynamic Flows," AIAA Paper 92-0439, 30th Aerospace Sciences Meeting,
  Reno, NV, January 1992. Also published as: Spalart, P.R. & Allmaras, S.R.
  (1994), "A One-Equation Turbulence Model for Aerodynamic Flows," *La
  Recherche Aérospatiale*, No. 1, pp. 5-21 (a journal reprint of the same
  paper, not a distinct/later revision). This project implements the
  "SA-noft2" variant (no trip term `ft2`, fully-turbulent flow assumed
  everywhere — a deliberate simplification, see
  [docs/archive/rans-spalart-allmaras-tracker.md](archive/rans-spalart-allmaras-tracker.md)'s
  architecture decision).
- **Negative-`S~` robustness fix** (`cv2=0.7`, `cv3=0.9`, the piecewise `s_tilde`
  formula in `compute_sa_source_terms()`): Allmaras, S.R., Johnson, F.T., &
  Spalart, P.R. (2012), "Modifications and Clarifications for the
  Implementation of the Spalart-Allmaras Turbulence Model," 7th International
  Conference on Computational Fluid Dynamics (ICCFD7), Big Island, Hawaii,
  ICCFD7-1902. **Correction**: this fix was previously attributed in this
  codebase's own comments (`SpalartAllmaras.h`/`.cpp`,
  `docs/archive/rans-spalart-allmaras-tracker.md`) to "Spalart & Allmaras
  1994" — verified against NASA's Turbulence Modeling Resource and the 2012
  paper itself to be incorrect; the exact `cv2`/`cv3` formulation is specific
  to the 2012 paper, not the original 1992/1994 model definition. Fixed at
  all three sites as part of compiling this document.

**Validation**: `--verify-sa-source` checks `compute_sa_source_terms()` in
isolation — a zero-vorticity/`nut=0` state is an exact fixed point (all three
source terms exactly 0, confirmed over 1000 steps), and a manufactured smooth
`omega`/`nut` field produces finite, non-negative source terms everywhere.
`--verify-rans-stability` couples this into a real `RANSTurbulenceSASolver` run on a
sheared/turbulent Couette-like setup and confirms both the mean-flow and
`nut` residuals settle rather than diverge or persistently oscillate
(measured: `rho_u` residual 3.6e-12 -> 2.4e-16, `nut` residual 1.3e-8 ->
3.0e-13 from step 10,000 to the final step). `--verify-flat-plate` is the
closest thing to a literature comparison for the coupled solver, but compares
against the **laminar** Pohlhausen profile (see
[Boundary-layer / wall diagnostics](#boundary-layer--wall-diagnostics)
below), not a turbulent SA reference — **SA's actual turbulence-sustaining
behavior has never been validated against a real turbulent flow in this
project**; see that tracker's Phase 4 for the full disclosure. This remains
the single largest open validation gap in the RANS solver.

## k-omega SST turbulence model

Implemented in [include/SSTKOmega.h](../include/SSTKOmega.h) /
[src/SSTKOmega.cpp](../src/SSTKOmega.cpp), coupled into a real solver in
[RANSTurbulenceSSTSolver](../include/RANSTurbulenceSSTSolver.h). A second,
independent RANS closure alongside Spalart-Allmaras above, not a revision of
it.

- **Base model** (production/destruction/cross-diffusion terms, `F1`/`F2`
  blending functions, the eddy-viscosity limiter): Menter, F.R. (1994),
  "Two-Equation Eddy-Viscosity Turbulence Models for Engineering
  Applications," *AIAA Journal*, 32(8), pp. 1598-1605. Revised constant set:
  Menter, F.R., Kuntz, M., Langtry, R. (2003), "Ten Years of Industrial
  Experience with the SST Turbulence Model," *Turbulence, Heat and Mass
  Transfer 4*.
- **Exact numeric constants** (`beta_star=0.09`, `kappa=0.41`, `a1=0.31`,
  `sigma_k1=0.85`, `sigma_k2=1.0`, `sigma_omega1=0.5`, `sigma_omega2=0.856`,
  `beta1=0.075`, `beta2=0.0828`; `gamma1`/`gamma2` derived, not independent):
  NASA's Turbulence Modeling Resource, currently at
  `tmbwg.github.io/turbmodels/sst.html` (the historical
  `turbmodels.larc.nasa.gov` URL 301-redirects there via `nasa.gov`) —
  chosen as the single source for every constant specifically because
  published SST implementations disagree at the 3rd-4th decimal depending on
  which paper/revision they trace to (e.g. `gamma1` derived as `~0.553` here,
  vs. `5/9 ~ 0.556` some codes hardcode instead) — see
  [docs/sst-komega-tracker.md](sst-komega-tracker.md)'s "On literature
  constants" section for the full reasoning against mixing sources.
- **Eddy-viscosity limiter variant**: BOTH documented variants are
  implemented and user-selectable (`SSTLimiterVariant::Vorticity` = Menter's
  1994 original, using vorticity magnitude `Omega`; `::StrainRate` = the
  2003 revision, using strain-rate magnitude `S`), each paired with its own
  production-clip coefficient (20 vs. 10) — a deliberate choice not to pick
  one, since NASA TMR documents both as distinct verification cases with
  different reference results.
- **Kato-Launder production limiter** (replaces production's `S^2` with
  `S*Omega`, optional, off by default): Kato, M. & Launder, B.E. (1993),
  "The Modelling of Turbulent Flow Around Stationary and Vibrating Square
  Cylinders," Proc. 9th Symposium on Turbulent Shear Flows, Kyoto, Japan,
  August 1993, pp. 10.4.1-10.4.6.
- **Wall boundary condition** (`omega_wall = 60*nu/(beta1*d1^2)`): the
  standard near-wall analytic asymptotic value for `omega` in the log/viscous
  sublayer, following the same convention Menter's own papers and Wilcox,
  D.C. (2006), *Turbulence Modeling for CFD*, 3rd ed., use — no single
  additional citation beyond the base-model papers above.

**Validation**: `--verify-sst-source` checks `compute_sst_source_terms()` in
isolation — a zero-`k`/`omega` state is an exact fixed point; a manufactured
linear mean-flow velocity field (whose gradient `GradientCalculator`
reconstructs exactly) gives machine-precision-exact strain-rate/vorticity
magnitudes to check against; `F1`/`F2` are checked against their known
near-wall (`->1`)/far-field (`->0`) asymptotic limits; production/destruction/
cross-diffusion are checked both for sane finite/non-negative values and via
an independent re-aggregation from the returned intermediates (exact, `0`
error) — for both `SSTLimiterVariant` options.
`--verify-sst-stability` couples this into a real `RANSTurbulenceSSTSolver`
run on the identical sheared/turbulent Couette-like setup `--verify-rans-stability`
uses, for both `SSTLimiterVariant` options plus a third scenario exercising a
live `Farfield` `k`/`omega` boundary, and confirms every residual settles
rather than diverges. `--verify-sst-flat-plate` mirrors `--verify-flat-plate`'s
exact setup (same mesh/`Re_L`) but decides its comparison target dynamically
by monitoring `nu_t/nu`'s own trend during the run rather than assuming the
outcome: it was observed monotonically decaying (`1.6e-5` -> `1.4e-6` over
20000 steps), which routed the code to the same laminar Pohlhausen comparison
SA uses (L2 error `7.5%`, matching SA's own `7.5%` on the identical mesh/BCs —
corroborating evidence this is a property of the tractable Reynolds-number
regime, not an SST-specific defect) — **SST's actual turbulence-sustaining
behavior has never been validated against a real turbulent flow in this
project either**; see docs/sst-komega-tracker.md's Phase 4 for the full
disclosure. Same open validation gap as Spalart-Allmaras, for the same
underlying (acoustic-CFL time-stepping) reason.

## Numerical flux schemes (Euler / Navier-Stokes / RANS inviscid term)

Implemented in [include/RusanovFlux.h](../include/RusanovFlux.h),
[include/HllcFlux.h](../include/HllcFlux.h),
[include/ExactRiemannFlux.h](../include/ExactRiemannFlux.h), selected via the
case-file `flux_scheme` key. Full design/validation history in
[docs/hllc-and-exact-riemann-plan.md](hllc-and-exact-riemann-plan.md) — this
section summarizes it; that document has the complete numeric results.

- **Rusanov (local Lax-Friedrichs)**: no external citation in this codebase's
  comments — implemented as the textbook-standard scheme (also covered in
  Toro's book below, among many others); the most robust and cheapest of the
  three, and the most dissipative.
- **HLLC**: Toro, E.F. (2009), *Riemann Solvers and Numerical Methods for
  Fluid Dynamics: A Practical Introduction*, 3rd ed., Springer, Ch. 9 (HLL
  family, wave-speed estimates) and Ch. 10 (HLLC); Toro, E.F., Spruce, M., &
  Speares, W. (1994), "Restoration of the contact surface in the HLL-Riemann
  solver," *Shock Waves*, 4(1), 25-34 (the original HLLC paper); Harten, A.,
  Lax, P.D., & van Leer, B. (1983), "On upstream differencing and Godunov-type
  schemes for hyperbolic conservation laws," *SIAM Review*, 25(1), 35-61 (the
  HLL framework HLLC extends); Davis, S.F. (1988), "Simplified second-order
  Godunov-type methods," *SIAM Journal on Scientific and Statistical
  Computing*, 9(3), 445-473 (the `SL`/`SR` wave-speed estimate actually used).
  Einfeldt, B. (1988), "On Godunov-type methods for gas dynamics," *SIAM
  Journal on Numerical Analysis*, 25(2), 294-318, gives a more refined
  wave-speed estimate (HLLE) — **considered but not implemented**; Davis'
  simpler estimate proved adequate during validation.
- **Exact Riemann solver (Godunov's method)**: Toro (2009) above, Ch. 4.

**Validation** (full numbers in
[docs/hllc-and-exact-riemann-plan.md](hllc-and-exact-riemann-plan.md)'s
"Status" section): consistency (`F*(U,U)=F(U)`) and reflection-symmetry
checks pass at machine precision for all three schemes. Sod's shock tube
(Sod, G.A. (1978), "A survey of several finite difference methods for
systems of nonlinear hyperbolic conservation laws," *Journal of Computational
Physics*, 27(1), 1-31 — the test case itself, already this project's
existing Euler example) against the exact analytic solution: L1 density
error Rusanov 6.59e-3 > HLLC 4.66e-3 > Exact 4.47e-3, the expected ordering;
the exact solver's own analytic reference reproduced textbook Sod values
(`p*=0.30313`, `u*=0.92745`, `rho*_L=0.42632`, `rho*_R=0.26557`) to 5 decimal
places. Contact-discontinuity sharpness, strong-shock/vacuum-forming
robustness, and performance overhead (+9% HLLC, +23% Exact vs. Rusanov,
Release build) were all measured directly, not assumed from the textbook.

## Gradient reconstruction

Implemented in
[include/GradientReconstruction.h](../include/GradientReconstruction.h)/[.cpp](../src/GradientReconstruction.cpp).
The non-orthogonality-corrected face-normal derivative (`face_gradient()`)
implements the standard **over-relaxed decomposition**: Jasak, H. (1996),
"Error Analysis and Estimation for the Finite Volume Method with
Applications to Fluid Flows," PhD Thesis, Imperial College of Science,
Technology and Medicine, University of London (available at
[spiral.imperial.ac.uk/handle/10044/1/8335](https://spiral.imperial.ac.uk/handle/10044/1/8335)).

**Validation**: `--verify-gradient` checks both `GradientCalculator` schemes
against a known linear analytic field on a self-generated skewed mesh —
Least-Squares reproduces the exact gradient near machine precision (max
error ~9.3e-15), while Green-Gauss shows a real ~7.3% error from
non-orthogonality (expected; Green-Gauss is only exact for a linear field
when every face midpoint lies on the line joining its two cell centroids).
The corrected face-normal derivative itself reproduces the linear field's
exact normal derivative to ~1.8e-14, versus a naive uncorrected two-point
difference's 1.74 (absolute) error on the same skewed mesh — the concrete
demonstration of why the Jasak correction is load-bearing here, not
cosmetic.

## Time-step stability estimate

`compute_dt()` in
[NavierStokesFVMSolver](../include/NavierStokesFVMSolver.h)/[RANSTurbulenceSASolver](../include/RANSTurbulenceSASolver.h)/[RANSTurbulenceSSTSolver](../include/RANSTurbulenceSSTSolver.h)
combines the inviscid CFL condition with the standard explicit-viscous-
diffusion stability limit (`dt < length^2/(2*nu)`) — described in this
codebase's comments only as "Blazek-style," referring to: Blazek, J. (2015),
*Computational Fluid Dynamics: Principles and Applications*, 3rd ed.,
Butterworth-Heinemann/Elsevier (1st ed. 2001, Elsevier Science) — a standard
CFD graduate textbook, not a specific paper; no single equation number is
cited in this codebase's comments, so this is a general-methodology
attribution rather than a line-by-line one.
`RANSTurbulenceSSTSolver::compute_dt()` adds one further term with no
external citation at all — an explicit-SOURCE stability limit
(`cfl / (max(beta_star, beta1, beta2) * omega)`) guarding against `omega`'s
own destruction term overshooting in a single step near the wall, distinct
from (and in addition to) the diffusion limit above. This is this project's
own addition, not drawn from Blazek or any other cited source — see
[docs/sst-komega-tracker.md](sst-komega-tracker.md) Phase 2 for the
derivation and the divergence it was found to actually prevent.

**Validation**: this is a stability heuristic, not a result with a
literature-quoted numeric target, so it's validated indirectly by whether
runs stay stable rather than by matching a published number.
`--verify-ns-stretched-cfl` specifically stress-tests it on a genuinely
anisotropic (boundary-layer-clustered) mesh at `cfl=2` and confirms the run
completes without diverging; every other Navier-Stokes/RANS `--verify-*`
gate (`--verify-couette`, `--verify-ns-uniform`, `--verify-rans-stability`,
`--verify-flat-plate`, `--verify-sst-stability`, `--verify-sst-flat-plate`)
also exercises it implicitly on every step.

## Boundary-layer / wall diagnostics

Implemented in [include/WallTraction.h](../include/WallTraction.h)/[.cpp](../src/WallTraction.cpp),
exercised by `run_verify_flat_plate_boundary_layer()`/`run_verify_sst_flat_plate()`
in [src/main.cpp](../src/main.cpp) (`--verify-flat-plate`/`--verify-sst-flat-plate`).

- **Blasius laminar skin friction**, `Cf(x) = 0.664/sqrt(Re_x)`: Blasius, H.
  (1908), "Grenzschichten in Flüssigkeiten mit kleiner Reibung," *Zeitschrift
  für Mathematik und Physik*, 56, 1-37 (English translation: NACA Technical
  Memorandum 1256). The `0.664` coefficient is the standard textbook result
  derived from Blasius's exact similarity solution, not stated in that exact
  form in the 1908 paper itself.
- **Pohlhausen quartic velocity profile** (`u/U = 2*eta - 2*eta^3 + eta^4`)
  and its closed-form displacement/momentum-thickness ratios
  (`delta*/delta_99 = 3/10`, `theta/delta_99 = 37/315`): Pohlhausen, K.
  (1921), "Zur näherungsweisen Integration der Differentialgleichung der
  laminaren Grenzschicht," *Zeitschrift für Angewandte Mathematik und
  Mechanik*, 1(4), 252-290 — the paper that established the Kármán-Pohlhausen
  integral method.
- **Schlichting turbulent flat-plate correlation**, `Cf_x = 0.0592*Re_x^-0.2`
  (the 1/7-power-law, valid `5e5 < Re_x < 1e7`): Schlichting, H. (1979),
  *Boundary-Layer Theory*, 7th ed., McGraw-Hill, ISBN 0-07-055334-3, eq.
  21.12. **Used only as an a priori mesh-sizing heuristic** (e.g. to compute
  a target first-cell height for a given `y+`) — not implemented as a
  runtime formula inside any solver, so it has no corresponding `--verify-*`
  gate; it's design guidance, not code.

**Validation**: `--verify-flat-plate` runs `RANSTurbulenceSASolver` on a boundary-
layer-clustered flat-plate mesh at a sub-transition Reynolds number (where
SA correctly predicts negligible `nu_t`) and compares the resulting laminar
mean-flow profile against Pohlhausen's quartic approximation: L2 error 7.5%
(under a 15% tolerance). A second, independent check reruns the same exit
station through the general `compute_boundary_layer_profiles()`/
`compute_wall_traction()` machinery (not the hand-indexed column walk) and
compares against Blasius's `Cf` (33% relative error, under a looser 40%
tolerance — attributed to this test's domain acting more like a shallow
duct-entrance flow than a true isolated flat plate, a disclosed confinement
effect, not a code defect) and Pohlhausen's exact thickness ratios (within
~2-5%). `--verify-wall-forces` separately validates `WallTraction.h`'s
`tau_wall`/`Cf`/`Cp`/`y+`/moment against the **exact analytic** planar
Couette flow solution (a first-principles result, not literature) to near
machine precision on `tau_wall`/`y+` and 5% on `Cp`/moment (a disclosed,
understood secular-drift artifact of that specific test's boundary
conditions, not a Blasius/Pohlhausen comparison).

## Considered but not implemented

- Einfeldt, B. (1988), "On Godunov-type methods for gas dynamics," *SIAM
  Journal on Numerical Analysis*, 25(2), 294-318 — a more robust HLLC
  wave-speed estimate (HLLE); Davis' simpler estimate (above) proved adequate
  during validation, so this was never implemented.
- Barth, T.J. & Jespersen, D.C. (1989), "The Design and Application of
  Upwind Schemes on Unstructured Meshes," AIAA Paper 89-0366, 27th Aerospace
  Sciences Meeting, Reno, NV — an unstructured-mesh slope limiter, named in
  passing in [docs/dns-higher-order-scheme-plan.md](dns-higher-order-scheme-plan.md)
  as a candidate if this project ever moves beyond its current first-order,
  no-limiter Euler/Navier-Stokes scheme (see
  [Discontinuity handling](MANUAL.md#discontinuity-handling-shock-capturing)).
  Not implemented; no validation exists for it here.

## No external literature — original/in-house implementations

Documented here explicitly so this file is a complete map, not just a list
of the parts that cite something:

- **`WallDistance.h`/`.cpp`** — a brute-force `O(cells * wall_faces)`
  point-to-segment distance, no spatial acceleration structure. This is this
  project's own straightforward geometric implementation, not drawn from a
  specific cited method. Validated via `--verify-wall-distance` against the
  exact analytic distance to a flat horizontal wall (`max |distance -
  y_centroid| = 0`, i.e. exact to floating-point precision on that geometry).
- **`VtkWriter.cpp`** — implements the well-known legacy VTK ASCII
  `DATASET UNSTRUCTURED_GRID` grammar. No specific format specification
  document is cited in-code; readers wanting the canonical grammar reference
  can consult Kitware's *VTK File Formats* document, but that isn't something
  this codebase's comments themselves point to.
- **FVMESH format** ([docs/fvmesh-format.md](fvmesh-format.md)) — this
  project's own original mesh format, by design not based on any external
  specification (see that document's "Non-goals").
- **Checkpoint binary format** ([include/Checkpoint.h](../include/Checkpoint.h)) —
  an in-house fixed-layout format, not modeled on any external
  checkpoint/restart convention.
