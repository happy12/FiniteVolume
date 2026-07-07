// SPDX-License-Identifier: GPL-3.0-only
#ifndef SPALARTALLMARAS_H_INCLUDED
#define SPALARTALLMARAS_H_INCLUDED

#include "GradientReconstruction.h" // Gradient2

// SA-noft2 model constants (no trip term ft2, fully-turbulent flow assumed
// everywhere -- the architecture decision in
// docs/archive/rans-spalart-allmaras-tracker.md), collected into a struct
// (rather than fixed constexpr globals) so a caller (e.g. RANSTurbulenceSASolver, via
// a case-file override) can substitute non-default values -- cw1 is left out
// since it's derived from cb1/kappa/sigma, not independent; see cw1() below.
struct SAModelConstants {
    double cb1 = 0.1355;
    double cb2 = 0.622;
    double sigma = 2.0 / 3.0;
    double kappa = 0.41;
    double cw2 = 0.3;
    double cw3 = 2.0;
    double cv1 = 7.1;
    double cv2 = 0.7; // negative-S~ robustness fix (Allmaras, Johnson & Spalart 2012, ICCFD7-1902) -- see compute_sa_source_terms()
    double cv3 = 0.9;
};

// cw1 = cb1/kappa^2 + (1+cb2)/sigma -- the model's own derivation, not an
// independently tunable constant (Spalart & Allmaras 1992).
double sa_cw1(const SAModelConstants& c);

// The nut ("nu-tilde") transport equation's three source terms at a single
// cell, evaluated from local field values only (see compute_sa_source_terms()
// below). Kept separate rather than pre-summed so a caller/test can inspect
// each independently, matching this project's convention of returning
// structured intermediate results (see GradientReconstruction.h's FaceGradient).
struct SASourceTerms {
    double production = 0.0;      // Cb1 * S_tilde * nut
    double destruction = 0.0;     // Cw1 * fw * (nut/d)^2
    double cross_diffusion = 0.0; // (Cb2/sigma) * |grad(nut)|^2 -- the model's
                                    // non-conservative diffusion term only. The
                                    // conservative div((nu+nu_t)/sigma * grad(nut))
                                    // part needs mesh face geometry and is
                                    // deferred to Phase 3's coupling into an
                                    // actual solver, reusing GradientReconstruction.h's
                                    // already-verified face_gradient()/GradientCalculator
                                    // unchanged -- the same way NavierStokesFVMSolver
                                    // reuses them for the viscous stress tensor.
};

// Computes the SA eddy-viscosity-response function fv1(chi), chi = nut/nu.
// Returns 0 at nut <= 0 (chi <= 0) -- this codebase implements only the
// standard (nut >= 0) SA-noft2 formulation, not the separate "negative SA"
// extension some solvers add for nut < 0.
// 'c' defaults to the standard SA-noft2 constants; pass a caller-configured
// SAModelConstants (e.g. from a case file) to override them.
double sa_fv1(double nut, double nu, const SAModelConstants& c = SAModelConstants{});

// Computes nu_t = nut * fv1(nut, nu), the derived turbulent eddy viscosity --
// the one quantity Phase 3's coupling needs to substitute into
// mu_eff = mu + rho*nu_t (see the tracker's "What SA adds..." section).
double sa_eddy_viscosity(double nut, double nu, const SAModelConstants& c = SAModelConstants{});

// Computes this cell's production/destruction/cross-diffusion source terms.
//
// Methodology (SA-noft2, standard form, with the standard negative-S~
// robustness fix -- Allmaras, Johnson & Spalart 2012, ICCFD7-1902; added in
// Phase 3 of docs/archive/rans-spalart-allmaras-tracker.md after a real
// coupled RANSTurbulenceSASolver run hit exactly the failure mode this fix exists
// for -- see that phase's notes):
//   chi   = nut / nu (0 if nut <= 0)
//   fv1   = chi^3 / (chi^3 + Cv1^3)
//   fv2   = 1 - chi / (1 + chi*fv1)
//   Sbar  = (nut / (kappa^2 * d^2)) * fv2   -- can be negative (fv2 is
//           negative for a wide range of chi, roughly 1 to 18)
//   S~    = omega + Sbar,                          if Sbar >= -Cv2*omega
//         = omega + omega*(Cv2^2*omega + Cv3*Sbar)
//                   / ((Cv3 - 2*Cv2)*omega - Sbar), otherwise
//           (the second form is bounded away from 0 whenever omega > 0, and
//           reduces to exactly 0 when omega == 0 too, avoiding the plain
//           formula's sign flip -- see "Known limitation" below for what
//           this does NOT fully resolve)
//   production = Cb1 * S~ * nut
//   r     = min(nut / (S~ * kappa^2 * d^2), 10)   -- the standard SA clip
//           (Spalart & Allmaras 1992), part of the model's own definition,
//           not a numerical-safety patch; only evaluated when nut > 0 (see
//           "Known limitation" below)
//   g     = r + Cw2*(r^6 - r)
//   fw    = g * ((1 + Cw3^6) / (g^6 + Cw3^6))^(1/6)
//   destruction = Cw1 * fw * (nut/d)^2
//   cross_diffusion = (Cb2/sigma) * |grad(nut)|^2
//
// Known limitation: the fix above guarantees S~ doesn't flip sign or divide
// by zero, but it does NOT make destruction's own r/g/fw chain well-behaved
// at every combination of nut/omega/d -- r's denominator (S~ * kappa^2 * d^2)
// can still be extremely small (though never exactly 0 when omega > 0) when
// omega is tiny and d is small, driving r toward its 10-clip very
// aggressively and destruction toward a very large, stiff value relative to
// this solver's explicit dt. Phase 3's own verification (a moderate,
// literature-typical freestream nut/nu ~ 3) exercises this fix directly and
// stays stable, but an even larger freestream nut/nu, or a much thinner
// near-wall cell than this tracker's test meshes use, has not been tried.
//
// Input:
//   nut           - local nu-tilde value (the transported scalar)
//   nu            - local molecular kinematic viscosity (mu/rho)
//   omega         - local vorticity magnitude, |dv/dx - du/dy| in 2D
//                    (computable directly from NavierStokesFVMSolver's
//                    existing grad_u/grad_v cell gradients -- see the
//                    tracker's "What SA adds..." section)
//   wall_distance - distance to the nearest no-slip wall (see WallDistance.h);
//                    must be > 0 (no guard here -- consistent with this
//                    project's "no mesh-quality safeguards" stance)
//   grad_nut      - reconstructed cell gradient of nut (see GradientReconstruction.h)
//   c             - model constants; defaults to the standard SA-noft2 set
// Returns: the three source terms above. destruction is forced to exactly
//          0.0 when nut <= 0 rather than evaluated through the r/g/fw chain
//          -- both because destruction is physically ~0 there, and because
//          r's denominator (S~ * kappa^2 * d^2) is a literal 0/0 at
//          nut == 0 AND omega == 0 (S~ == 0 too in that state); evaluating
//          it anyway would produce nut * NaN = NaN even though the true
//          limit is exactly 0.
SASourceTerms compute_sa_source_terms(double nut, double nu, double omega, double wall_distance,
                                        const Gradient2& grad_nut, const SAModelConstants& c = SAModelConstants{});

#endif // SPALARTALLMARAS_H_INCLUDED
