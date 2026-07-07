// SPDX-License-Identifier: GPL-3.0-only
#ifndef SSTKOMEGA_H_INCLUDED
#define SSTKOMEGA_H_INCLUDED

#include "GradientReconstruction.h" // Gradient2

// Menter k-omega SST model constants, values from NASA's Turbulence Modeling
// Resource (currently at tmbwg.github.io/turbmodels/sst.html -- the historical
// turbmodels.larc.nasa.gov URL now redirects there via nasa.gov; see
// docs/sst-komega-tracker.md's "On literature constants" section for why this
// single source was chosen over Menter (1994) or Menter/Kuntz/Langtry (2003)
// directly). gamma1/gamma2 are derived from beta1/beta2/beta_star/sigma_omega1/
// sigma_omega2/kappa (see sst_gamma1()/sst_gamma2() below), not independently
// tunable -- same rationale as SpalartAllmaras.h's sa_cw1(). Collected into a
// struct rather than fixed globals so a caller (e.g. RANSTurbulenceSSTSolver, via
// a future case-file override) can substitute non-default values, matching
// SAModelConstants's precedent.
struct SSTModelConstants {
    double beta_star = 0.09;
    double kappa = 0.41;
    double a1 = 0.31;
    double sigma_k1 = 0.85;
    double sigma_k2 = 1.0;
    double sigma_omega1 = 0.5;
    double sigma_omega2 = 0.856;
    double beta1 = 0.075;
    double beta2 = 0.0828;
};

// Which eddy-viscosity limiter (and paired k-production clip coefficient)
// nu_t/production use -- see the tracker's "On literature constants" section:
// NASA TMR documents both as distinct verification cases with different
// reference results, and this implementation exposes both rather than
// picking one, per Mathieu's 2026-07-06 direction. Vorticity is Menter's
// original (1994) SST; StrainRate is the Menter/Kuntz/Langtry (2003)
// revision. Each variant's production-clip coefficient (20 vs. 10, both
// multiplying beta_star*rho*omega*k) is paired with its limiter choice below
// rather than independently configurable, so this never silently mixes a
// clip coefficient from one revision with a limiter from the other -- the
// same "do not mix constants from different sources" discipline the tracker
// applies to the rest of this model.
enum class SSTLimiterVariant {
    Vorticity,  // original SST (Menter 1994): nu_t = a1*k / max(a1*omega, Omega*F2), Pk clipped at 20*beta_star*rho*omega*k
    StrainRate  // SST-2003 (Menter/Kuntz/Langtry): nu_t = a1*k / max(a1*omega, S*F2), Pk clipped at 10*beta_star*rho*omega*k
};

// gamma1 = beta1/beta_star - sigma_omega1*kappa^2/sqrt(beta_star) (Menter's own
// derivation -- NASA TMR reports this as ~0.553, NOT the 5/9 ~ 0.556 some
// codes hardcode instead; see the tracker's "verify before coding" note).
double sst_gamma1(const SSTModelConstants& c);

// gamma2 = beta2/beta_star - sigma_omega2*kappa^2/sqrt(beta_star) (~0.44 per NASA TMR).
double sst_gamma2(const SSTModelConstants& c);

// 2D strain-rate magnitude S = sqrt(2*Sij*Sij), Sij = 0.5*(dui/dxj + duj/dxi):
//   S = sqrt(2*(du/dx)^2 + 2*(dv/dy)^2 + (du/dy + dv/dx)^2)
// Kept as a separate, directly callable function (not fused into any source-
// term computation) per the tracker's "Why SST matters for later work" note:
// a future gamma-Re_theta transition model needs S available on its own.
double sst_strain_rate_magnitude(const Gradient2& grad_u, const Gradient2& grad_v);

// 2D vorticity magnitude Omega = |dv/dx - du/dy| -- identical definition and
// rationale to RANSTurbulenceSASolver's own "omega" input (see
// SpalartAllmaras.h), kept as its own function here for the same
// keep-it-separate reason as sst_strain_rate_magnitude() above.
double sst_vorticity_magnitude(const Gradient2& grad_u, const Gradient2& grad_v);

// F1 blending function (inner k-omega region -> 1, outer k-epsilon-equivalent
// region -> 0), per NASA TMR:
//   CDkomega = max(2*rho*sigma_omega2*(grad_k . grad_omega)/omega, 1e-20)
//   arg1 = min(max(sqrt(k)/(beta_star*omega*d), 500*nu/(d^2*omega)),
//              4*rho*sigma_omega2*k/(CDkomega*d^2))
//   F1 = tanh(arg1^4)
// Returns 0 when k <= 0 or omega <= 0 (both terms feeding arg1's first max()
// are then singular or physically meaningless) -- mirrors compute_sa_source_terms()'s
// convention of a defined, safe value at a degenerate input rather than NaN.
double sst_F1(double rho, double k, double omega, double nu, double wall_distance, const Gradient2& grad_k,
              const Gradient2& grad_omega, const SSTModelConstants& c = SSTModelConstants{});

// F2 blending function (used only by the eddy-viscosity limiter, not by any
// constant blending), per NASA TMR:
//   arg2 = max(2*sqrt(k)/(beta_star*omega*d), 500*nu/(d^2*omega))
//   F2 = tanh(arg2^2)
// Returns 0 when k <= 0 or omega <= 0, same rationale as sst_F1() above.
double sst_F2(double k, double omega, double nu, double wall_distance, const SSTModelConstants& c = SSTModelConstants{});

// Eddy viscosity nu_t = a1*k / max(a1*omega, limiter*F2), limiter = Omega or S
// per 'variant' (see SSTLimiterVariant above). Returns 0 when k <= 0.
double sst_eddy_viscosity(double k, double omega, double S, double Omega, double F2, SSTLimiterVariant variant,
                           const SSTModelConstants& c = SSTModelConstants{});

// This cell's k/omega transport source terms, plus the intermediate
// quantities (S, Omega, F1, F2, nu_t) the tracker's "Why SST matters for
// later work" section asks to keep separately inspectable for a future
// gamma-Re_theta extension. Kept in one struct rather than pre-summed into
// two totals, matching SASourceTerms's precedent.
//
// Unlike SASourceTerms, these terms are rho-weighted (rho appears explicitly
// in every term below) rather than following SA's convention of omitting
// rho entirely -- SA's nut is not a mass-specific quantity in its classical
// (Spalart & Allmaras 1992) form, but k (turbulent kinetic energy per unit
// mass) and omega (specific dissipation rate) are conventionally transported
// in Menter's compressible rho-weighted form (see the NASA TMR equations
// this module follows); this is a genuine difference in the two models'
// standard forms, not an inconsistency to reconcile.
//
// The k-equation's own production term structurally requires the eddy
// viscosity nu_t (Pk = mu_t*S^2, mu_t = rho*nu_t) even though nu_t is
// nominally introduced in the tracker's Phase 2 ("Eddy-viscosity closure").
// This module computes it here regardless, since Pk cannot be evaluated
// without it; Phase 2's remaining job is coupling this SAME nu_t formula
// into the mean-flow viscous flux and compute_dt(), not deriving nu_t for
// the first time. See docs/sst-komega-tracker.md's Phase 1 section for this
// resolution recorded in one place.
struct SSTSourceTerms {
    double k_production = 0.0;      // min(mu_t*S^2, clip_coef*beta_star*rho*omega*k) -- clip_coef per SSTLimiterVariant
    double k_destruction = 0.0;     // beta_star*rho*omega*k
    double omega_production = 0.0;  // gamma*rho*S^2 (gamma blended via F1) -- equals (gamma/nu_t)*Pk_raw
                                     // simplified via mu_t = rho*nu_t; NASA TMR's own (gamma/nu_t)*P form,
                                     // algebraically independent of any clipping applied to Pk in the k-equation
    double omega_destruction = 0.0; // beta*rho*omega^2 (beta blended via F1)
    double cross_diffusion = 0.0;   // 2*(1-F1)*rho*sigma_omega2*(grad_k . grad_omega)/omega -- omega-equation-only term
    double nu_t = 0.0;               // derived eddy viscosity (kinematic), see sst_eddy_viscosity() above
    double F1 = 0.0;
    double F2 = 0.0;
    double S = 0.0;     // strain-rate magnitude, see sst_strain_rate_magnitude() above
    double Omega = 0.0; // vorticity magnitude, see sst_vorticity_magnitude() above
};

// Computes this cell's k/omega source terms and intermediate quantities.
//
// Input:
//   rho           - local density, mesh-consistent units
//   k, omega      - local transported k/omega values
//   nu            - local molecular kinematic viscosity (mu/rho)
//   wall_distance - distance to the nearest no-slip wall (WallDistance.h);
//                    must be > 0 (no guard here, consistent with this
//                    project's "no mesh-quality safeguards" stance)
//   grad_u, grad_v - reconstructed cell gradients of the mean-flow velocity
//                    components (GradientReconstruction.h), used to derive
//                    both S and Omega
//   grad_k, grad_omega - reconstructed cell gradients of k/omega
//   variant       - which eddy-viscosity limiter/production-clip pairing to use
//   c             - model constants; defaults to the NASA TMR set above
//   kato_launder  - if true, replaces the raw production P = mu_t*S^2 with
//                    the Kato-Launder form P = mu_t*S*Omega in BOTH
//                    k_production (before its clip) and omega_production
//                    (which is P's own gamma*rho*S^2 simplification,
//                    replaced consistently with gamma*rho*S*Omega) --
//                    docs/sst-komega-tracker.md Phase 3's stagnation-point
//                    over-production limiter, SST-specific, no SA analogue.
//                    S*Omega == S^2 in pure shear (where production is
//                    already correct), so this only changes anything in
//                    irrotational/stagnation-point regions where Omega << S;
//                    defaults to false (this tracker's Phase 1/2 gates never
//                    enabled it, and every existing call site keeps their
//                    prior behavior unchanged).
// Returns: the source terms and intermediates above, all forced to exactly
//          0.0 (rather than evaluated through F1/F2's singular-at-omega==0
//          formulas) when k <= 0 or omega <= 0 -- both because every term is
//          physically ~0 in that state, and because F1/F2/nu_t's formulas
//          divide by omega, a literal division-by-zero at omega == 0.
SSTSourceTerms compute_sst_source_terms(double rho, double k, double omega, double nu, double wall_distance,
                                          const Gradient2& grad_u, const Gradient2& grad_v, const Gradient2& grad_k,
                                          const Gradient2& grad_omega, SSTLimiterVariant variant,
                                          const SSTModelConstants& c = SSTModelConstants{}, bool kato_launder = false);

#endif // SSTKOMEGA_H_INCLUDED
