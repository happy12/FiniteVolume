// SPDX-License-Identifier: GPL-3.0-only
#include "SSTKOmega.h"

#include <algorithm>
#include <cmath>

// See SSTKOmega.h for the input/output contract.
double sst_gamma1(const SSTModelConstants& c) {
    return c.beta1 / c.beta_star - c.sigma_omega1 * c.kappa * c.kappa / std::sqrt(c.beta_star);
}

// See SSTKOmega.h for the input/output contract.
double sst_gamma2(const SSTModelConstants& c) {
    return c.beta2 / c.beta_star - c.sigma_omega2 * c.kappa * c.kappa / std::sqrt(c.beta_star);
}

// See SSTKOmega.h for the input/output contract.
double sst_strain_rate_magnitude(const Gradient2& grad_u, const Gradient2& grad_v) {
    double s11 = grad_u.dphidx;
    double s22 = grad_v.dphidy;
    double s12 = 0.5 * (grad_u.dphidy + grad_v.dphidx);
    return std::sqrt(2.0 * (s11 * s11 + s22 * s22 + 2.0 * s12 * s12));
}

// See SSTKOmega.h for the input/output contract.
double sst_vorticity_magnitude(const Gradient2& grad_u, const Gradient2& grad_v) {
    return std::fabs(grad_v.dphidx - grad_u.dphidy);
}

// See SSTKOmega.h for the input/output contract and full methodology.
double sst_F1(double rho, double k, double omega, double nu, double wall_distance, const Gradient2& grad_k,
              const Gradient2& grad_omega, const SSTModelConstants& c) {
    if (k <= 0.0 || omega <= 0.0) return 0.0;

    double d = wall_distance;
    double d2 = d * d;
    double grad_dot = grad_k.dphidx * grad_omega.dphidx + grad_k.dphidy * grad_omega.dphidy;
    double CDkomega = std::max(2.0 * rho * c.sigma_omega2 * grad_dot / omega, 1e-20);

    double term1 = std::max(std::sqrt(k) / (c.beta_star * omega * d), 500.0 * nu / (d2 * omega));
    double term2 = 4.0 * rho * c.sigma_omega2 * k / (CDkomega * d2);
    double arg1 = std::min(term1, term2);

    double arg1_4 = arg1 * arg1 * arg1 * arg1;
    return std::tanh(arg1_4);
}

// See SSTKOmega.h for the input/output contract and full methodology.
double sst_F2(double k, double omega, double nu, double wall_distance, const SSTModelConstants& c) {
    if (k <= 0.0 || omega <= 0.0) return 0.0;

    double d = wall_distance;
    double d2 = d * d;
    double arg2 = std::max(2.0 * std::sqrt(k) / (c.beta_star * omega * d), 500.0 * nu / (d2 * omega));
    return std::tanh(arg2 * arg2);
}

// See SSTKOmega.h for the input/output contract and full methodology.
double sst_eddy_viscosity(double k, double omega, double S, double Omega, double F2, SSTLimiterVariant variant,
                           const SSTModelConstants& c) {
    if (k <= 0.0) return 0.0;
    double limiter = (variant == SSTLimiterVariant::Vorticity) ? Omega : S;
    return c.a1 * k / std::max(c.a1 * omega, limiter * F2);
}

// See SSTKOmega.h for the input/output contract and full methodology.
SSTSourceTerms compute_sst_source_terms(double rho, double k, double omega, double nu, double wall_distance,
                                          const Gradient2& grad_u, const Gradient2& grad_v, const Gradient2& grad_k,
                                          const Gradient2& grad_omega, SSTLimiterVariant variant,
                                          const SSTModelConstants& c, bool kato_launder) {
    SSTSourceTerms result;

    result.S = sst_strain_rate_magnitude(grad_u, grad_v);
    result.Omega = sst_vorticity_magnitude(grad_u, grad_v);

    if (k <= 0.0 || omega <= 0.0) return result;

    result.F1 = sst_F1(rho, k, omega, nu, wall_distance, grad_k, grad_omega, c);
    result.F2 = sst_F2(k, omega, nu, wall_distance, c);
    result.nu_t = sst_eddy_viscosity(k, omega, result.S, result.Omega, result.F2, variant, c);

    // Kato-Launder (1993): replaces production's S^2 with S*Omega -- see
    // SSTKOmega.h's own note on this parameter.
    double production_rate_sq = kato_launder ? (result.S * result.Omega) : (result.S * result.S);

    double mu_t = rho * result.nu_t;
    double beta_star_rho_omega_k = c.beta_star * rho * omega * k;
    double clip_coef = (variant == SSTLimiterVariant::Vorticity) ? 20.0 : 10.0;
    result.k_production = std::min(mu_t * production_rate_sq, clip_coef * beta_star_rho_omega_k);
    result.k_destruction = beta_star_rho_omega_k;

    double gamma = result.F1 * sst_gamma1(c) + (1.0 - result.F1) * sst_gamma2(c);
    double beta = result.F1 * c.beta1 + (1.0 - result.F1) * c.beta2;
    result.omega_production = gamma * rho * production_rate_sq;
    result.omega_destruction = beta * rho * omega * omega;

    double grad_dot = grad_k.dphidx * grad_omega.dphidx + grad_k.dphidy * grad_omega.dphidy;
    result.cross_diffusion = 2.0 * (1.0 - result.F1) * rho * c.sigma_omega2 * grad_dot / omega;

    return result;
}
