// SPDX-License-Identifier: GPL-3.0-only
#include "SpalartAllmaras.h"

#include <algorithm>
#include <cmath>

// See SpalartAllmaras.h for the input/output contract.
double sa_cw1(const SAModelConstants& c) {
    return c.cb1 / (c.kappa * c.kappa) + (1.0 + c.cb2) / c.sigma;
}

// See SpalartAllmaras.h for the input/output contract.
double sa_fv1(double nut, double nu, const SAModelConstants& c) {
    if (nut <= 0.0) return 0.0;
    double chi = nut / nu;
    double chi3 = chi * chi * chi;
    return chi3 / (chi3 + c.cv1 * c.cv1 * c.cv1);
}

// See SpalartAllmaras.h for the input/output contract.
double sa_eddy_viscosity(double nut, double nu, const SAModelConstants& c) {
    return nut * sa_fv1(nut, nu, c);
}

// See SpalartAllmaras.h for the input/output contract and full methodology.
SASourceTerms compute_sa_source_terms(double nut, double nu, double omega, double wall_distance,
                                        const Gradient2& grad_nut, const SAModelConstants& c) {
    SASourceTerms result;

    double chi = (nut > 0.0) ? nut / nu : 0.0;
    double fv1 = sa_fv1(nut, nu, c);
    double fv2 = 1.0 - chi / (1.0 + chi * fv1);

    double d2 = wall_distance * wall_distance;
    double s_bar = (nut / (c.kappa * c.kappa * d2)) * fv2;

    // Negative-S~ robustness fix (Spalart & Allmaras 1994): s_bar is
    // negative for a wide range of chi (fv2 itself is negative roughly for
    // chi in [1, 18]), and the plain omega + s_bar formula can drive S~
    // negative or even flip production's sign when omega is small -- see
    // SpalartAllmaras.h's methodology comment for the exact case this was
    // added to fix (a real divergence in a coupled RANSFVMSolver run).
    double s_tilde;
    if (s_bar >= -c.cv2 * omega) {
        s_tilde = omega + s_bar;
    } else {
        s_tilde = omega + omega * (c.cv2 * c.cv2 * omega + c.cv3 * s_bar) / ((c.cv3 - 2.0 * c.cv2) * omega - s_bar);
    }

    result.production = c.cb1 * s_tilde * nut;

    if (nut > 0.0) {
        // s_tilde_safe floors the denominator, not nut, since destruction's
        // true physical value in the nut>0/S~->0 limit is governed by nut,
        // not by this division -- see SpalartAllmaras.h's note on why r's
        // formula is a literal 0/0 at nut == 0 AND omega == 0.
        double s_tilde_safe = std::max(s_tilde, 1e-10);
        double r = std::min(nut / (s_tilde_safe * c.kappa * c.kappa * d2), 10.0);
        double g = r + c.cw2 * (std::pow(r, 6.0) - r);
        double fw = g * std::pow((1.0 + std::pow(c.cw3, 6.0)) / (std::pow(g, 6.0) + std::pow(c.cw3, 6.0)), 1.0 / 6.0);
        double nut_over_d = nut / wall_distance;
        result.destruction = sa_cw1(c) * fw * nut_over_d * nut_over_d;
    }

    result.cross_diffusion =
        (c.cb2 / c.sigma) * (grad_nut.dphidx * grad_nut.dphidx + grad_nut.dphidy * grad_nut.dphidy);

    return result;
}
