// SPDX-License-Identifier: GPL-3.0-only
#ifndef EXACTRIEMANNFLUX_H_INCLUDED
#define EXACTRIEMANNFLUX_H_INCLUDED

#include <cmath>
#include <algorithm>

#include "EulerState.h"

// Exact Riemann (Godunov) solver: one of the numerical fluxes selectable via
// EulerFVMSolver's NumericalFluxScheme (see EulerFVMSolver.h). This is the
// exact solution to the local 1D Riemann problem rather than an
// approximation -- highest possible accuracy for a first-order scheme, at
// the cost of a per-face Newton-Raphson iteration (see the cost caveat in
// docs/hllc-and-exact-riemann-plan.md). See Toro, "Riemann Solvers and
// Numerical Methods for Fluid Dynamics" 3rd ed., Ch. 4 -- equation numbers
// referenced below are from there.
//
// Kept separate from EulerState.h (state type + underlying physics only)
// and from RusanovFlux.h/HllcFlux.h (the other numerical schemes), since
// this scheme's internals (Newton-Raphson pressure solve, wave-pattern
// sampling, vacuum handling) are substantially larger and self-contained.
namespace ExactRiemannDetail {

// The pressure function f_K(p) (Toro eq. 4.6/4.7): the velocity jump across
// wave K (left or right) implied by a trial star pressure p, either a shock
// relation (p > p_K) or an isentropic rarefaction relation (p <= p_K).
inline double pressure_function(double p, double rho_K, double p_K, double c_K, double gamma) {
    if (p > p_K) {
        double A_K = 2.0 / ((gamma + 1.0) * rho_K);
        double B_K = (gamma - 1.0) / (gamma + 1.0) * p_K;
        return (p - p_K) * std::sqrt(A_K / (p + B_K));
    }
    return (2.0 * c_K / (gamma - 1.0)) * (std::pow(p / p_K, (gamma - 1.0) / (2.0 * gamma)) - 1.0);
}

// d(pressure_function)/dp, used by the Newton-Raphson iteration.
inline double pressure_function_deriv(double p, double rho_K, double p_K, double c_K, double gamma) {
    if (p > p_K) {
        double A_K = 2.0 / ((gamma + 1.0) * rho_K);
        double B_K = (gamma - 1.0) / (gamma + 1.0) * p_K;
        return std::sqrt(A_K / (p + B_K)) * (1.0 - 0.5 * (p - p_K) / (p + B_K));
    }
    return (1.0 / (rho_K * c_K)) * std::pow(p / p_K, -(gamma + 1.0) / (2.0 * gamma));
}

// Initial guess for the Newton-Raphson star-pressure solve: the two-shock
// estimate (Toro eq. 4.47) seeded by the primitive-variable estimate p_pv
// (Toro eq. 4.46), floored above zero since a trial pressure must be
// positive. Chosen over a plain PVRS guess for faster/more robust
// convergence on strong waves.
inline double guess_pressure(double rho_L, double u_L, double p_L, double c_L,
                              double rho_R, double u_R, double p_R, double c_R, double gamma) {
    double p_pv = std::max(1e-6, 0.5 * (p_L + p_R) - 0.125 * (u_R - u_L) * (rho_L + rho_R) * (c_L + c_R));
    double A_L = 2.0 / ((gamma + 1.0) * rho_L), B_L = (gamma - 1.0) / (gamma + 1.0) * p_L;
    double A_R = 2.0 / ((gamma + 1.0) * rho_R), B_R = (gamma - 1.0) / (gamma + 1.0) * p_R;
    double g_L = std::sqrt(A_L / (p_pv + B_L));
    double g_R = std::sqrt(A_R / (p_pv + B_R));
    return std::max(1e-6, (g_L * p_L + g_R * p_R - (u_R - u_L)) / (g_L + g_R));
}

// Newton-Raphson solve for the star-region pressure p* (Toro eq. 4.85/4.86):
// the root of f_L(p) + f_R(p) + (u_R - u_L) = 0. rel_tol/max_iter are
// user-facing (see NumericalFluxScheme/EulerFVMSolver.h's
// exact_riemann_tol/exact_riemann_max_iter, and CaseInput.h's
// exact_riemann_tol/exact_riemann_max_iter case-file keys) but are not
// meant to be tuned in normal use -- they're exposed for transparency,
// with the same default values (1e-6, 20) this solver always used before
// they became configurable.
// Input:   rho_L/u_L/p_L/c_L, rho_R/u_R/p_R/c_R - primitive left/right
//          states (normal-frame velocity) and their sound speeds; gamma
//          rel_tol   - relative pressure change that stops the iteration
//          max_iter  - iteration cap, in case rel_tol is never reached
// Returns: p*, the pressure in the star region (uniform across both star
//          states, by construction of the Riemann problem)
inline double solve_star_pressure(double rho_L, double u_L, double p_L, double c_L,
                                   double rho_R, double u_R, double p_R, double c_R, double gamma,
                                   double rel_tol, int max_iter) {
    double p = guess_pressure(rho_L, u_L, p_L, c_L, rho_R, u_R, p_R, c_R, gamma);
    for (int iter = 0; iter < max_iter; ++iter) {
        double f = pressure_function(p, rho_L, p_L, c_L, gamma)
                 + pressure_function(p, rho_R, p_R, c_R, gamma)
                 + (u_R - u_L);
        double fd = pressure_function_deriv(p, rho_L, p_L, c_L, gamma)
                  + pressure_function_deriv(p, rho_R, p_R, c_R, gamma);
        double p_new = p - f / fd;
        if (p_new < 1e-10) p_new = 1e-10; // trial pressure must stay positive
        double change = 2.0 * std::fabs(p_new - p) / (p_new + p);
        p = p_new;
        if (change < rel_tol) break;
    }
    return p;
}

// Samples the state inside the left rarefaction fan at speed S (Toro
// eq. 4.56), i.e. between the fan's head (u_L - c_L) and tail.
inline void left_fan_interior(double rho_L, double u_L, double p_L, double c_L, double gamma, double S,
                               double& rho_out, double& u_out, double& p_out) {
    double c_fan = (2.0 / (gamma + 1.0)) * (c_L + 0.5 * (gamma - 1.0) * (u_L - S));
    u_out = (2.0 / (gamma + 1.0)) * (c_L + 0.5 * (gamma - 1.0) * u_L + S);
    rho_out = rho_L * std::pow(c_fan / c_L, 2.0 / (gamma - 1.0));
    p_out = p_L * std::pow(c_fan / c_L, 2.0 * gamma / (gamma - 1.0));
}

// Mirror of left_fan_interior for the right rarefaction fan (Toro eq. 4.63).
inline void right_fan_interior(double rho_R, double u_R, double p_R, double c_R, double gamma, double S,
                                double& rho_out, double& u_out, double& p_out) {
    double c_fan = (2.0 / (gamma + 1.0)) * (c_R - 0.5 * (gamma - 1.0) * (u_R - S));
    u_out = (2.0 / (gamma + 1.0)) * (-c_R + 0.5 * (gamma - 1.0) * u_R + S);
    rho_out = rho_R * std::pow(c_fan / c_R, 2.0 / (gamma - 1.0));
    p_out = p_R * std::pow(c_fan / c_R, 2.0 * gamma / (gamma - 1.0));
}

// Samples the solution at speed S on the left of the contact (S <= u*),
// across whichever of the left shock (Toro eq. 4.50/4.52) or left
// rarefaction fan (Toro eq. 4.55/4.56) the star pressure p* implies.
inline void sample_left(double rho_L, double u_L, double p_L, double c_L, double gamma,
                         double p_star, double u_star, double S,
                         double& rho_out, double& u_out, double& p_out) {
    if (p_star > p_L) {
        double S_L = u_L - c_L * std::sqrt((gamma + 1.0) / (2.0 * gamma) * (p_star / p_L) + (gamma - 1.0) / (2.0 * gamma));
        if (S <= S_L) { rho_out = rho_L; u_out = u_L; p_out = p_L; return; }
        double ratio = p_star / p_L;
        rho_out = rho_L * (ratio + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * ratio + 1.0);
        u_out = u_star; p_out = p_star;
    } else {
        double c_star_L = c_L * std::pow(p_star / p_L, (gamma - 1.0) / (2.0 * gamma));
        double S_HL = u_L - c_L, S_TL = u_star - c_star_L;
        if (S <= S_HL) { rho_out = rho_L; u_out = u_L; p_out = p_L; }
        else if (S >= S_TL) { rho_out = rho_L * std::pow(p_star / p_L, 1.0 / gamma); u_out = u_star; p_out = p_star; }
        else left_fan_interior(rho_L, u_L, p_L, c_L, gamma, S, rho_out, u_out, p_out);
    }
}

// Mirror of sample_left for the right of the contact (S >= u*).
inline void sample_right(double rho_R, double u_R, double p_R, double c_R, double gamma,
                          double p_star, double u_star, double S,
                          double& rho_out, double& u_out, double& p_out) {
    if (p_star > p_R) {
        double S_R = u_R + c_R * std::sqrt((gamma + 1.0) / (2.0 * gamma) * (p_star / p_R) + (gamma - 1.0) / (2.0 * gamma));
        if (S >= S_R) { rho_out = rho_R; u_out = u_R; p_out = p_R; return; }
        double ratio = p_star / p_R;
        rho_out = rho_R * (ratio + (gamma - 1.0) / (gamma + 1.0)) / ((gamma - 1.0) / (gamma + 1.0) * ratio + 1.0);
        u_out = u_star; p_out = p_star;
    } else {
        double c_star_R = c_R * std::pow(p_star / p_R, (gamma - 1.0) / (2.0 * gamma));
        double S_HR = u_R + c_R, S_TR = u_star + c_star_R;
        if (S >= S_HR) { rho_out = rho_R; u_out = u_R; p_out = p_R; }
        else if (S <= S_TR) { rho_out = rho_R * std::pow(p_star / p_R, 1.0 / gamma); u_out = u_star; p_out = p_star; }
        else right_fan_interior(rho_R, u_R, p_R, c_R, gamma, S, rho_out, u_out, p_out);
    }
}

} // namespace ExactRiemannDetail

// Exact Riemann (Godunov) numerical flux: solves the local 1D Riemann
// problem exactly (in the face-normal/tangential frame, same rotation
// technique as HLLC -- see HllcFlux.h) and samples it at x/t = 0 (a
// stationary face), which *is* Godunov's flux by definition once that
// sampled state's physical flux is taken.
//
// Methodology:
//   1. Vacuum check (Toro Ch. 4.7): if the two states are separating faster
//      than 2*(c_L+c_R)/(gamma-1) <= u_R-u_L, a vacuum forms between them --
//      sampled explicitly via the closed-form vacuum-front speeds
//      (u_L + 2*c_L/(gamma-1), u_R - 2*c_R/(gamma-1)) rather than running
//      the (inapplicable) Newton iteration below.
//   2. Otherwise, Newton-Raphson for the star pressure p*
//      (solve_star_pressure), then the star velocity u* (Toro eq. 4.9):
//      u* = 0.5*(u_L+u_R) + 0.5*(f_R(p*) - f_L(p*)).
//   3. Sample the solution at S = x/t = 0: left of the contact (S <= u*)
//      via sample_left, right of it via sample_right. Tangential velocity
//      is carried through unchanged from whichever side (L/R) the sampled
//      point falls on, since a contact wave doesn't affect it.
//   4. The physical flux of that sampled state is the numerical flux.
//
// Input:   U_L, U_R  - conserved states on the left/right of the face;
//          nx, ny    - unit normal pointing from left to right; gamma
//          rel_tol, max_iter - see solve_star_pressure above; user-facing
//                      via EulerFVMSolver's exact_riemann_tol/
//                      exact_riemann_max_iter but not meant to be tuned in
//                      normal use
// Returns: the numerical flux crossing the face (same convention as flux())
inline EulerState exact_riemann_flux(const EulerState& U_L, const EulerState& U_R, double nx, double ny, double gamma,
                                      double rel_tol, int max_iter) {
    double tx = -ny, ty = nx;

    double rho_L = U_L.rho, rho_R = U_R.rho;
    double ux_L = U_L.rho_u / rho_L, uy_L = U_L.rho_v / rho_L;
    double ux_R = U_R.rho_u / rho_R, uy_R = U_R.rho_v / rho_R;
    double u_L = ux_L * nx + uy_L * ny, v_L = ux_L * tx + uy_L * ty;
    double u_R = ux_R * nx + uy_R * ny, v_R = ux_R * tx + uy_R * ty;
    double p_L = pressure(U_L, gamma), p_R = pressure(U_R, gamma);
    double c_L = sound_speed(U_L, gamma), c_R = sound_speed(U_R, gamma);

    double rho_s, u_s, p_s, v_s;

    if (2.0 * (c_L + c_R) / (gamma - 1.0) <= u_R - u_L) {
        double S_L_vac = u_L + 2.0 * c_L / (gamma - 1.0);
        double S_R_vac = u_R - 2.0 * c_R / (gamma - 1.0);
        if (0.0 <= u_L - c_L) {
            rho_s = rho_L; u_s = u_L; p_s = p_L; v_s = v_L;
        } else if (0.0 <= S_L_vac) {
            ExactRiemannDetail::left_fan_interior(rho_L, u_L, p_L, c_L, gamma, 0.0, rho_s, u_s, p_s);
            v_s = v_L;
        } else if (0.0 < S_R_vac) {
            return EulerState{0.0, 0.0, 0.0, 0.0}; // vacuum: rho = p = 0, no mass/momentum/energy flux
        } else if (0.0 <= u_R + c_R) {
            ExactRiemannDetail::right_fan_interior(rho_R, u_R, p_R, c_R, gamma, 0.0, rho_s, u_s, p_s);
            v_s = v_R;
        } else {
            rho_s = rho_R; u_s = u_R; p_s = p_R; v_s = v_R;
        }
    } else {
        double p_star = ExactRiemannDetail::solve_star_pressure(rho_L, u_L, p_L, c_L, rho_R, u_R, p_R, c_R, gamma,
                                                                  rel_tol, max_iter);
        double u_star = 0.5 * (u_L + u_R)
                       + 0.5 * (ExactRiemannDetail::pressure_function(p_star, rho_R, p_R, c_R, gamma)
                              - ExactRiemannDetail::pressure_function(p_star, rho_L, p_L, c_L, gamma));

        if (0.0 <= u_star) {
            ExactRiemannDetail::sample_left(rho_L, u_L, p_L, c_L, gamma, p_star, u_star, 0.0, rho_s, u_s, p_s);
            v_s = v_L;
        } else {
            ExactRiemannDetail::sample_right(rho_R, u_R, p_R, c_R, gamma, p_star, u_star, 0.0, rho_s, u_s, p_s);
            v_s = v_R;
        }
    }

    double ux_s = u_s * nx - v_s * ny;
    double uy_s = u_s * ny + v_s * nx;
    EulerState U_s = from_primitive(rho_s, ux_s, uy_s, p_s, gamma);
    return flux(U_s, nx, ny, gamma);
}

#endif // EXACTRIEMANNFLUX_H_INCLUDED
