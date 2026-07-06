// SPDX-License-Identifier: GPL-3.0-only
#ifndef RUSANOVFLUX_H_INCLUDED
#define RUSANOVFLUX_H_INCLUDED

#include <cmath>
#include <algorithm>

#include "EulerState.h"

// Rusanov (local Lax-Friedrichs) approximate Riemann solver: one of the
// numerical fluxes selectable via EulerFVMSolver's NumericalFluxScheme (see
// EulerFVMSolver.h). Kept separate from EulerState.h, which holds only the
// conserved-state type and its underlying physics (pressure, sound speed,
// the analytic flux()) -- this is a numerical scheme built on top of that
// physics, not the physics itself.
//
// The numerical flux combines the physical fluxes on each side with an
// upwind-biased dissipation term sized to the fastest wave speed present.
// This is what makes the scheme stable for a hyperbolic (convection-
// dominated) system, unlike the pure central differencing that is
// sufficient for the scalar diffusion solver. Cheapest and most robust of
// the schemes here, but the most dissipative -- it applies the same amount
// of dissipation to every wave family, including slow-moving contact/shear
// waves that don't need it.
//
// Methodology: F* = 0.5*(F(U_L) + F(U_R)) - 0.5*S_max*(U_R - U_L), where
// S_max = max(|Vn_L| + c_L, |Vn_R| + c_R) is a conservative bound on the
// fastest signal speed crossing the face in either direction.
//
// Input:   U_L, U_R - conserved states on the left/right of the face;
//          nx, ny   - unit normal pointing from left to right; gamma
// Returns: the numerical flux crossing the face (same convention as flux())
inline EulerState rusanov_flux(const EulerState& U_L, const EulerState& U_R, double nx, double ny, double gamma) {
    double Vn_L = (U_L.rho_u * nx + U_L.rho_v * ny) / U_L.rho;
    double Vn_R = (U_R.rho_u * nx + U_R.rho_v * ny) / U_R.rho;
    double S_max = std::max(std::fabs(Vn_L) + sound_speed(U_L, gamma),
                             std::fabs(Vn_R) + sound_speed(U_R, gamma));

    EulerState F_L = flux(U_L, nx, ny, gamma);
    EulerState F_R = flux(U_R, nx, ny, gamma);

    return 0.5 * (F_L + F_R) - 0.5 * S_max * (U_R - U_L);
}

#endif // RUSANOVFLUX_H_INCLUDED
