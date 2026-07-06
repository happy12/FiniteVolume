// SPDX-License-Identifier: GPL-3.0-only
#ifndef HLLCFLUX_H_INCLUDED
#define HLLCFLUX_H_INCLUDED

#include <cmath>
#include <algorithm>

#include "EulerState.h"

// HLLC (Harten-Lax-van Leer-Contact) approximate Riemann solver: one of the
// numerical fluxes selectable via EulerFVMSolver's NumericalFluxScheme (see
// EulerFVMSolver.h). Restores the middle (contact/shear) wave that plain
// HLL/Rusanov collapse into a single jump, giving much sharper resolution
// of contact discontinuities and shear layers while remaining robust and
// positivity-preserving. See Toro, "Riemann Solvers and Numerical Methods
// for Fluid Dynamics" 3rd ed., Ch. 10, and Toro/Spruce/Speares (1994),
// "Restoration of the contact surface in the HLL-Riemann solver," Shock
// Waves 4(1).
//
// Methodology: the state is rotated into the face-normal/tangential frame
// (Vn/Vt) since the 1D HLLC derivation is along a single direction. Wave
// speeds SL/SR use Davis' simple estimate (Davis 1988): SL = min(Vn_L - c_L,
// Vn_R - c_R), SR = max(Vn_L + c_L, Vn_R + c_R). If 0 lies outside [SL, SR],
// the flux is just the upwind physical flux (supersonic case). Otherwise the
// contact speed S* (Toro eq. 10.37) locates which star region (left or
// right of the contact) the face sits in, and the star-side flux is built
// directly from the conservative jump condition F*_K = F_K + S_K*(U*_K -
// U_K) (Toro eq. 10.38), with U*_K from Toro eq. 10.39: tangential momentum
// is carried through unchanged from side K, since a contact wave doesn't
// affect tangential velocity.
//
// Input:   U_L, U_R - conserved states on the left/right of the face;
//          nx, ny   - unit normal pointing from left to right; gamma
// Returns: the numerical flux crossing the face (same convention as flux())
inline EulerState hllc_flux(const EulerState& U_L, const EulerState& U_R, double nx, double ny, double gamma) {
    // Tangent direction: 90 degrees counter-clockwise from (nx, ny).
    double tx = -ny, ty = nx;

    double rho_L = U_L.rho, rho_R = U_R.rho;
    double u_L = U_L.rho_u / rho_L, v_L = U_L.rho_v / rho_L;
    double u_R = U_R.rho_u / rho_R, v_R = U_R.rho_v / rho_R;
    double Vn_L = u_L * nx + v_L * ny, Vn_R = u_R * nx + v_R * ny;
    double Vt_L = u_L * tx + v_L * ty, Vt_R = u_R * tx + v_R * ty;
    double p_L = pressure(U_L, gamma), p_R = pressure(U_R, gamma);
    double c_L = sound_speed(U_L, gamma), c_R = sound_speed(U_R, gamma);

    double SL = std::min(Vn_L - c_L, Vn_R - c_R);
    double SR = std::max(Vn_L + c_L, Vn_R + c_R);

    if (SL >= 0.0) return flux(U_L, nx, ny, gamma);
    if (SR <= 0.0) return flux(U_R, nx, ny, gamma);

    double S_star = (p_R - p_L + rho_L * Vn_L * (SL - Vn_L) - rho_R * Vn_R * (SR - Vn_R))
                   / (rho_L * (SL - Vn_L) - rho_R * (SR - Vn_R));

    bool left = (S_star >= 0.0);
    double S_K = left ? SL : SR;
    double rho_K = left ? rho_L : rho_R;
    double Vn_K = left ? Vn_L : Vn_R;
    double Vt_K = left ? Vt_L : Vt_R;
    double p_K = left ? p_L : p_R;
    const EulerState& U_K = left ? U_L : U_R;

    double factor = rho_K * (S_K - Vn_K) / (S_K - S_star);
    double mom_n_star = factor * S_star;
    double mom_t_star = factor * Vt_K;
    double E_star = factor * (U_K.E / rho_K + (S_star - Vn_K) * (S_star + p_K / (rho_K * (S_K - Vn_K))));

    EulerState U_star;
    U_star.rho = factor;
    U_star.rho_u = mom_n_star * nx - mom_t_star * ny;
    U_star.rho_v = mom_n_star * ny + mom_t_star * nx;
    U_star.E = E_star;

    return flux(U_K, nx, ny, gamma) + S_K * (U_star - U_K);
}

#endif // HLLCFLUX_H_INCLUDED
