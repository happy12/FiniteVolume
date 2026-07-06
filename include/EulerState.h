// SPDX-License-Identifier: GPL-3.0-only
#ifndef EULERSTATE_H_INCLUDED
#define EULERSTATE_H_INCLUDED

#include <cmath>

// Conserved-variable state for the 2D compressible Euler equations, one
// instance per mesh cell. Units are left generic (consistent length/mass/time
// units chosen by the case setup), matching UnstructuredMesh's convention.
//
//   rho   - density
//   rho_u - x-momentum per volume (rho * u)
//   rho_v - y-momentum per volume (rho * v)
//   E     - total energy per volume (internal + kinetic)
struct EulerState {
    double rho = 0.0;
    double rho_u = 0.0;
    double rho_v = 0.0;
    double E = 0.0;
};

// Elementwise arithmetic, needed to accumulate/scale fluxes and residuals the
// same way plain doubles were combined in the scalar diffusion solver.
inline EulerState operator+(const EulerState& a, const EulerState& b) {
    return {a.rho + b.rho, a.rho_u + b.rho_u, a.rho_v + b.rho_v, a.E + b.E};
}
inline EulerState operator-(const EulerState& a, const EulerState& b) {
    return {a.rho - b.rho, a.rho_u - b.rho_u, a.rho_v - b.rho_v, a.E - b.E};
}
inline EulerState operator*(double s, const EulerState& a) {
    return {s * a.rho, s * a.rho_u, s * a.rho_v, s * a.E};
}
inline EulerState& operator+=(EulerState& a, const EulerState& b) {
    a.rho += b.rho; a.rho_u += b.rho_u; a.rho_v += b.rho_v; a.E += b.E;
    return a;
}

// Builds a conserved EulerState from primitive flow variables -- the natural
// way to specify initial/farfield conditions.
// Input:   rho, u, v - density and velocity components; p - pressure;
//          gamma     - ratio of specific heats (1.4 for air)
// Returns: the equivalent EulerState, using
//          E = p/(gamma-1) + 0.5*rho*(u^2+v^2)
//          (ideal-gas internal energy plus kinetic energy, per unit volume)
inline EulerState from_primitive(double rho, double u, double v, double p, double gamma) {
    EulerState s;
    s.rho = rho;
    s.rho_u = rho * u;
    s.rho_v = rho * v;
    s.E = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);
    return s;
}

// Recovers static pressure from a conserved state via the ideal-gas equation
// of state: p = (gamma - 1) * (E - kinetic energy per volume).
// Input:   s - conserved state; gamma - ratio of specific heats
// Returns: pressure
inline double pressure(const EulerState& s, double gamma) {
    double u = s.rho_u / s.rho;
    double v = s.rho_v / s.rho;
    return (gamma - 1.0) * (s.E - 0.5 * s.rho * (u * u + v * v));
}

// Recovers static temperature from a conserved state via the ideal-gas
// equation of state p = rho*R*T. R (the specific gas constant) is not part
// of EulerState itself -- the Euler equations never need it (pressure and
// sound speed are R-independent) -- only a caller doing viscous heat
// conduction (Fourier's law needs an actual temperature, not just pressure)
// does, so it's a parameter here rather than a stored field.
// Input:   s - conserved state; gamma - ratio of specific heats;
//          gas_constant - specific gas constant, R
// Returns: temperature
inline double temperature(const EulerState& s, double gamma, double gas_constant) {
    return pressure(s, gamma) / (s.rho * gas_constant);
}

// Local speed of sound: c = sqrt(gamma * p / rho).
// Input:   s - conserved state; gamma - ratio of specific heats
// Returns: sound speed (same length/time units as velocity)
inline double sound_speed(const EulerState& s, double gamma) {
    return std::sqrt(gamma * pressure(s, gamma) / s.rho);
}

// Physical (analytic) Euler flux in the direction of unit normal (nx, ny):
//   F(U)*n = [rho*Vn, rho*u*Vn + p*nx, rho*v*Vn + p*ny, (E+p)*Vn]
// where Vn = u*nx + v*ny is the velocity component normal to the face.
// Input:   s - conserved state; nx, ny - unit normal components; gamma
// Returns: the flux vector crossing a unit-length face with that normal
inline EulerState flux(const EulerState& s, double nx, double ny, double gamma) {
    double u = s.rho_u / s.rho;
    double v = s.rho_v / s.rho;
    double Vn = u * nx + v * ny;
    double p = pressure(s, gamma);

    EulerState f;
    f.rho = s.rho * Vn;
    f.rho_u = s.rho_u * Vn + p * nx;
    f.rho_v = s.rho_v * Vn + p * ny;
    f.E = (s.E + p) * Vn;
    return f;
}

#endif // EULERSTATE_H_INCLUDED
