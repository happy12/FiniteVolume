// SPDX-License-Identifier: GPL-3.0-only
#include "NavierStokesFVMSolver.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "ExactRiemannFlux.h"
#include "HllcFlux.h"
#include "RusanovFlux.h"

namespace {

// compute_geometry() must run before GradientCalculator's constructor
// inspects the mesh, so it can't be deferred to the constructor body (see
// AdvectionDiffusionFVMSolver.cpp for the same pattern) -- gradient_calc is
// a member initialized in the same initializer list as mesh, before the
// body runs.
UnstructuredMesh mesh_with_geometry(UnstructuredMesh m) {
    m.compute_geometry();
    return m;
}

} // namespace

// See NavierStokesFVMSolver.h for the input/output contract.
NavierStokesFVMSolver::NavierStokesFVMSolver(const UnstructuredMesh& input_mesh,
                                              const std::vector<NSBoundaryCondition>& boundary_conditions,
                                              double gamma, double gas_constant, double mu, double prandtl,
                                              double cfl, const EulerInitialCondition& initial_condition,
                                              NumericalFluxScheme flux_scheme, GradientScheme gradient_scheme,
                                              double exact_riemann_tol, int exact_riemann_max_iter)
    : mesh(mesh_with_geometry(input_mesh)), bcs(boundary_conditions), gamma(gamma), gas_constant(gas_constant),
      mu(mu), prandtl(prandtl), cfl(cfl), flux_scheme(flux_scheme), gradient_calc(mesh, gradient_scheme),
      exact_riemann_tol(exact_riemann_tol), exact_riemann_max_iter(exact_riemann_max_iter)
{
    U.resize(mesh.cells.size());

    for (size_t i = 0; i < mesh.cells.size(); ++i) {
        if (initial_condition.mode == EulerICMode::Freestream) {
            U[i] = from_primitive(initial_condition.rho, initial_condition.u, initial_condition.v,
                                   initial_condition.p, gamma);
        } else if (mesh.cells[i].x_centroid < initial_condition.x0) {
            U[i] = from_primitive(initial_condition.rho_l, initial_condition.u_l, initial_condition.v_l,
                                   initial_condition.p_l, gamma);
        } else {
            U[i] = from_primitive(initial_condition.rho_r, initial_condition.u_r, initial_condition.v_r,
                                   initial_condition.p_r, gamma);
        }
    }
}

// See NavierStokesFVMSolver.h for the input/output contract.
void NavierStokesFVMSolver::build_boundary_fields(const std::vector<double>& u, const std::vector<double>& v,
                                                    const std::vector<double>& T, std::vector<double>& boundary_u,
                                                    std::vector<double>& boundary_v,
                                                    std::vector<double>& boundary_T) const
{
    boundary_u.assign(mesh.faces.size(), 0.0);
    boundary_v.assign(mesh.faces.size(), 0.0);
    boundary_T.assign(mesh.faces.size(), 0.0);

    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        if (face.cell_right != -1) continue;
        int cl = face.cell_left;
        const NSBoundaryCondition& bc = bcs.at(face.patch_id);

        if (bc.type == NSBoundaryType::NoSlipWall) {
            boundary_u[i] = bc.wall_u;
            boundary_v[i] = bc.wall_v;
            boundary_T[i] = bc.is_isothermal_wall ? bc.wall_temperature : T[cl];
        } else {
            // Farfield/Outflow: no viscous flux is used at these faces (see
            // step()), so these values only complete the gradient stencil
            // for neighboring interior cells -- a zero-order extrapolation
            // is good enough.
            boundary_u[i] = u[cl];
            boundary_v[i] = v[cl];
            boundary_T[i] = T[cl];
        }
    }
}

// See NavierStokesFVMSolver.h for the input/output contract.
//
// Methodology: cell gradients of u, v, T are reconstructed once up front
// (needed by every face's viscous term), then flux assembly follows the
// same two-OpenMP-pass pattern as EulerFVMSolver::step() (pass 1: per-face
// flux, no shared writes; pass 2: per-cell residual + integration, no
// shared writes). Each face's total flux is (inviscid - viscous), per the
// standard compressible Navier-Stokes flux-form derivation (see class
// comment); viscous is left at all-zero for Farfield/Outflow boundary faces.
void NavierStokesFVMSolver::step()
{
    double dt = compute_dt();

    std::vector<double> u(mesh.cells.size()), v(mesh.cells.size()), T(mesh.cells.size());
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        u[c] = U[c].rho_u / U[c].rho;
        v[c] = U[c].rho_v / U[c].rho;
        T[c] = temperature(U[c], gamma, gas_constant);
    }

    std::vector<double> boundary_u, boundary_v, boundary_T;
    build_boundary_fields(u, v, T, boundary_u, boundary_v, boundary_T);

    std::vector<Gradient2> grad_u = gradient_calc.compute(mesh, u, boundary_u);
    std::vector<Gradient2> grad_v = gradient_calc.compute(mesh, v, boundary_v);
    std::vector<Gradient2> grad_T = gradient_calc.compute(mesh, T, boundary_T);

    double cp = gamma * gas_constant / (gamma - 1.0);
    double k = mu * cp / prandtl;

    std::vector<EulerState> face_flux(mesh.faces.size());

    // --- Pass 1: per-face flux (parallel over faces) ---
    #pragma omp parallel for
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        int cl = face.cell_left, cr = face.cell_right;
        bool is_boundary = (cr == -1);
        const EulerState& U_L = U[cl];
        EulerState U_R = is_boundary ? ghost_state(face, U_L) : U[cr];

        EulerState inviscid;
        switch (flux_scheme) {
            case NumericalFluxScheme::HLLC:
                inviscid = hllc_flux(U_L, U_R, face.nx, face.ny, gamma);
                break;
            case NumericalFluxScheme::Exact:
                inviscid = exact_riemann_flux(U_L, U_R, face.nx, face.ny, gamma, exact_riemann_tol,
                                                exact_riemann_max_iter);
                break;
            case NumericalFluxScheme::Rusanov:
            default:
                inviscid = rusanov_flux(U_L, U_R, face.nx, face.ny, gamma);
                break;
        }

        EulerState viscous; // all-zero by default (Farfield/Outflow boundary faces)
        bool has_viscous_term = !is_boundary || bcs.at(face.patch_id).type == NSBoundaryType::NoSlipWall;
        if (has_viscous_term) {
            double phi_R_u = is_boundary ? boundary_u[i] : u[cr];
            double phi_R_v = is_boundary ? boundary_v[i] : v[cr];
            double phi_R_T = is_boundary ? boundary_T[i] : T[cr];
            Gradient2 gr_u = is_boundary ? Gradient2{} : grad_u[cr];
            Gradient2 gr_v = is_boundary ? Gradient2{} : grad_v[cr];
            Gradient2 gr_T = is_boundary ? Gradient2{} : grad_T[cr];

            FaceGradient fu = face_gradient(mesh, face, u[cl], phi_R_u, grad_u[cl], gr_u, boundary_u[i]);
            FaceGradient fv = face_gradient(mesh, face, v[cl], phi_R_v, grad_v[cl], gr_v, boundary_v[i]);
            FaceGradient fT = face_gradient(mesh, face, T[cl], phi_R_T, grad_T[cl], gr_T, boundary_T[i]);

            Gradient2 gu = corrected_face_gradient_vector(fu, face.nx, face.ny);
            Gradient2 gv = corrected_face_gradient_vector(fv, face.nx, face.ny);
            Gradient2 gT = corrected_face_gradient_vector(fT, face.nx, face.ny);

            double u_face = is_boundary ? boundary_u[i] : 0.5 * (u[cl] + u[cr]);
            double v_face = is_boundary ? boundary_v[i] : 0.5 * (v[cl] + v[cr]);

            double div_v = gu.dphidx + gv.dphidy;
            double tau_xx = mu * (2.0 * gu.dphidx - (2.0 / 3.0) * div_v);
            double tau_yy = mu * (2.0 * gv.dphidy - (2.0 / 3.0) * div_v);
            double tau_xy = mu * (gu.dphidy + gv.dphidx);

            double qx = -k * gT.dphidx;
            double qy = -k * gT.dphidy;
            if (is_boundary) {
                const NSBoundaryCondition& bc = bcs.at(face.patch_id);
                if (!bc.is_isothermal_wall) {
                    // Adiabatic wall: the heat flux is prescribed directly
                    // (zero), not reconstructed from a gradient -- same
                    // "prescribed flux bypasses gradient math" convention as
                    // a Neumann BC elsewhere in this codebase.
                    qx = 0.0;
                    qy = 0.0;
                }
            }

            viscous.rho_u = tau_xx * face.nx + tau_xy * face.ny;
            viscous.rho_v = tau_xy * face.nx + tau_yy * face.ny;
            viscous.E = (u_face * tau_xx + v_face * tau_xy - qx) * face.nx +
                        (u_face * tau_xy + v_face * tau_yy - qy) * face.ny;
        }

        face_flux[i] = face.area * (inviscid - viscous);
    }

    // --- Pass 2: per-cell residual + integration (parallel over cells) ---
    double sum_sq_rho = 0.0, sum_sq_rho_u = 0.0, sum_sq_rho_v = 0.0, sum_sq_E = 0.0;
    #pragma omp parallel for reduction(+:sum_sq_rho, sum_sq_rho_u, sum_sq_rho_v, sum_sq_E)
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        EulerState residual_c;
        for (int face_idx : mesh.cells[c].faces) {
            const Face& face = mesh.faces[face_idx];
            if ((size_t)face.cell_left == c) residual_c = residual_c - face_flux[face_idx];
            else residual_c += face_flux[face_idx];
        }

        sum_sq_rho += residual_c.rho * residual_c.rho;
        sum_sq_rho_u += residual_c.rho_u * residual_c.rho_u;
        sum_sq_rho_v += residual_c.rho_v * residual_c.rho_v;
        sum_sq_E += residual_c.E * residual_c.E;

        U[c] += (dt / mesh.cells[c].volume) * residual_c;
    }

    last_residual.rho = std::sqrt(sum_sq_rho);
    last_residual.rho_u = std::sqrt(sum_sq_rho_u);
    last_residual.rho_v = std::sqrt(sum_sq_rho_v);
    last_residual.E = std::sqrt(sum_sq_E);
}

// See NavierStokesFVMSolver.h for the input/output contract.
void NavierStokesFVMSolver::run(int total_steps)
{
    for (int t = 0; t < total_steps; ++t) {
        step();
    }
    std::cout << "Simulation completed across " << total_steps << " steps.\n";
}

// See NavierStokesFVMSolver.h for the input/output contract.
double NavierStokesFVMSolver::compute_dt() const
{
    double dt = std::numeric_limits<double>::max();

    #pragma omp parallel for reduction(min:dt)
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        const EulerState& U_L = U[face.cell_left];
        double u_L = U_L.rho_u / U_L.rho;
        double v_L = U_L.rho_v / U_L.rho;
        double Vn_L = u_L * face.nx + v_L * face.ny;
        double Smax = std::fabs(Vn_L) + sound_speed(U_L, gamma);
        double rho_min = U_L.rho;

        if (face.cell_right != -1) {
            const EulerState& U_R = U[face.cell_right];
            double u_R = U_R.rho_u / U_R.rho;
            double v_R = U_R.rho_v / U_R.rho;
            double Vn_R = u_R * face.nx + v_R * face.ny;
            Smax = std::max(Smax, std::fabs(Vn_R) + sound_speed(U_R, gamma));
            rho_min = std::min(rho_min, U_R.rho);
        }
        // Face-normal-projected cell-centroid separation (see
        // GradientReconstruction.h), not volume/face.area: on an
        // anisotropic (e.g. boundary-layer-clustered) cell, volume/area
        // averages over both cell dimensions and understates the true
        // wall-normal stiffness the viscous term below is estimating.
        double length = face_normal_distance(mesh, face);

        double nu = (mu > 0.0) ? mu / rho_min : 0.0; // kinematic viscosity
        double Smax_total = Smax + 2.0 * nu / length;

        dt = std::min(dt, cfl * length / Smax_total);
    }

    return dt;
}

// See NavierStokesFVMSolver.h for the input/output contract.
EulerState NavierStokesFVMSolver::ghost_state(const Face& face, const EulerState& U_L) const
{
    const NSBoundaryCondition& bc = bcs.at(face.patch_id);

    if (bc.type == NSBoundaryType::Farfield) {
        return bc.farfield_state;
    }
    if (bc.type == NSBoundaryType::Outflow) {
        return U_L;
    }

    // NoSlipWall: mirror BOTH velocity components about the wall's own
    // velocity (not just the normal component mirrored about zero, unlike
    // EulerFVMSolver's slip wall), so the inviscid flux's own
    // central-average term comes out to exactly (wall_u, wall_v) at the wall.
    double u_L = U_L.rho_u / U_L.rho;
    double v_L = U_L.rho_v / U_L.rho;
    double p_L = pressure(U_L, gamma);
    return from_primitive(U_L.rho, 2.0 * bc.wall_u - u_L, 2.0 * bc.wall_v - v_L, p_L, gamma);
}

// See NavierStokesFVMSolver.h for the input/output contract and methodology.
ResolutionDiagnostics NavierStokesFVMSolver::compute_resolution_diagnostics() const
{
    std::vector<double> u(mesh.cells.size()), v(mesh.cells.size()), T(mesh.cells.size());
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        u[c] = U[c].rho_u / U[c].rho;
        v[c] = U[c].rho_v / U[c].rho;
        T[c] = temperature(U[c], gamma, gas_constant);
    }

    std::vector<double> boundary_u, boundary_v, boundary_T;
    build_boundary_fields(u, v, T, boundary_u, boundary_v, boundary_T);

    std::vector<Gradient2> grad_u = gradient_calc.compute(mesh, u, boundary_u);
    std::vector<Gradient2> grad_v = gradient_calc.compute(mesh, v, boundary_v);

    const double dissipation_floor = 1e-12; // below this, treat as "no local shear to resolve"

    ResolutionDiagnostics diag;
    diag.min_ratio = std::numeric_limits<double>::max();
    double sum_ratio = 0.0;

    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        double nu = mu / U[c].rho;
        double S11 = grad_u[c].dphidx;
        double S22 = grad_v[c].dphidy;
        double S12 = 0.5 * (grad_u[c].dphidy + grad_v[c].dphidx);
        double dissipation = 2.0 * nu * (S11 * S11 + S22 * S22 + 2.0 * S12 * S12);
        if (dissipation < dissipation_floor) continue;

        double eta = std::pow(nu * nu * nu / dissipation, 0.25);
        double h = std::sqrt(mesh.cells[c].volume);
        double ratio = h / eta;

        diag.min_ratio = std::min(diag.min_ratio, ratio);
        diag.max_ratio = std::max(diag.max_ratio, ratio);
        sum_ratio += ratio;
        ++diag.n_active;
    }

    diag.mean_ratio = (diag.n_active > 0) ? sum_ratio / diag.n_active : 0.0;
    if (diag.n_active == 0) diag.min_ratio = 0.0;

    return diag;
}

// See NavierStokesFVMSolver.h for the input/output contract.
std::vector<WallFaceSample> NavierStokesFVMSolver::compute_wall_traction_samples(const std::vector<int>& wall_faces) const
{
    std::vector<double> u(mesh.cells.size()), v(mesh.cells.size()), T(mesh.cells.size());
    std::vector<double> p(mesh.cells.size()), rho(mesh.cells.size());
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        u[c] = U[c].rho_u / U[c].rho;
        v[c] = U[c].rho_v / U[c].rho;
        T[c] = temperature(U[c], gamma, gas_constant);
        p[c] = pressure(U[c], gamma);
        rho[c] = U[c].rho;
    }

    std::vector<double> boundary_u, boundary_v, boundary_T;
    build_boundary_fields(u, v, T, boundary_u, boundary_v, boundary_T);

    std::vector<Gradient2> grad_u = gradient_calc.compute(mesh, u, boundary_u);
    std::vector<Gradient2> grad_v = gradient_calc.compute(mesh, v, boundary_v);

    std::vector<double> effective_viscosity(mesh.cells.size(), mu);

    return compute_wall_traction(mesh, wall_faces, u, v, p, rho, effective_viscosity, boundary_u, boundary_v,
                                  grad_u, grad_v);
}

// See NavierStokesFVMSolver.h for the input/output contract.
std::vector<BoundaryLayerProfile> NavierStokesFVMSolver::compute_boundary_layer_profile_samples(
    const std::vector<int>& wall_faces, double u_edge, int max_cells_per_march) const
{
    std::vector<double> u(mesh.cells.size()), v(mesh.cells.size()), T(mesh.cells.size());
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        u[c] = U[c].rho_u / U[c].rho;
        v[c] = U[c].rho_v / U[c].rho;
        T[c] = temperature(U[c], gamma, gas_constant);
    }

    std::vector<double> boundary_u, boundary_v, boundary_T;
    build_boundary_fields(u, v, T, boundary_u, boundary_v, boundary_T);

    return compute_boundary_layer_profiles(mesh, wall_faces, u, v, boundary_u, boundary_v, u_edge,
                                            max_cells_per_march);
}

// See NavierStokesFVMSolver.h for the input/output contract.
std::vector<BoundaryLayerProfile> NavierStokesFVMSolver::compute_boundary_layer_profile_samples_point_location(
    const std::vector<int>& wall_faces, double u_edge, double max_distance, int n_samples) const
{
    std::vector<double> u(mesh.cells.size()), v(mesh.cells.size()), T(mesh.cells.size());
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        u[c] = U[c].rho_u / U[c].rho;
        v[c] = U[c].rho_v / U[c].rho;
        T[c] = temperature(U[c], gamma, gas_constant);
    }

    std::vector<double> boundary_u, boundary_v, boundary_T;
    build_boundary_fields(u, v, T, boundary_u, boundary_v, boundary_T);

    return compute_boundary_layer_profiles_point_location(mesh, wall_faces, u, v, boundary_u, boundary_v, u_edge,
                                                            max_distance, n_samples);
}
