// SPDX-License-Identifier: GPL-3.0-only
#include "RANSTurbulenceSASolver.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "ExactRiemannFlux.h"
#include "HllcFlux.h"
#include "RusanovFlux.h"
#include "WallDistance.h"

namespace {

// See NavierStokesFVMSolver.cpp for why this can't be deferred to the
// constructor body -- gradient_calc is initialized in the same initializer
// list as mesh, before it, and needs compute_geometry() to have already run.
UnstructuredMesh mesh_with_geometry(UnstructuredMesh m) {
    m.compute_geometry();
    return m;
}

// Collects the indices of every boundary face whose patch's NSBoundaryType
// is NoSlipWall, for compute_wall_distance() (WallDistance.h).
std::vector<int> collect_wall_faces(const UnstructuredMesh& mesh, const std::vector<RANSBoundaryConditionSA>& bcs) {
    std::vector<int> wall_faces;
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        if (face.cell_right == -1 && bcs.at(face.patch_id).ns.type == NSBoundaryType::NoSlipWall) {
            wall_faces.push_back((int)i);
        }
    }
    return wall_faces;
}

} // namespace

// See RANSTurbulenceSASolver.h for the input/output contract.
RANSTurbulenceSASolver::RANSTurbulenceSASolver(const UnstructuredMesh& input_mesh, const std::vector<RANSBoundaryConditionSA>& boundary_conditions,
                               double gamma, double gas_constant, double mu, double prandtl, double prandtl_t,
                               double cfl, const EulerInitialCondition& initial_condition, double initial_nut,
                               NumericalFluxScheme flux_scheme, GradientScheme gradient_scheme,
                               double exact_riemann_tol, int exact_riemann_max_iter,
                               const SAModelConstants& sa_constants)
    : mesh(mesh_with_geometry(input_mesh)), bcs(boundary_conditions), gamma(gamma), gas_constant(gas_constant),
      mu(mu), prandtl(prandtl), prandtl_t(prandtl_t), cfl(cfl), flux_scheme(flux_scheme),
      gradient_calc(mesh, gradient_scheme), exact_riemann_tol(exact_riemann_tol),
      exact_riemann_max_iter(exact_riemann_max_iter), sa_constants(sa_constants)
{
    U.resize(mesh.cells.size());
    nut.assign(mesh.cells.size(), initial_nut);

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

    wall_distance = compute_wall_distance(mesh, collect_wall_faces(mesh, bcs));
}

// See RANSTurbulenceSASolver.h for the input/output contract.
void RANSTurbulenceSASolver::build_boundary_fields(const std::vector<double>& u, const std::vector<double>& v,
                                            const std::vector<double>& T, std::vector<double>& boundary_u,
                                            std::vector<double>& boundary_v, std::vector<double>& boundary_T) const
{
    boundary_u.assign(mesh.faces.size(), 0.0);
    boundary_v.assign(mesh.faces.size(), 0.0);
    boundary_T.assign(mesh.faces.size(), 0.0);

    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        if (face.cell_right != -1) continue;
        int cl = face.cell_left;
        const NSBoundaryCondition& bc = bcs.at(face.patch_id).ns;

        if (bc.type == NSBoundaryType::NoSlipWall) {
            boundary_u[i] = bc.wall_u;
            boundary_v[i] = bc.wall_v;
            boundary_T[i] = bc.is_isothermal_wall ? bc.wall_temperature : T[cl];
        } else {
            boundary_u[i] = u[cl];
            boundary_v[i] = v[cl];
            boundary_T[i] = T[cl];
        }
    }
}

// See RANSTurbulenceSASolver.h for the input/output contract.
void RANSTurbulenceSASolver::build_boundary_nut(const std::vector<double>& nut_field, std::vector<double>& boundary_nut) const
{
    boundary_nut.assign(mesh.faces.size(), 0.0);

    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        if (face.cell_right != -1) continue;
        const RANSBoundaryConditionSA& bc = bcs.at(face.patch_id);

        if (bc.ns.type == NSBoundaryType::NoSlipWall) {
            boundary_nut[i] = 0.0; // the SA model's own wall condition, not case-configurable
        } else if (bc.ns.type == NSBoundaryType::Farfield) {
            boundary_nut[i] = bc.farfield_nut;
        } else {
            boundary_nut[i] = nut_field[face.cell_left]; // Outflow: zero-order extrapolation
        }
    }
}

// See RANSTurbulenceSASolver.h for the input/output contract and methodology.
void RANSTurbulenceSASolver::step()
{
    double dt = compute_dt();
    const size_t N = mesh.cells.size();

    std::vector<double> u(N), v(N), T(N), nu_lam(N), nu_t(N);
    for (size_t c = 0; c < N; ++c) {
        u[c] = U[c].rho_u / U[c].rho;
        v[c] = U[c].rho_v / U[c].rho;
        T[c] = temperature(U[c], gamma, gas_constant);
        nu_lam[c] = mu / U[c].rho;
        nu_t[c] = sa_eddy_viscosity(nut[c], nu_lam[c], sa_constants);
    }

    std::vector<double> boundary_u, boundary_v, boundary_T, boundary_nut;
    build_boundary_fields(u, v, T, boundary_u, boundary_v, boundary_T);
    build_boundary_nut(nut, boundary_nut);

    std::vector<Gradient2> grad_u = gradient_calc.compute(mesh, u, boundary_u);
    std::vector<Gradient2> grad_v = gradient_calc.compute(mesh, v, boundary_v);
    std::vector<Gradient2> grad_T = gradient_calc.compute(mesh, T, boundary_T);
    std::vector<Gradient2> grad_nut = gradient_calc.compute(mesh, nut, boundary_nut);

    // Per-cell SA source terms (volumetric, not face fluxes -- see class
    // comment), from the already-reconstructed grad_u/grad_v (vorticity) and
    // grad_nut, plus the wall-distance field precomputed at construction.
    std::vector<SASourceTerms> source(N);
    #pragma omp parallel for
    for (size_t c = 0; c < N; ++c) {
        double omega = std::fabs(grad_v[c].dphidx - grad_u[c].dphidy);
        source[c] = compute_sa_source_terms(nut[c], nu_lam[c], omega, wall_distance[c], grad_nut[c], sa_constants);
    }

    double cp = gamma * gas_constant / (gamma - 1.0);

    std::vector<EulerState> face_flux(mesh.faces.size());
    std::vector<double> nut_face_flux(mesh.faces.size());

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

        // nut's advective flux uses the same ghost/interior state pairing as
        // the inviscid mean-flow flux, so a Farfield boundary's real
        // prescribed velocity (not an extrapolated one) sets the convecting
        // velocity there -- see class comment.
        double u_R_ghost = U_R.rho_u / U_R.rho, v_R_ghost = U_R.rho_v / U_R.rho;
        double Vn = 0.5 * ((u[cl] + u_R_ghost) * face.nx + (v[cl] + v_R_ghost) * face.ny);
        double nut_R_value = is_boundary ? boundary_nut[i] : nut[cr];
        double nut_upwind = (Vn >= 0.0) ? nut[cl] : nut_R_value;
        double nut_advective = Vn * nut_upwind * face.area;

        EulerState viscous; // all-zero by default (Farfield/Outflow boundary faces)
        double nut_diffusive = 0.0;
        bool has_viscous_term = !is_boundary || bcs.at(face.patch_id).ns.type == NSBoundaryType::NoSlipWall;
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

            // Face-interpolated (arithmetic-mean) turbulent dynamic
            // viscosity, feeding both mu_eff and k_eff below. At a
            // NoSlipWall boundary face, the "R" side's nu_t is evaluated
            // from boundary_nut[i] (exactly 0 there), so this correctly
            // blends toward zero turbulent viscosity right at the wall
            // rather than retaining the interior cell's full value.
            double nu_lam_R = is_boundary ? nu_lam[cl] : nu_lam[cr];
            double nu_t_R = is_boundary ? sa_eddy_viscosity(boundary_nut[i], nu_lam_R, sa_constants) : nu_t[cr];
            double rho_R = is_boundary ? U[cl].rho : U[cr].rho;
            double turb_dyn_visc_face = 0.5 * (U[cl].rho * nu_t[cl] + rho_R * nu_t_R);

            double mu_eff = mu + turb_dyn_visc_face;
            double k_eff = cp * mu / prandtl + cp * turb_dyn_visc_face / prandtl_t;

            double u_face = is_boundary ? boundary_u[i] : 0.5 * (u[cl] + u[cr]);
            double v_face = is_boundary ? boundary_v[i] : 0.5 * (v[cl] + v[cr]);

            double div_v = gu.dphidx + gv.dphidy;
            double tau_xx = mu_eff * (2.0 * gu.dphidx - (2.0 / 3.0) * div_v);
            double tau_yy = mu_eff * (2.0 * gv.dphidy - (2.0 / 3.0) * div_v);
            double tau_xy = mu_eff * (gu.dphidy + gv.dphidx);

            double qx = -k_eff * gT.dphidx;
            double qy = -k_eff * gT.dphidy;
            if (is_boundary) {
                const NSBoundaryCondition& bc = bcs.at(face.patch_id).ns;
                if (!bc.is_isothermal_wall) {
                    qx = 0.0;
                    qy = 0.0;
                }
            }

            viscous.rho_u = tau_xx * face.nx + tau_xy * face.ny;
            viscous.rho_v = tau_xy * face.nx + tau_yy * face.ny;
            viscous.E = (u_face * tau_xx + v_face * tau_xy - qx) * face.nx +
                        (u_face * tau_xy + v_face * tau_yy - qy) * face.ny;

            // nut's own diffusive flux: (nu_lam + nu_t)/sigma * dphidn,
            // forced to zero at Farfield/Outflow exactly like the mean-flow
            // viscous flux above, per the class comment's stated convention.
            double nu_lam_face = 0.5 * (nu_lam[cl] + nu_lam_R);
            double nu_t_face = 0.5 * (nu_t[cl] + nu_t_R);
            Gradient2 gr_nut = is_boundary ? Gradient2{} : grad_nut[cr];
            FaceGradient fnut = face_gradient(mesh, face, nut[cl], nut_R_value, grad_nut[cl], gr_nut, boundary_nut[i]);
            double diffusivity = (nu_lam_face + nu_t_face) / sa_constants.sigma;
            nut_diffusive = -diffusivity * fnut.dphidn * face.area;
        }

        face_flux[i] = face.area * (inviscid - viscous);
        nut_face_flux[i] = nut_advective + nut_diffusive;
    }

    // --- Pass 2: per-cell residual + integration (parallel over cells) ---
    double sum_sq_rho = 0.0, sum_sq_rho_u = 0.0, sum_sq_rho_v = 0.0, sum_sq_E = 0.0, sum_sq_nut = 0.0;
    #pragma omp parallel for reduction(+:sum_sq_rho, sum_sq_rho_u, sum_sq_rho_v, sum_sq_E, sum_sq_nut)
    for (size_t c = 0; c < N; ++c) {
        EulerState residual_c;
        double nut_residual_c = 0.0;
        for (int face_idx : mesh.cells[c].faces) {
            const Face& face = mesh.faces[face_idx];
            if ((size_t)face.cell_left == c) {
                residual_c = residual_c - face_flux[face_idx];
                nut_residual_c -= nut_face_flux[face_idx];
            } else {
                residual_c += face_flux[face_idx];
                nut_residual_c += nut_face_flux[face_idx];
            }
        }

        // Volumetric SA source term: production - destruction + cross-diffusion.
        nut_residual_c += mesh.cells[c].volume * (source[c].production - source[c].destruction + source[c].cross_diffusion);

        sum_sq_rho += residual_c.rho * residual_c.rho;
        sum_sq_rho_u += residual_c.rho_u * residual_c.rho_u;
        sum_sq_rho_v += residual_c.rho_v * residual_c.rho_v;
        sum_sq_E += residual_c.E * residual_c.E;
        sum_sq_nut += nut_residual_c * nut_residual_c;

        U[c] += (dt / mesh.cells[c].volume) * residual_c;
        nut[c] += (dt / mesh.cells[c].volume) * nut_residual_c;
    }

    last_residual.rho = std::sqrt(sum_sq_rho);
    last_residual.rho_u = std::sqrt(sum_sq_rho_u);
    last_residual.rho_v = std::sqrt(sum_sq_rho_v);
    last_residual.E = std::sqrt(sum_sq_E);
    last_nut_residual = std::sqrt(sum_sq_nut);
}

// See RANSTurbulenceSASolver.h for the input/output contract.
void RANSTurbulenceSASolver::run(int total_steps)
{
    for (int t = 0; t < total_steps; ++t) {
        step();
    }
    std::cout << "Simulation completed across " << total_steps << " steps.\n";
}

// See RANSTurbulenceSASolver.h for the input/output contract.
double RANSTurbulenceSASolver::compute_dt() const
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
        double nu_t_L = sa_eddy_viscosity(nut[face.cell_left], mu / U_L.rho, sa_constants);
        double nu_t_max = nu_t_L;

        if (face.cell_right != -1) {
            const EulerState& U_R = U[face.cell_right];
            double u_R = U_R.rho_u / U_R.rho;
            double v_R = U_R.rho_v / U_R.rho;
            double Vn_R = u_R * face.nx + v_R * face.ny;
            Smax = std::max(Smax, std::fabs(Vn_R) + sound_speed(U_R, gamma));
            rho_min = std::min(rho_min, U_R.rho);

            double nu_t_R = sa_eddy_viscosity(nut[face.cell_right], mu / U_R.rho, sa_constants);
            nu_t_max = std::max(nu_t_max, nu_t_R);
        }
        // Face-normal-projected cell-centroid separation (see
        // GradientReconstruction.h), not volume/face.area: on an
        // anisotropic (e.g. boundary-layer-clustered) cell, volume/area
        // averages over both cell dimensions and understates the true
        // wall-normal stiffness the viscous term below is estimating.
        double length = face_normal_distance(mesh, face);

        double nu_lam = (mu > 0.0) ? mu / rho_min : 0.0; // molecular kinematic viscosity
        double nu_eff = nu_lam + nu_t_max;                // + turbulent kinematic viscosity -- see class comment
        double Smax_total = Smax + 2.0 * nu_eff / length;

        dt = std::min(dt, cfl * length / Smax_total);
    }

    return dt;
}

// See RANSTurbulenceSASolver.h for the input/output contract.
std::vector<WallFaceSample> RANSTurbulenceSASolver::compute_wall_traction_samples(const std::vector<int>& wall_faces) const
{
    const size_t N = mesh.cells.size();
    std::vector<double> u(N), v(N), T(N), p(N), rho(N), effective_viscosity(N);
    for (size_t c = 0; c < N; ++c) {
        u[c] = U[c].rho_u / U[c].rho;
        v[c] = U[c].rho_v / U[c].rho;
        T[c] = temperature(U[c], gamma, gas_constant);
        p[c] = pressure(U[c], gamma);
        rho[c] = U[c].rho;
        double nu_t_c = sa_eddy_viscosity(nut[c], mu / U[c].rho, sa_constants);
        effective_viscosity[c] = mu + U[c].rho * nu_t_c;
    }

    std::vector<double> boundary_u, boundary_v, boundary_T;
    build_boundary_fields(u, v, T, boundary_u, boundary_v, boundary_T);

    std::vector<Gradient2> grad_u = gradient_calc.compute(mesh, u, boundary_u);
    std::vector<Gradient2> grad_v = gradient_calc.compute(mesh, v, boundary_v);

    return compute_wall_traction(mesh, wall_faces, u, v, p, rho, effective_viscosity, boundary_u, boundary_v,
                                  grad_u, grad_v);
}

// See RANSTurbulenceSASolver.h for the input/output contract.
std::vector<BoundaryLayerProfile> RANSTurbulenceSASolver::compute_boundary_layer_profile_samples(
    const std::vector<int>& wall_faces, double u_edge, int max_cells_per_march) const
{
    const size_t N = mesh.cells.size();
    std::vector<double> u(N), v(N), T(N);
    for (size_t c = 0; c < N; ++c) {
        u[c] = U[c].rho_u / U[c].rho;
        v[c] = U[c].rho_v / U[c].rho;
        T[c] = temperature(U[c], gamma, gas_constant);
    }

    std::vector<double> boundary_u, boundary_v, boundary_T;
    build_boundary_fields(u, v, T, boundary_u, boundary_v, boundary_T);

    return compute_boundary_layer_profiles(mesh, wall_faces, u, v, boundary_u, boundary_v, u_edge,
                                            max_cells_per_march);
}

// See RANSTurbulenceSASolver.h for the input/output contract.
std::vector<BoundaryLayerProfile> RANSTurbulenceSASolver::compute_boundary_layer_profile_samples_point_location(
    const std::vector<int>& wall_faces, double u_edge, double max_distance, int n_samples) const
{
    const size_t N = mesh.cells.size();
    std::vector<double> u(N), v(N), T(N);
    for (size_t c = 0; c < N; ++c) {
        u[c] = U[c].rho_u / U[c].rho;
        v[c] = U[c].rho_v / U[c].rho;
        T[c] = temperature(U[c], gamma, gas_constant);
    }

    std::vector<double> boundary_u, boundary_v, boundary_T;
    build_boundary_fields(u, v, T, boundary_u, boundary_v, boundary_T);

    return compute_boundary_layer_profiles_point_location(mesh, wall_faces, u, v, boundary_u, boundary_v, u_edge,
                                                            max_distance, n_samples);
}

// See RANSTurbulenceSASolver.h for the input/output contract.
EulerState RANSTurbulenceSASolver::ghost_state(const Face& face, const EulerState& U_L) const
{
    const NSBoundaryCondition& bc = bcs.at(face.patch_id).ns;

    if (bc.type == NSBoundaryType::Farfield) {
        return bc.farfield_state;
    }
    if (bc.type == NSBoundaryType::Outflow) {
        return U_L;
    }

    double u_L = U_L.rho_u / U_L.rho;
    double v_L = U_L.rho_v / U_L.rho;
    double p_L = pressure(U_L, gamma);
    return from_primitive(U_L.rho, 2.0 * bc.wall_u - u_L, 2.0 * bc.wall_v - v_L, p_L, gamma);
}
