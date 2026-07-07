// SPDX-License-Identifier: GPL-3.0-only
#include "RANSTurbulenceSSTSolver.h"

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
std::vector<int> collect_wall_faces(const UnstructuredMesh& mesh, const std::vector<RANSBoundaryConditionSST>& bcs) {
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

// See RANSTurbulenceSSTSolver.h for the input/output contract.
RANSTurbulenceSSTSolver::RANSTurbulenceSSTSolver(const UnstructuredMesh& input_mesh,
                               const std::vector<RANSBoundaryConditionSST>& boundary_conditions, double gamma,
                               double gas_constant, double mu, double prandtl, double prandtl_t, double cfl,
                               const EulerInitialCondition& initial_condition, double initial_k, double initial_omega,
                               NumericalFluxScheme flux_scheme, GradientScheme gradient_scheme,
                               double exact_riemann_tol, int exact_riemann_max_iter, SSTLimiterVariant variant,
                               const SSTModelConstants& sst_constants, bool kato_launder)
    : mesh(mesh_with_geometry(input_mesh)), bcs(boundary_conditions), gamma(gamma), gas_constant(gas_constant),
      mu(mu), prandtl(prandtl), prandtl_t(prandtl_t), cfl(cfl), flux_scheme(flux_scheme),
      gradient_calc(mesh, gradient_scheme), exact_riemann_tol(exact_riemann_tol),
      exact_riemann_max_iter(exact_riemann_max_iter), variant(variant), sst_constants(sst_constants),
      kato_launder(kato_launder)
{
    U.resize(mesh.cells.size());
    k.assign(mesh.cells.size(), initial_k);
    omega.assign(mesh.cells.size(), initial_omega);

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

// See RANSTurbulenceSSTSolver.h for the input/output contract.
void RANSTurbulenceSSTSolver::build_boundary_fields(const std::vector<double>& u, const std::vector<double>& v,
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

// See RANSTurbulenceSSTSolver.h for the input/output contract.
void RANSTurbulenceSSTSolver::build_boundary_k_omega(const std::vector<double>& k_field,
                                                       const std::vector<double>& omega_field,
                                                       std::vector<double>& boundary_k,
                                                       std::vector<double>& boundary_omega,
                                                       std::vector<double>& boundary_rho_k,
                                                       std::vector<double>& boundary_rho_omega) const
{
    boundary_k.assign(mesh.faces.size(), 0.0);
    boundary_omega.assign(mesh.faces.size(), 0.0);
    boundary_rho_k.assign(mesh.faces.size(), 0.0);
    boundary_rho_omega.assign(mesh.faces.size(), 0.0);

    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        if (face.cell_right != -1) continue;
        int cl = face.cell_left;
        const RANSBoundaryConditionSST& bc = bcs.at(face.patch_id);

        if (bc.ns.type == NSBoundaryType::NoSlipWall) {
            boundary_k[i] = 0.0; // the SST model's own wall condition, not case-configurable
            double nu_lam_cl = mu / U[cl].rho;
            double d1 = face_normal_distance(mesh, face);
            boundary_omega[i] = 60.0 * nu_lam_cl / (sst_constants.beta1 * d1 * d1);
            boundary_rho_k[i] = 0.0;
            boundary_rho_omega[i] = U[cl].rho * boundary_omega[i];
        } else if (bc.ns.type == NSBoundaryType::Farfield) {
            boundary_k[i] = bc.farfield_k;
            boundary_omega[i] = bc.farfield_omega;
            double rho_farfield = bc.ns.farfield_state.rho;
            boundary_rho_k[i] = rho_farfield * bc.farfield_k;
            boundary_rho_omega[i] = rho_farfield * bc.farfield_omega;
        } else { // Outflow: zero-order extrapolation
            boundary_k[i] = k_field[cl];
            boundary_omega[i] = omega_field[cl];
            boundary_rho_k[i] = U[cl].rho * k_field[cl];
            boundary_rho_omega[i] = U[cl].rho * omega_field[cl];
        }
    }
}

// See RANSTurbulenceSSTSolver.h for the input/output contract and methodology.
void RANSTurbulenceSSTSolver::step()
{
    double dt = compute_dt();
    const size_t N = mesh.cells.size();

    std::vector<double> u(N), v(N), T(N), nu_lam(N);
    for (size_t c = 0; c < N; ++c) {
        u[c] = U[c].rho_u / U[c].rho;
        v[c] = U[c].rho_v / U[c].rho;
        T[c] = temperature(U[c], gamma, gas_constant);
        nu_lam[c] = mu / U[c].rho;
    }

    std::vector<double> boundary_u, boundary_v, boundary_T;
    build_boundary_fields(u, v, T, boundary_u, boundary_v, boundary_T);
    std::vector<double> boundary_k, boundary_omega, boundary_rho_k, boundary_rho_omega;
    build_boundary_k_omega(k, omega, boundary_k, boundary_omega, boundary_rho_k, boundary_rho_omega);

    std::vector<Gradient2> grad_u = gradient_calc.compute(mesh, u, boundary_u);
    std::vector<Gradient2> grad_v = gradient_calc.compute(mesh, v, boundary_v);
    std::vector<Gradient2> grad_T = gradient_calc.compute(mesh, T, boundary_T);
    std::vector<Gradient2> grad_k = gradient_calc.compute(mesh, k, boundary_k);
    std::vector<Gradient2> grad_omega = gradient_calc.compute(mesh, omega, boundary_omega);

    // Per-cell SST source terms and intermediates (S, Omega, F1, F2, nu_t),
    // from the already-reconstructed gradients plus the wall-distance field
    // precomputed at construction -- see SSTKOmega.h.
    std::vector<SSTSourceTerms> source(N);
    #pragma omp parallel for
    for (size_t c = 0; c < N; ++c) {
        source[c] = compute_sst_source_terms(U[c].rho, k[c], omega[c], nu_lam[c], wall_distance[c], grad_u[c],
                                               grad_v[c], grad_k[c], grad_omega[c], variant, sst_constants,
                                               kato_launder);
    }

    double cp = gamma * gas_constant / (gamma - 1.0);

    std::vector<EulerState> face_flux(mesh.faces.size());
    std::vector<double> k_face_flux(mesh.faces.size()), omega_face_flux(mesh.faces.size());

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

        // k/omega's advective flux upwinds the CONSERVED quantity (rho*k,
        // rho*omega) directly, using the same face-normal velocity Vn as the
        // inviscid mean-flow flux -- see class comment on why k/omega are
        // rho-weighted here unlike RANSTurbulenceSASolver's nut.
        double u_R_ghost = U_R.rho_u / U_R.rho, v_R_ghost = U_R.rho_v / U_R.rho;
        double Vn = 0.5 * ((u[cl] + u_R_ghost) * face.nx + (v[cl] + v_R_ghost) * face.ny);

        double rho_k_L = U[cl].rho * k[cl];
        double rho_k_R_val = is_boundary ? boundary_rho_k[i] : U[cr].rho * k[cr];
        double k_upwind_flux = (Vn >= 0.0) ? rho_k_L : rho_k_R_val;
        double k_advective = Vn * k_upwind_flux * face.area;

        double rho_omega_L = U[cl].rho * omega[cl];
        double rho_omega_R_val = is_boundary ? boundary_rho_omega[i] : U[cr].rho * omega[cr];
        double omega_upwind_flux = (Vn >= 0.0) ? rho_omega_L : rho_omega_R_val;
        double omega_advective = Vn * omega_upwind_flux * face.area;

        EulerState viscous; // all-zero by default (Farfield/Outflow boundary faces)
        double k_diffusive = 0.0, omega_diffusive = 0.0;
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

            // is_boundary here can only mean NoSlipWall (the only boundary
            // type with has_viscous_term true), where boundary_k[i] is
            // exactly 0 -- so nu_t there is exactly 0 too
            // (sst_eddy_viscosity()'s own k<=0 guard), no need to actually
            // evaluate it.
            double nu_t_L = source[cl].nu_t;
            double nu_t_R = is_boundary ? 0.0 : source[cr].nu_t;
            double rho_R = is_boundary ? U[cl].rho : U[cr].rho;
            double turb_dyn_visc_face = 0.5 * (U[cl].rho * nu_t_L + rho_R * nu_t_R);

            double mu_eff = mu + turb_dyn_visc_face;
            double k_eff_meanflow = cp * mu / prandtl + cp * turb_dyn_visc_face / prandtl_t;

            double u_face = is_boundary ? boundary_u[i] : 0.5 * (u[cl] + u[cr]);
            double v_face = is_boundary ? boundary_v[i] : 0.5 * (v[cl] + v[cr]);

            double div_v = gu.dphidx + gv.dphidy;
            double tau_xx = mu_eff * (2.0 * gu.dphidx - (2.0 / 3.0) * div_v);
            double tau_yy = mu_eff * (2.0 * gv.dphidy - (2.0 / 3.0) * div_v);
            double tau_xy = mu_eff * (gu.dphidy + gv.dphidx);

            double qx = -k_eff_meanflow * gT.dphidx;
            double qy = -k_eff_meanflow * gT.dphidy;
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

            // k/omega's own diffusive flux: (mu + sigma*mu_t)*dphidn, sigma_k/
            // sigma_omega blended via this face's F1 (arithmetic-mean of the
            // two adjoining cells', matching turb_dyn_visc_face's own
            // interpolation convention) -- see class comment for why sigma
            // scales ONLY the turbulent part here, unlike SA's nut diffusion
            // which divides the WHOLE (nu_lam+nu_t) sum by a single sigma.
            double F1_R = is_boundary ? source[cl].F1 : source[cr].F1;
            double F1_face = is_boundary ? source[cl].F1 : 0.5 * (source[cl].F1 + F1_R);
            double sigma_k_face = F1_face * sst_constants.sigma_k1 + (1.0 - F1_face) * sst_constants.sigma_k2;
            double sigma_omega_face =
                F1_face * sst_constants.sigma_omega1 + (1.0 - F1_face) * sst_constants.sigma_omega2;
            double mu_eff_k = mu + sigma_k_face * turb_dyn_visc_face;
            double mu_eff_omega = mu + sigma_omega_face * turb_dyn_visc_face;

            double phi_R_k = is_boundary ? boundary_k[i] : k[cr];
            double phi_R_omega = is_boundary ? boundary_omega[i] : omega[cr];
            Gradient2 gr_k = is_boundary ? Gradient2{} : grad_k[cr];
            Gradient2 gr_omega = is_boundary ? Gradient2{} : grad_omega[cr];
            FaceGradient fk = face_gradient(mesh, face, k[cl], phi_R_k, grad_k[cl], gr_k, boundary_k[i]);
            FaceGradient fomega =
                face_gradient(mesh, face, omega[cl], phi_R_omega, grad_omega[cl], gr_omega, boundary_omega[i]);

            k_diffusive = -mu_eff_k * fk.dphidn * face.area;
            omega_diffusive = -mu_eff_omega * fomega.dphidn * face.area;
        }

        face_flux[i] = face.area * (inviscid - viscous);
        k_face_flux[i] = k_advective + k_diffusive;
        omega_face_flux[i] = omega_advective + omega_diffusive;
    }

    // --- Pass 2: per-cell residual + integration (parallel over cells) ---
    double sum_sq_rho = 0.0, sum_sq_rho_u = 0.0, sum_sq_rho_v = 0.0, sum_sq_E = 0.0, sum_sq_k = 0.0,
           sum_sq_omega = 0.0;
    #pragma omp parallel for reduction(+:sum_sq_rho, sum_sq_rho_u, sum_sq_rho_v, sum_sq_E, sum_sq_k, sum_sq_omega)
    for (size_t c = 0; c < N; ++c) {
        EulerState residual_c;
        double k_residual_c = 0.0, omega_residual_c = 0.0;
        for (int face_idx : mesh.cells[c].faces) {
            const Face& face = mesh.faces[face_idx];
            if ((size_t)face.cell_left == c) {
                residual_c = residual_c - face_flux[face_idx];
                k_residual_c -= k_face_flux[face_idx];
                omega_residual_c -= omega_face_flux[face_idx];
            } else {
                residual_c += face_flux[face_idx];
                k_residual_c += k_face_flux[face_idx];
                omega_residual_c += omega_face_flux[face_idx];
            }
        }

        // Volumetric SST source terms (already rho-weighted -- see SSTKOmega.h).
        k_residual_c += mesh.cells[c].volume * (source[c].k_production - source[c].k_destruction);
        omega_residual_c +=
            mesh.cells[c].volume * (source[c].omega_production - source[c].omega_destruction + source[c].cross_diffusion);

        sum_sq_rho += residual_c.rho * residual_c.rho;
        sum_sq_rho_u += residual_c.rho_u * residual_c.rho_u;
        sum_sq_rho_v += residual_c.rho_v * residual_c.rho_v;
        sum_sq_E += residual_c.E * residual_c.E;
        sum_sq_k += k_residual_c * k_residual_c;
        sum_sq_omega += omega_residual_c * omega_residual_c;

        // Conservative rho*k/rho*omega bookkeeping done as local scratch
        // values (see class comment): update the mean flow first (so
        // U[c].rho is the NEW density), then recover primitive k/omega by
        // dividing the updated rho*k/rho*omega by that new density.
        double rho_k_old = U[c].rho * k[c];
        double rho_omega_old = U[c].rho * omega[c];

        U[c] += (dt / mesh.cells[c].volume) * residual_c;

        double rho_k_new = rho_k_old + (dt / mesh.cells[c].volume) * k_residual_c;
        double rho_omega_new = rho_omega_old + (dt / mesh.cells[c].volume) * omega_residual_c;
        k[c] = rho_k_new / U[c].rho;
        omega[c] = rho_omega_new / U[c].rho;
    }

    last_residual.rho = std::sqrt(sum_sq_rho);
    last_residual.rho_u = std::sqrt(sum_sq_rho_u);
    last_residual.rho_v = std::sqrt(sum_sq_rho_v);
    last_residual.E = std::sqrt(sum_sq_E);
    last_k_residual = std::sqrt(sum_sq_k);
    last_omega_residual = std::sqrt(sum_sq_omega);
}

// See RANSTurbulenceSSTSolver.h for the input/output contract.
void RANSTurbulenceSSTSolver::run(int total_steps)
{
    for (int t = 0; t < total_steps; ++t) {
        step();
    }
    std::cout << "Simulation completed across " << total_steps << " steps.\n";
}

// See RANSTurbulenceSSTSolver.h for the input/output contract and methodology.
double RANSTurbulenceSSTSolver::compute_dt() const
{
    double dt = std::numeric_limits<double>::max();
    double beta_max = std::max(sst_constants.beta1, sst_constants.beta2);
    double source_coef = std::max(sst_constants.beta_star, beta_max);

    #pragma omp parallel for reduction(min:dt)
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        int cl = face.cell_left, cr = face.cell_right;
        const EulerState& U_L = U[cl];
        double u_L = U_L.rho_u / U_L.rho;
        double v_L = U_L.rho_v / U_L.rho;
        double Vn_L = u_L * face.nx + v_L * face.ny;
        double Smax = std::fabs(Vn_L) + sound_speed(U_L, gamma);
        double rho_min = U_L.rho;

        // nu_t <= k/omega always (nu_t = a1*k/max(a1*omega, limiter*F2), and
        // the max()'s denominator is always >= a1*omega since limiter*F2 >= 0)
        // -- a cheap, always-conservative upper bound requiring no velocity
        // gradients, unlike the true nu_t (which needs S or Omega). Using
        // this bound here (rather than recomputing fresh gradients just for
        // compute_dt(), or reusing stale ones) can only make dt MORE
        // conservative, never less.
        double nu_t_bound_L = (k[cl] > 0.0 && omega[cl] > 0.0) ? k[cl] / omega[cl] : 0.0;
        double nu_t_bound_max = nu_t_bound_L;
        double omega_max_cell = omega[cl];

        if (cr != -1) {
            const EulerState& U_R = U[cr];
            double u_R = U_R.rho_u / U_R.rho;
            double v_R = U_R.rho_v / U_R.rho;
            double Vn_R = u_R * face.nx + v_R * face.ny;
            Smax = std::max(Smax, std::fabs(Vn_R) + sound_speed(U_R, gamma));
            rho_min = std::min(rho_min, U_R.rho);

            double nu_t_bound_R = (k[cr] > 0.0 && omega[cr] > 0.0) ? k[cr] / omega[cr] : 0.0;
            nu_t_bound_max = std::max(nu_t_bound_max, nu_t_bound_R);
            omega_max_cell = std::max(omega_max_cell, omega[cr]);
        }

        double length = face_normal_distance(mesh, face);

        double nu_lam = (mu > 0.0) ? mu / rho_min : 0.0;
        double nu_eff = nu_lam + nu_t_bound_max; // sigma_k/sigma_omega <= 1.0 always -- see class comment
        double Smax_total = Smax + 2.0 * nu_eff / length;
        double dt_diffusion = cfl * length / Smax_total;

        // Explicit-source stability limit for k/omega's own destruction
        // terms (per-unit-mass rates beta_star*omega, beta*omega -- see
        // class comment); beta_max is a safe upper bound on the true
        // F1-blended beta without needing fresh gradients here.
        double dt_source = cfl / (source_coef * std::max(omega_max_cell, 1e-30));

        dt = std::min(dt, std::min(dt_diffusion, dt_source));
    }

    return dt;
}

// See RANSTurbulenceSSTSolver.h for the input/output contract.
std::vector<WallFaceSample> RANSTurbulenceSSTSolver::compute_wall_traction_samples(const std::vector<int>& wall_faces) const
{
    const size_t N = mesh.cells.size();
    std::vector<double> u(N), v(N), T(N), p(N), rho(N), effective_viscosity(N);
    for (size_t c = 0; c < N; ++c) {
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

    for (size_t c = 0; c < N; ++c) {
        double S = sst_strain_rate_magnitude(grad_u[c], grad_v[c]);
        double Omega = sst_vorticity_magnitude(grad_u[c], grad_v[c]);
        double nu_lam_c = mu / U[c].rho;
        double F2c = sst_F2(k[c], omega[c], nu_lam_c, wall_distance[c], sst_constants);
        double nu_t_c = sst_eddy_viscosity(k[c], omega[c], S, Omega, F2c, variant, sst_constants);
        effective_viscosity[c] = mu + U[c].rho * nu_t_c;
    }

    return compute_wall_traction(mesh, wall_faces, u, v, p, rho, effective_viscosity, boundary_u, boundary_v,
                                  grad_u, grad_v);
}

// See RANSTurbulenceSSTSolver.h for the input/output contract.
std::vector<BoundaryLayerProfile> RANSTurbulenceSSTSolver::compute_boundary_layer_profile_samples(
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

// See RANSTurbulenceSSTSolver.h for the input/output contract.
std::vector<BoundaryLayerProfile> RANSTurbulenceSSTSolver::compute_boundary_layer_profile_samples_point_location(
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

// See RANSTurbulenceSSTSolver.h for the input/output contract.
EulerState RANSTurbulenceSSTSolver::ghost_state(const Face& face, const EulerState& U_L) const
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
