// SPDX-License-Identifier: GPL-3.0-only
#include "EulerFVMSolver.h"

#include <algorithm>
#include <limits>

#include "RusanovFlux.h"
#include "HllcFlux.h"
#include "ExactRiemannFlux.h"

// See EulerFVMSolver.h for the input/output contract.
EulerFVMSolver::EulerFVMSolver(const UnstructuredMesh& input_mesh,
                                const std::vector<EulerBoundaryCondition>& boundary_conditions,
                                double gamma, double cfl,
                                const EulerInitialCondition& initial_condition,
                                NumericalFluxScheme flux_scheme,
                                double exact_riemann_tol, int exact_riemann_max_iter)
    : mesh(input_mesh), bcs(boundary_conditions), gamma(gamma), cfl(cfl), flux_scheme(flux_scheme),
      exact_riemann_tol(exact_riemann_tol), exact_riemann_max_iter(exact_riemann_max_iter)
{
    mesh.compute_geometry();
    U.resize(mesh.cells.size());

    // Initial Condition: either a uniform free-stream state, or a left/right
    // state split at x = x0 (e.g. a Sod shock tube), by cell centroid.
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

// See EulerFVMSolver.h for the input/output contract.
//
// Methodology: a single explicit forward-Euler finite-volume step, using
// this solver's selected numerical flux (see NumericalFluxScheme/
// RiemannSolvers.h) at every face, split into two OpenMP-parallel passes to
// avoid a scatter-add race (two faces writing to the same cell's residual
// concurrently):
//   1. Recompute a stable dt from the CFL condition (see compute_dt()).
//   2. For every face (parallel, no shared writes -- each face only writes
//      its own face_flux[i] entry), compute the numerical flux crossing it
//      (against the real neighbor cell for internal faces, or a
//      boundary-condition-derived ghost state for boundary faces).
//   3. For every cell (parallel, no shared writes -- each cell only writes
//      its own U[c]), sum the fluxes of its own bounding faces (via
//      Cell::faces), with a sign of -1 if the cell is that face's cell_left
//      and +1 if it's cell_right -- equivalent to face-major scatter-add,
//      just re-ordered to be race-free -- then integrate:
//      U_new = U_old + dt/volume * residual (explicit Euler).
void EulerFVMSolver::step()
{
    double dt = compute_dt();
    std::vector<EulerState> face_flux(mesh.faces.size());

    // --- Pass 1: per-face flux (parallel over faces) ---
    #pragma omp parallel for
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        const EulerState& U_L = U[face.cell_left];
        EulerState U_R = (face.cell_right != -1) ? U[face.cell_right] : ghost_state(face, U_L);

        switch (flux_scheme) {
            case NumericalFluxScheme::HLLC:
                face_flux[i] = face.area * hllc_flux(U_L, U_R, face.nx, face.ny, gamma);
                break;
            case NumericalFluxScheme::Exact:
                face_flux[i] = face.area * exact_riemann_flux(U_L, U_R, face.nx, face.ny, gamma,
                                                                exact_riemann_tol, exact_riemann_max_iter);
                break;
            case NumericalFluxScheme::Rusanov:
                face_flux[i] = face.area * rusanov_flux(U_L, U_R, face.nx, face.ny, gamma);
                break;
        }
    }

    // --- Pass 2: per-cell residual + integration (parallel over cells) ---
    // Reduction variables must be plain scalars, so the four conserved
    // variables' sums-of-squares are tracked separately and only assembled
    // into an EulerResidualNorms afterward.
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

    // Track the L2 norm of each conserved variable's flux residual computed
    // above, for callers monitoring convergence via residual().
    last_residual.rho = std::sqrt(sum_sq_rho);
    last_residual.rho_u = std::sqrt(sum_sq_rho_u);
    last_residual.rho_v = std::sqrt(sum_sq_rho_v);
    last_residual.E = std::sqrt(sum_sq_E);
}

// See EulerFVMSolver.h for the input/output contract.
void EulerFVMSolver::run(int total_steps)
{
    for (int t = 0; t < total_steps; ++t) {
        step();
    }
    std::cout << "Simulation completed across " << total_steps << " steps.\n";
}

// See EulerFVMSolver.h for the input/output contract.
//
// Methodology: for each face, the fastest signal speed crossing it is
// max(|Vn|+c) over the cell(s) touching it (boundary faces reuse the
// interior cell's own state on both sides -- cheap and sufficient for a
// stability estimate; the real ghost state is only needed for flux accuracy
// in step()). Dividing a characteristic cell length (mesh length units) by
// that speed (mesh length/time units) and scaling by the dimensionless CFL
// number gives a per-face stable dt, in time units; the smallest such dt
// over the whole mesh is the step size that keeps every face stable.
double EulerFVMSolver::compute_dt() const
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

        double length;
        if (face.cell_right != -1) {
            const EulerState& U_R = U[face.cell_right];
            double u_R = U_R.rho_u / U_R.rho;
            double v_R = U_R.rho_v / U_R.rho;
            double Vn_R = u_R * face.nx + v_R * face.ny;
            Smax = std::max(Smax, std::fabs(Vn_R) + sound_speed(U_R, gamma));
            length = std::min(mesh.cells[face.cell_left].volume, mesh.cells[face.cell_right].volume) / face.area;
        } else {
            length = mesh.cells[face.cell_left].volume / face.area;
        }

        dt = std::min(dt, cfl * length / Smax);
    }

    return dt;
}

// See EulerFVMSolver.h for the input/output contract.
//
// Methodology:
//   - Wall: mirror the normal velocity component (Vn_ghost = -Vn_L) while
//     keeping density and pressure equal to the interior cell, so the
//     Rusanov flux across the face carries no net normal mass/energy flux.
//   - Farfield: the prescribed state is returned directly, already converted
//     to conserved form once (at BC-setup time, not per-step).
//   - Outflow: zero-gradient extrapolation -- the ghost state equals the
//     interior cell's own state.
EulerState EulerFVMSolver::ghost_state(const Face& face, const EulerState& U_L) const
{
    const EulerBoundaryCondition& bc = bcs.at(face.patch_id);

    if (bc.type == EulerBoundaryType::Farfield) {
        return bc.farfield_state;
    }
    if (bc.type == EulerBoundaryType::Outflow) {
        return U_L;
    }

    // Wall
    double u_L = U_L.rho_u / U_L.rho;
    double v_L = U_L.rho_v / U_L.rho;
    double Vn = u_L * face.nx + v_L * face.ny;
    double u_ghost = u_L - 2.0 * Vn * face.nx;
    double v_ghost = v_L - 2.0 * Vn * face.ny;
    double p_L = pressure(U_L, gamma);
    return from_primitive(U_L.rho, u_ghost, v_ghost, p_L, gamma);
}
