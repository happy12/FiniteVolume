// SPDX-License-Identifier: GPL-3.0-only
#include "AdvectionDiffusionFVMSolver.h"

#include <cmath>
#include <iostream>

namespace {

// compute_geometry() must run before GradientCalculator's constructor
// inspects the mesh (it precomputes Least-Squares' per-cell matrix from
// cell/face geometry), so it can't be deferred to the constructor body the
// way UnstructuredFVMSolver/EulerFVMSolver do it -- gradient_calc is a
// member initialized in the same initializer list as mesh, before the body
// runs.
UnstructuredMesh mesh_with_geometry(UnstructuredMesh m) {
    m.compute_geometry();
    return m;
}

} // namespace

// See AdvectionDiffusionFVMSolver.h for the input/output contract.
AdvectionDiffusionFVMSolver::AdvectionDiffusionFVMSolver(const UnstructuredMesh& input_mesh, double alpha,
                                                           double u_adv, double v_adv, double dt,
                                                           GradientScheme gradient_scheme, double initial_value,
                                                           double initial_radius)
    : mesh(mesh_with_geometry(input_mesh)), alpha(alpha), u_adv(u_adv), v_adv(v_adv), dt(dt),
      gradient_calc(mesh, gradient_scheme)
{
    phi.resize(mesh.cells.size(), 0.0);

    // Initial Condition: set phi = initial_value in a disc of radius
    // initial_radius centered on the origin, and phi = 0 everywhere else.
    for (size_t i = 0; i < mesh.cells.size(); ++i) {
        if (std::hypot(mesh.cells[i].x_centroid, mesh.cells[i].y_centroid) < initial_radius) {
            phi[i] = initial_value;
        }
    }
}

AdvectionDiffusionFVMSolver::~AdvectionDiffusionFVMSolver()
{
    //dtor
}

// See AdvectionDiffusionFVMSolver.h for the input/output contract.
std::vector<double> AdvectionDiffusionFVMSolver::build_boundary_phi() const
{
    std::vector<double> boundary_phi(mesh.faces.size(), 0.0);
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        if (face.cell_right != -1) continue;
        const BoundaryPatch& patch = mesh.patches.at(face.patch_id);
        boundary_phi[i] = (patch.type == BoundaryType::Dirichlet) ? patch.value : phi[face.cell_left];
    }
    return boundary_phi;
}

// See AdvectionDiffusionFVMSolver.h for the input/output contract.
//
// Methodology: a single explicit forward-Euler finite-volume step. Cell
// gradients of phi are reconstructed once up front (needed by every
// boundary/internal face's diffusive term), then the flux assembly follows
// the same two-OpenMP-pass pattern as UnstructuredFVMSolver::step() (pass 1:
// per-face flux, no shared writes; pass 2: per-cell residual + integration,
// no shared writes) to avoid a scatter-add race between the two cells
// sharing a face.
void AdvectionDiffusionFVMSolver::step()
{
    std::vector<double> boundary_phi = build_boundary_phi();
    std::vector<Gradient2> grad = gradient_calc.compute(mesh, phi, boundary_phi);

    std::vector<double> face_flux(mesh.faces.size());

    // --- Pass 1: per-face flux (parallel over faces) ---
    #pragma omp parallel for
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        double phi_L = phi[face.cell_left];

        if (face.cell_right != -1) {
            double phi_R = phi[face.cell_right];

            FaceGradient fg = face_gradient(mesh, face, phi_L, phi_R, grad[face.cell_left],
                                             grad[face.cell_right], 0.0);
            double diffusive = -alpha * fg.dphidn * face.area;

            double Vn = u_adv * face.nx + v_adv * face.ny;
            double phi_upwind = (Vn >= 0.0) ? phi_L : phi_R;
            double advective = Vn * phi_upwind * face.area;

            face_flux[i] = advective + diffusive;
        } else {
            face_flux[i] = boundary_flux(face, phi_L, grad[face.cell_left]);
        }
    }

    // --- Pass 2: per-cell residual + integration (parallel over cells) ---
    double sum_sq = 0.0;
    #pragma omp parallel for reduction(+:sum_sq)
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        double residual_c = 0.0;
        for (int face_idx : mesh.cells[c].faces) {
            const Face& face = mesh.faces[face_idx];
            if ((size_t)face.cell_left == c) residual_c -= face_flux[face_idx];
            else residual_c += face_flux[face_idx];
        }

        sum_sq += residual_c * residual_c;
        phi[c] += (dt / mesh.cells[c].volume) * residual_c;
    }

    last_residual = std::sqrt(sum_sq);
}

// See AdvectionDiffusionFVMSolver.h for the input/output contract.
void AdvectionDiffusionFVMSolver::run(int total_steps)
{
    for (int t = 0; t < total_steps; ++t) {
        step();
    }
    std::cout << "Simulation completed across " << total_steps << " steps.\n";
}

// See AdvectionDiffusionFVMSolver.h for the input/output contract.
//
// Methodology:
//   - Neumann: the diffusive flux itself is prescribed directly (patch.value
//     = -alpha * dPhi/dn, i.e. already in flux form, same convention as
//     UnstructuredFVMSolver's boundary_flux()), bypassing gradient-based
//     reconstruction for this term entirely. The advective term still needs
//     an upwind value; since no exterior value is prescribed for Neumann,
//     the owning cell's own value is used both for outflow (correct) and,
//     as a zero-order approximation, for any inflow (no better option is
//     available without an exterior state).
//   - Dirichlet: phi = patch.value is pinned at the face midpoint, giving
//     both the diffusive term's boundary value (via face_gradient(),
//     treating the face midpoint as the "R" point -- see
//     GradientReconstruction.h) and the advective term's upwind value on
//     inflow.
double AdvectionDiffusionFVMSolver::boundary_flux(const Face& face, double phi_L, const Gradient2& grad_L) const
{
    const BoundaryPatch& patch = mesh.patches.at(face.patch_id);
    double Vn = u_adv * face.nx + v_adv * face.ny;

    if (patch.type == BoundaryType::Neumann) {
        double advective = Vn * phi_L * face.area;
        return advective + patch.value * face.area;
    }

    FaceGradient fg = face_gradient(mesh, face, phi_L, 0.0, grad_L, Gradient2{}, patch.value);
    double diffusive = -alpha * fg.dphidn * face.area;

    double phi_upwind = (Vn >= 0.0) ? phi_L : patch.value;
    double advective = Vn * phi_upwind * face.area;

    return advective + diffusive;
}
