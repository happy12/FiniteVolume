// SPDX-License-Identifier: GPL-3.0-only
#include "UnstructuredFVMSolver.h"

// See UnstructuredFVMSolver.h for the input/output contract.
UnstructuredFVMSolver::UnstructuredFVMSolver(const UnstructuredMesh& input_mesh, double alpha, double dt,
                                              double initial_value, double initial_radius)
    : mesh(input_mesh), alpha(alpha), dt(dt)
{
    mesh.compute_geometry();
    phi.resize(mesh.cells.size(), 0.0);

    // Initial Condition: set phi = initial_value in a disc of radius
    // initial_radius centered on the origin, and phi = 0 everywhere else.
    for (size_t i = 0; i < mesh.cells.size(); ++i) {
        if (std::sqrt(mesh.cells[i].x_centroid * mesh.cells[i].x_centroid +
                      mesh.cells[i].y_centroid * mesh.cells[i].y_centroid) < initial_radius) {
            phi[i] = initial_value;
        }
    }
}
UnstructuredFVMSolver::~UnstructuredFVMSolver()
{
    //dtor
}

// See UnstructuredFVMSolver.h for the input/output contract.
//
// Methodology: a single explicit forward-Euler finite-volume step, split
// into two OpenMP-parallel passes to avoid a scatter-add race (two faces
// writing to the same cell's residual concurrently):
//   1. For every face (parallel, no shared writes -- each face only writes
//      its own face_flux[i] entry), compute the diffusive flux crossing it
//      (Fick's law: flux = -alpha * dPhi/dn * area, using a two-point
//      central difference for the normal gradient between the two cells it
//      connects, or the boundary condition if it's a boundary face).
//   2. For every cell (parallel, no shared writes -- each cell only writes
//      its own phi[c]), sum the fluxes of its own bounding faces (via
//      Cell::faces), with a sign of -1 if the cell is that face's cell_left
//      and +1 if it's cell_right -- equivalent to face-major scatter-add,
//      just re-ordered to be race-free -- then integrate:
//      phi_new = phi_old + dt/volume * residual (explicit Euler).
// Stability requires dt small enough relative to alpha and the mesh spacing
// (the classic explicit-diffusion CFL-type limit); this is not checked here.
void UnstructuredFVMSolver::step()
{
    std::vector<double> face_flux(mesh.faces.size());

    // --- Pass 1: per-face flux (parallel over faces) ---
    #pragma omp parallel for
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        double phi_L = phi[face.cell_left];

        // Track geometry paths depending on whether the face is internal or a boundary wall
        if (face.cell_right != -1) {
            double phi_R = phi[face.cell_right];

            double dx = mesh.cells[face.cell_right].x_centroid - mesh.cells[face.cell_left].x_centroid;
            double dy = mesh.cells[face.cell_right].y_centroid - mesh.cells[face.cell_left].y_centroid;
            double d_LR = std::sqrt(dx*dx + dy*dy);

            // Normal derivative approximation: dPhi/dn ~= (Phi_R - Phi_L) / d_LR
            double grad_phi_n = (phi_R - phi_L) / d_LR;

            // Diffusive Flux = -Alpha * Spatial Gradient * Area of interface
            face_flux[i] = -alpha * grad_phi_n * face.area;
        } else {
            face_flux[i] = boundary_flux(face, phi_L);
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

        // New State = Old State + dt/Volume * Net Face Flux
        phi[c] += (dt / mesh.cells[c].volume) * residual_c;
    }

    // Track the L2 norm of the flux residual computed above, for callers
    // monitoring convergence via residual_norm().
    last_residual = std::sqrt(sum_sq);
}

// See UnstructuredFVMSolver.h for the input/output contract.
void UnstructuredFVMSolver::run(int total_steps)
{
    for (int t = 0; t < total_steps; ++t) {
        step();
    }
    std::cout << "Simulation completed across " << total_steps << " steps.\n";
}

// See UnstructuredFVMSolver.h for the input/output contract.
//
// Methodology:
//   - Neumann: the flux itself is prescribed (patch.value = -alpha * dPhi/dn,
//     i.e. already in flux form), so it's applied directly, scaled by face area.
//   - Dirichlet: phi is pinned to patch.value AT THE FACE MIDPOINT, so the
//     normal gradient is approximated as a one-point difference between the
//     face midpoint and cell_left's centroid (rather than the two-cell
//     centroid-to-centroid difference used for internal faces, since there is
//     no neighbor cell on the other side of a boundary face).
double UnstructuredFVMSolver::boundary_flux(const Face& face, double phi_L) const
{
    const BoundaryPatch& patch = mesh.patches.at(face.patch_id);

    if (patch.type == BoundaryType::Neumann) {
        // Directly prescribed diffusive flux (patch.value = -alpha * dPhi/dn)
        return patch.value * face.area;
    }

    // Dirichlet: phi = patch.value enforced at the face midpoint
    double dx = face.x_mid - mesh.cells[face.cell_left].x_centroid;
    double dy = face.y_mid - mesh.cells[face.cell_left].y_centroid;
    double d_LR = std::sqrt(dx*dx + dy*dy);

    double grad_phi_n = (patch.value - phi_L) / d_LR;
    return -alpha * grad_phi_n * face.area;
}
