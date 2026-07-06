// SPDX-License-Identifier: GPL-3.0-only
#ifndef ADVECTIONDIFFUSIONFVMSOLVER_H_INCLUDED
#define ADVECTIONDIFFUSIONFVMSOLVER_H_INCLUDED

#include <cassert>
#include <vector>

#include "GradientReconstruction.h"
#include "UnstructuredMesh.h"

// Explicit finite-volume solver for the 2D advection-diffusion equation of a
// passive scalar on an arbitrary unstructured polygon mesh:
//
//   dPhi/dt + div(u * Phi) = alpha * Laplacian(Phi)
//
// where u = (u_adv, v_adv) is a uniform, prescribed (not solved-for)
// advection velocity. Reuses UnstructuredMesh's BoundaryType/BoundaryPatch
// (Dirichlet/Neumann on a scalar), since that BC semantics is identical to
// the pure-diffusion solver's.
//
// Methodology: for each cell, integrate the net advective + diffusive flux
// over the cell's bounding faces (divergence theorem), then advance phi one
// explicit forward-Euler step: phi_new = phi_old + dt/volume * (net flux
// in). Per face:
//   - Advective flux = (u . n) * phi_upwind * area, first-order upwind
//     (phi_upwind is whichever side the face-normal velocity flows FROM).
//   - Diffusive flux = -alpha * dphidn * area, where dphidn is the
//     non-orthogonality-corrected face-normal derivative from
//     GradientReconstruction.h's face_gradient() -- a naive two-point
//     difference has a real, large error on a non-orthogonal mesh (see
//     docs/navier-stokes-tracker.md, Phase 1). Getting dphidn requires a
//     cell gradient of phi in every cell first, reconstructed once per step
//     via GradientCalculator (see gradient_scheme).
// Face fluxes are computed once per face and applied with opposite sign to
// its two neighboring cells, which keeps the scheme conservative, following
// the same two-OpenMP-pass pattern (no scatter-add race) as
// UnstructuredFVMSolver/EulerFVMSolver.
class AdvectionDiffusionFVMSolver {
public:
    // Input:
    //   input_mesh      - mesh topology/geometry to solve on (copied; its
    //                      compute_geometry() is called here, so the caller
    //                      does not need to call it first)
    //   alpha           - diffusion coefficient, in (mesh length units)^2 / (time units)
    //   u_adv, v_adv    - uniform advection velocity components, in mesh
    //                      length units / time unit
    //   dt              - time step size, in time units (fixed; stability
    //                      against both the advective CFL limit and the
    //                      explicit-diffusion limit is the caller's
    //                      responsibility, same as UnstructuredFVMSolver's dt)
    //   gradient_scheme - which GradientCalculator scheme reconstructs cell
    //                      gradients of phi every step (see GradientReconstruction.h)
    //   initial_value   - phi value assigned to cells within initial_radius of the origin
    //   initial_radius  - radius (from the origin) of the initial condition disc,
    //                      in mesh length units; cells outside it start at phi = 0
    // Output: constructs a solver with phi sized to mesh.cells.size() and the
    //         initial condition above applied.
    AdvectionDiffusionFVMSolver(const UnstructuredMesh& input_mesh, double alpha, double u_adv, double v_adv,
                                 double dt, GradientScheme gradient_scheme, double initial_value,
                                 double initial_radius);
    ~AdvectionDiffusionFVMSolver();

    // Advances the solution by a single explicit forward-Euler time step of size dt.
    // Input:  none (uses current phi/mesh state)
    // Output: phi is updated in place to its value at t + dt
    void step();

    // Repeatedly calls step() to advance the solution.
    // Input:  total_steps - number of time steps to advance (each of size dt)
    // Output: phi is advanced to t + total_steps * dt; prints a completion message
    void run(int total_steps);

    // Returns the current state field, one value per mesh cell (same order
    // and units as UnstructuredMesh::cells).
    const std::vector<double>& field() const { return phi; }

    // Input:  none
    // Returns: the L2 norm of the most recent step()'s per-cell flux residual
    //          (the net-flux-in vector, before the dt/volume update); 0.0
    //          before the first step() call.
    double residual_norm() const { return last_residual; }

    // Overwrites the current state field, e.g. to inject a checkpoint-resumed
    // state after construction (bypassing the constructor's initial-condition disc).
    // Input:  new_phi - must be sized exactly mesh.cells.size(); caller is
    //         responsible for that match (already validated against the
    //         checkpoint's stored cell count one call earlier)
    // Output: phi is replaced; residual_norm() resets to 0.0 (no step() has
    //         run against this field yet)
    void set_field(const std::vector<double>& new_phi) {
        assert(new_phi.size() == phi.size());
        phi = new_phi;
        last_residual = 0.0;
    }

private:
    UnstructuredMesh mesh;
    std::vector<double> phi; // State field: size maps exactly to mesh.cells.size()
    double alpha;             // Diffusion coefficient, in (mesh length units)^2 / (time units)
    double u_adv, v_adv;      // Uniform advection velocity components
    double dt;                // Time step size, in time units
    GradientCalculator gradient_calc; // Reconstructs cell gradients of phi every step
    double last_residual = 0.0; // L2 norm of the most recent step()'s flux residual

    // Builds the mesh.faces-sized boundary-value array GradientCalculator's
    // compute() needs: the prescribed value for a Dirichlet boundary face,
    // or a zero-order (equal to the owning cell's own current value)
    // extrapolation for a Neumann one, since only a flux (not a value) is
    // prescribed there. Entries for internal faces are left at 0.0 (unused).
    // Input:  none (uses the current phi)
    // Returns: a mesh.faces-sized array as described above
    std::vector<double> build_boundary_phi() const;

    // Evaluates the combined advective + diffusive flux crossing a boundary
    // face, applying its patch's BC.
    // Input:
    //   face   - a boundary face (face.cell_right == -1); its patch_id
    //            selects the BoundaryPatch (type + value) to apply
    //   phi_L  - current phi of face.cell_left (the face's owning cell)
    //   grad_L - reconstructed cell gradient of face.cell_left
    // Returns: the flux crossing the face, in the same convention as the
    //          internal-face flux in step() (positive = leaving cell_left)
    double boundary_flux(const Face& face, double phi_L, const Gradient2& grad_L) const;
};

#endif // ADVECTIONDIFFUSIONFVMSOLVER_H_INCLUDED
