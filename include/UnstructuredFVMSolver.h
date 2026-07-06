// SPDX-License-Identifier: GPL-3.0-only
#ifndef UNSTRUCTUREDFVMSOLVER_H
#define UNSTRUCTUREDFVMSOLVER_H

#include "UnstructuredMesh.h"
#include <cassert>
#include <numeric>

// Explicit finite-volume solver for the 2D scalar diffusion (heat/Fick's law)
// equation on an arbitrary unstructured polygon mesh:
//
//   dPhi/dt = alpha * Laplacian(Phi)
//
// Methodology: for each cell, integrate the diffusive flux -alpha * grad(Phi)
// over the cell's bounding faces (the divergence theorem turns the volume
// integral of the Laplacian into a sum of face fluxes), then advance phi one
// explicit forward-Euler step: phi_new = phi_old + dt/volume * (net flux in).
// Face fluxes are computed once per face and applied with opposite sign to
// its two neighboring cells, which guarantees the scheme is conservative
// (whatever leaves one cell through a face exactly enters its neighbor).
class UnstructuredFVMSolver {
public:
    // Input:
    //   input_mesh      - mesh topology/geometry to solve on (copied; its
    //                      compute_geometry() is called here, so the caller
    //                      does not need to call it first)
    //   alpha           - diffusion coefficient, in (mesh length units)^2 / (time units)
    //   dt              - time step size, in time units
    //   initial_value   - phi value assigned to cells within initial_radius of the origin
    //   initial_radius  - radius (from the origin) of the initial condition disc,
    //                      in mesh length units; cells outside it start at phi = 0
    // Output: constructs a solver with phi sized to mesh.cells.size() and the
    //         initial condition above applied.
    UnstructuredFVMSolver(const UnstructuredMesh& input_mesh, double alpha, double dt,
                           double initial_value, double initial_radius);
    ~UnstructuredFVMSolver();

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
    //          (the net-flux-in vector, before the dt/volume update) -- a
    //          measure of how far the discrete equations are from being
    //          satisfied; 0.0 before the first step() call. Units match phi's
    //          flux (phi * mesh length units / time unit, since it's an
    //          un-normalized face-integrated flux).
    double residual_norm() const { return last_residual; }

    // Overwrites the current state field, e.g. to inject a checkpoint-resumed
    // state after construction (bypassing the constructor's initial-condition
    // disc).
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
    double dt;                // Time step size, in time units
    double last_residual = 0.0; // L2 norm of the most recent step()'s flux residual

    // Evaluates the diffusive flux crossing a boundary face, applying its patch's BC.
    // Input:
    //   face  - a boundary face (face.cell_right == -1); its patch_id selects
    //           the BoundaryPatch (type + value) to apply
    //   phi_L - current phi of face.cell_left (the face's owning cell)
    // Returns: the flux crossing the face, in the same convention as the
    //          internal-face flux in step() (positive = leaving cell_left),
    //          so it can be subtracted from cell_left's residual the same way.
    double boundary_flux(const Face& face, double phi_L) const;

};

#endif // UNSTRUCTUREDFVMSOLVER_H
