// SPDX-License-Identifier: GPL-3.0-only
#ifndef EULERFVMSOLVER_H_INCLUDED
#define EULERFVMSOLVER_H_INCLUDED

#include <cassert>
#include <vector>

#include "UnstructuredMesh.h"
#include "EulerState.h"

// Boundary condition kind assigned to a patch for the Euler solver. Separate
// from the diffusion solver's BoundaryType (UnstructuredMesh.h) since the two
// solvers need different associated data per patch (a 4-component farfield
// state here, vs. a single scalar there) -- kept mesh-agnostic so
// UnstructuredMesh/BoundaryPatch stay shared, unmodified infrastructure for
// both solvers.
enum class EulerBoundaryType : int {
    Wall = 0,     // slip / no-penetration: reflect the normal velocity component
    Farfield = 1, // fixed prescribed state (farfield_state)
    Outflow = 2   // zero-gradient extrapolation of the interior cell's state
};

// Per-patch boundary condition for the Euler solver, indexed in parallel with
// UnstructuredMesh::patches (i.e. by Face::patch_id).
struct EulerBoundaryCondition {
    EulerBoundaryType type = EulerBoundaryType::Wall;
    EulerState farfield_state; // Prescribed conserved state; only meaningful when type == Farfield.
                                 // Units match EulerState (mesh-consistent mass/length/time units).
};

// Selects how EulerFVMSolver seeds its initial conserved state field.
enum class EulerICMode {
    Freestream, // uniform primitive state everywhere
    TwoRegion   // a left/right primitive state split at x = x0 (e.g. Sod shock tube)
};

// Selects which numerical flux (see RusanovFlux.h/HllcFlux.h/
// ExactRiemannFlux.h) EulerFVMSolver::step() uses at every face. Constant
// for the whole run.
enum class NumericalFluxScheme {
    Rusanov, // local Lax-Friedrichs; cheapest and most robust, most dissipative
    HLLC,    // restores the contact wave; sharper contact/shear resolution
    Exact    // exact Riemann (Godunov) solver; highest accuracy, most expensive
};

// L2 norms of the per-step flux residual, one per conserved variable (see
// EulerState.h) since they have very different physical scales and a single
// combined number would be misleading. Units match each variable's flux
// (e.g. rho's is a mass flux, un-normalized by dt/volume).
struct EulerResidualNorms {
    double rho = 0.0, rho_u = 0.0, rho_v = 0.0, E = 0.0;
};

// Initial condition for the Euler solver, in primitive variables (converted
// to conserved EulerState internally, using the solver's gamma). Units are
// generic/mesh-consistent, matching EulerState.h's convention: rho is a
// density (mass/length^2 for this 2D solver), u/v are velocity components
// (length/time), p is a pressure, x0 is in mesh length units.
struct EulerInitialCondition {
    EulerICMode mode = EulerICMode::Freestream;

    double rho = 1.0, u = 0.0, v = 0.0, p = 1.0; // Freestream: state everywhere

    double rho_l = 1.0, u_l = 0.0, v_l = 0.0, p_l = 1.0;     // TwoRegion: state where x_centroid < x0
    double rho_r = 0.125, u_r = 0.0, v_r = 0.0, p_r = 0.1;   // TwoRegion: state where x_centroid >= x0
    double x0 = 0.0;                                          // TwoRegion: split location, in mesh length units
};

// Explicit finite-volume solver for the 2D compressible Euler equations on an
// arbitrary unstructured polygon mesh, using a selectable numerical flux
// (see NumericalFluxScheme/RiemannSolvers.h) and a CFL-limited adaptive time
// step.
//
// Methodology: for each cell, integrate the Euler flux over the cell's
// bounding faces (divergence theorem), then advance the conserved state one
// explicit forward-Euler step: U_new = U_old + dt/volume * (net flux in).
// Face fluxes are computed once per face and applied with opposite sign to
// its two neighboring cells, which keeps the scheme conservative. Unlike the
// scalar diffusion solver, dt is not fixed: it is recomputed every step from
// a user-specified CFL number and the fastest wave speed currently present in
// the flow, since the stable time step for a hyperbolic system shrinks or
// grows as the solution evolves.
class EulerFVMSolver {
public:
    // Input:
    //   input_mesh          - mesh topology/geometry to solve on (copied; its
    //                          compute_geometry() is called here, so the
    //                          caller does not need to call it first)
    //   boundary_conditions - per-patch BC, indexed like mesh.patches
    //   gamma               - ratio of specific heats, dimensionless (1.4 for air)
    //   cfl                 - CFL number, dimensionless, used to size each
    //                          step's dt (e.g. 0.5; must be <= ~1 for stability,
    //                          smaller is more conservative)
    //   initial_condition   - primitive-variable initial condition, in
    //                          mesh-consistent units (see EulerInitialCondition above)
    //   flux_scheme         - which numerical flux step() uses at every face
    //                          (see NumericalFluxScheme above)
    //   exact_riemann_tol   - relative pressure tolerance for the exact
    //                          Riemann solver's Newton-Raphson iteration
    //                          (see ExactRiemannFlux.h); ignored unless
    //                          flux_scheme == Exact
    //   exact_riemann_max_iter - iteration cap for that same Newton-Raphson
    //                          solve; ignored unless flux_scheme == Exact
    // Output: constructs a solver with U sized to mesh.cells.size() and the
    //         initial condition above applied (converted to conserved form
    //         via from_primitive(), see EulerState.h).
    EulerFVMSolver(const UnstructuredMesh& input_mesh,
                    const std::vector<EulerBoundaryCondition>& boundary_conditions,
                    double gamma, double cfl,
                    const EulerInitialCondition& initial_condition,
                    NumericalFluxScheme flux_scheme,
                    double exact_riemann_tol, int exact_riemann_max_iter);

    // Advances the solution by one explicit step. The step size is recomputed
    // every call from compute_dt() (see below), not fixed, since it must
    // track the fastest wave speed currently present in the flow.
    // Input:  none (uses current U/mesh state)
    // Output: U is updated in place to its value at t + dt
    void step();

    // Repeatedly calls step() to advance the solution.
    // Input:  total_steps - number of time steps to advance (each of adaptive size)
    // Output: U is advanced total_steps steps; prints a completion message
    void run(int total_steps);

    // Input:  none
    // Returns: the current conserved-state field, one value per mesh cell
    //          (same order as UnstructuredMesh::cells; units as in EulerState.h)
    const std::vector<EulerState>& field() const { return U; }

    // Input:  none
    // Returns: the ratio of specific heats (dimensionless) this solver was
    //          constructed with
    double gas_gamma() const { return gamma; }

    // Input:  none
    // Returns: the L2 norm of the most recent step()'s per-cell flux residual
    //          (the net-flux-in vector, before the dt/volume update), one
    //          value per conserved variable; all zero before the first
    //          step() call.
    const EulerResidualNorms& residual() const { return last_residual; }

    // Overwrites the current conserved-state field, e.g. to inject a
    // checkpoint-resumed state after construction (bypassing the
    // constructor's initial condition).
    // Input:  new_U - must be sized exactly mesh.cells.size(); caller is
    //         responsible for that match (already validated against the
    //         checkpoint's stored cell count one call earlier)
    // Output: U is replaced; residual() resets to all-zero (no step() has
    //         run against this field yet)
    void set_field(const std::vector<EulerState>& new_U) {
        assert(new_U.size() == U.size());
        U = new_U;
        last_residual = EulerResidualNorms{};
    }

private:
    UnstructuredMesh mesh;
    std::vector<EulerState> U;                  // Conserved state field: one entry per mesh.cells
    std::vector<EulerBoundaryCondition> bcs;    // Per-patch BC, indexed like mesh.patches
    double gamma; // Ratio of specific heats, dimensionless
    double cfl;   // CFL number controlling the adaptive time step, dimensionless
    NumericalFluxScheme flux_scheme; // Which numerical flux step() uses at every face
    double exact_riemann_tol;    // Exact Riemann solver's Newton-Raphson relative pressure tolerance
    int exact_riemann_max_iter;  // Exact Riemann solver's Newton-Raphson iteration cap
    EulerResidualNorms last_residual; // L2 norms of the most recent step()'s flux residual

    // Computes a stable explicit time step from the CFL condition, scanning
    // every face for the fastest wave speed (|Vn| + speed of sound) crossing
    // it and a characteristic cell length, then taking the most restrictive
    // (smallest) resulting dt across the whole mesh.
    //
    // Methodology: for each face, a candidate dt is cfl * length / Smax,
    // where Smax is the largest |Vn|+c seen on either side of the face and
    // length is a characteristic cell size (volume/face.area); the true dt is
    // the minimum candidate over every face, so the whole mesh advances at
    // the pace its locally-fastest wave allows.
    // Input:  none (uses current U/mesh state)
    // Returns: dt, in time units
    double compute_dt() const;

    // Builds the state on the far side of a boundary face (the "ghost" state
    // used in place of a real neighbor cell), according to that face's
    // patch's EulerBoundaryType. See EulerFVMSolver.cpp for the per-type
    // methodology (wall reflection / farfield / zero-gradient outflow).
    // Input:  face  - a boundary face (face.cell_right == -1); its patch_id
    //                 selects the EulerBoundaryCondition to apply
    //         U_L   - current conserved state of face.cell_left
    // Returns: the ghost state to pair with U_L in the Rusanov flux
    EulerState ghost_state(const Face& face, const EulerState& U_L) const;
};

#endif // EULERFVMSOLVER_H_INCLUDED
