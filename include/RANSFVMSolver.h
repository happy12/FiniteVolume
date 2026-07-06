// SPDX-License-Identifier: GPL-3.0-only
#ifndef RANSFVMSOLVER_H_INCLUDED
#define RANSFVMSOLVER_H_INCLUDED

#include <cassert>
#include <vector>

#include "EulerFVMSolver.h" // reuses EulerICMode/EulerInitialCondition/NumericalFluxScheme/EulerResidualNorms
#include "EulerState.h"
#include "GradientReconstruction.h"
#include "NavierStokesFVMSolver.h" // reuses NSBoundaryType/NSBoundaryCondition unchanged
#include "SpalartAllmaras.h"
#include "UnstructuredMesh.h"
#include "WallTraction.h"

// Per-patch boundary condition for the RANS solver, indexed in parallel with
// UnstructuredMesh::patches. Wraps NSBoundaryCondition unchanged (a no-slip
// wall/farfield/outflow means exactly the same thing for the mean-flow
// equations here as it does for NavierStokesFVMSolver) and adds only what
// the nut transport equation needs beyond that: a prescribed freestream nut
// for a Farfield boundary. NoSlipWall's nut boundary value is NOT
// configurable here -- it is exactly 0 by the SA model's own definition
// (see build_boundary_nut() in RANSFVMSolver.cpp), and Outflow's is a
// zero-order extrapolation, matching every other field's Outflow treatment
// in this codebase.
struct RANSBoundaryCondition {
    NSBoundaryCondition ns;
    double farfield_nut = 0.0; // prescribed nut when ns.type == NSBoundaryType::Farfield; SA is documented as
                                 // sensitive to this choice (see docs/archive/rans-spalart-allmaras-tracker.md Phase 4's
                                 // "known setup difficulty")
};

// Explicit finite-volume solver for the 2D compressible RANS equations,
// closed with the Spalart-Allmaras (SA-noft2) one-equation turbulence model
// (see SpalartAllmaras.h): NavierStokesFVMSolver's exact mean-flow equations
// (inviscid flux, viscous stress tensor, Fourier heat conduction), plus one
// extra transported scalar nut ("nu-tilde") advected by the coupled mean
// flow velocity, with molecular viscosity/conductivity in the mean-flow
// viscous terms replaced by their turbulence-inclusive effective values:
//
//   mu_eff = mu + rho*nu_t                    (nu_t = nut * fv1(nut, mu/rho))
//   k_eff  = cp*(mu/Pr + rho*nu_t/Pr_t)        (Pr_t = turbulent Prandtl number)
//
// Per the 2026-07-05 architecture decision recorded in
// docs/archive/rans-spalart-allmaras-tracker.md, this is deliberately a THIRD
// solver class (after EulerFVMSolver/NavierStokesFVMSolver), duplicating
// NavierStokesFVMSolver's inviscid/viscous flux, ghost-state, and
// compute_dt() structure rather than adding a "turbulence_model" toggle to
// it -- consistent with this project's "one class per equation set" pattern.
//
// nut's own transport equation (structurally an advection-diffusion
// equation, like AdvectionDiffusionFVMSolver, but advected by the REAL
// coupled flow velocity rather than a prescribed one):
//
//   dnut/dt + div(u*nut) = production - destruction
//                        + (1/sigma)*div((nu+nut)*grad(nut))
//                        + (Cb2/sigma)*|grad(nut)|^2
//
// production/destruction/the (Cb2/sigma)*|grad(nut)|^2 cross-diffusion term
// are volumetric (per-cell) source terms from compute_sa_source_terms()
// (SpalartAllmaras.h, verified in isolation in this tracker's Phase 2); the
// div((nu+nut)*grad(nut)) part is a genuine face flux, assembled the same
// way as every other diffusive term in this codebase: a non-orthogonality-
// corrected face-normal derivative (GradientReconstruction.h's
// face_gradient()) times a face-interpolated (arithmetic-mean) diffusivity.
// nut's own advective flux uses the SAME ghost/interior state pairing as
// the inviscid mean-flow flux (so a Farfield boundary's real prescribed
// velocity, not an extrapolated one, sets the convecting velocity there),
// first-order upwind on the resulting face-normal velocity.
//
// nut's diffusive flux (unlike its advective flux) is forced to exactly
// zero at Farfield/Outflow boundaries, mirroring NavierStokesFVMSolver's own
// convention for the mean-flow viscous flux at open boundaries (assumed
// placed far enough from a wall that this doesn't matter) -- see
// has_viscous_term in RANSFVMSolver.cpp.
class RANSFVMSolver {
public:
    // Input:
    //   input_mesh          - mesh topology/geometry to solve on (copied;
    //                          its compute_geometry() is called here)
    //   boundary_conditions - per-patch BC, indexed like mesh.patches
    //   gamma               - ratio of specific heats, dimensionless (1.4 for air)
    //   gas_constant        - specific gas constant R (see EulerState.h's temperature())
    //   mu                  - molecular dynamic viscosity, mesh-consistent units
    //   prandtl             - molecular Prandtl number, dimensionless (0.72 for air)
    //   prandtl_t           - turbulent Prandtl number, dimensionless (~0.9 for air)
    //   cfl                 - CFL number, dimensionless, used to size each
    //                          step's dt (see compute_dt())
    //   initial_condition   - primitive-variable mean-flow initial condition
    //                          (see EulerInitialCondition, EulerFVMSolver.h)
    //   initial_nut         - uniform initial/freestream nut value applied to
    //                          every cell (SA is documented as sensitive to
    //                          this choice; see RANSBoundaryCondition::farfield_nut)
    //   flux_scheme         - which numerical flux the inviscid term uses at
    //                          every face (see NumericalFluxScheme, EulerFVMSolver.h)
    //   gradient_scheme     - which GradientCalculator scheme reconstructs
    //                          cell gradients of u/v/T/nut every step
    //   exact_riemann_tol, exact_riemann_max_iter - ignored unless flux_scheme == Exact
    //   sa_constants        - SA-noft2 model constants (see SpalartAllmaras.h);
    //                          defaults to the standard set
    // Output: constructs a solver with U/nut sized to mesh.cells.size() and
    //         the initial conditions above applied; also precomputes the
    //         wall-distance field (WallDistance.h) once from the faces whose
    //         patch's NSBoundaryType is NoSlipWall, since it depends only on
    //         mesh geometry and BC assignment, not on the flow field -- same
    //         reasoning as GradientCalculator's Least-Squares matrix.
    RANSFVMSolver(const UnstructuredMesh& input_mesh, const std::vector<RANSBoundaryCondition>& boundary_conditions,
                   double gamma, double gas_constant, double mu, double prandtl, double prandtl_t, double cfl,
                   const EulerInitialCondition& initial_condition, double initial_nut, NumericalFluxScheme flux_scheme,
                   GradientScheme gradient_scheme, double exact_riemann_tol, int exact_riemann_max_iter,
                   const SAModelConstants& sa_constants = SAModelConstants{});

    // Advances the solution by one explicit step, at a step size recomputed
    // every call from compute_dt() (the inviscid CFL limit and the
    // explicit-viscous-diffusion stability limit, the latter now using
    // nu + nu_t rather than just molecular nu -- see compute_dt()).
    void step();

    // Repeatedly calls step() to advance the solution.
    void run(int total_steps);

    const std::vector<EulerState>& field() const { return U; }
    const std::vector<double>& nut_field() const { return nut; }
    double gas_gamma() const { return gamma; }
    const EulerResidualNorms& residual() const { return last_residual; }
    double nut_residual() const { return last_nut_residual; }

    // Overwrites the current conserved-state AND nut fields, e.g. to inject a
    // checkpoint-resumed state after construction -- see
    // NavierStokesFVMSolver::set_field() for the identical mean-flow-only
    // precedent; RANS needs both fields resumed together since they're
    // coupled through the same step().
    void set_field(const std::vector<EulerState>& new_U, const std::vector<double>& new_nut) {
        assert(new_U.size() == U.size());
        assert(new_nut.size() == nut.size());
        U = new_U;
        nut = new_nut;
        last_residual = EulerResidualNorms{};
        last_nut_residual = 0.0;
    }

    // Wall traction diagnostics (see WallTraction.h) at every face in
    // 'wall_faces', from the CURRENT flow field -- identical on-demand,
    // nothing-retained-from-step() contract as
    // NavierStokesFVMSolver::compute_wall_traction_samples(), except
    // effective_viscosity is mu + rho*nu_t per cell (this solver's turbulence
    // closure), not uniform mu.
    //
    // Input:  wall_faces - indices into mesh.faces to sample; every entry
    //         should be a boundary face whose patch's NSBoundaryType is
    //         NoSlipWall (not checked here -- see WallTraction.h)
    // Returns: one WallFaceSample per entry of 'wall_faces', in the same order
    std::vector<WallFaceSample> compute_wall_traction_samples(const std::vector<int>& wall_faces) const;

    // Boundary-layer thickness diagnostics (see
    // compute_boundary_layer_profiles(), WallTraction.h) at every face in
    // 'wall_faces', from the CURRENT flow field -- identical contract and
    // parameters as NavierStokesFVMSolver::compute_boundary_layer_profile_samples().
    std::vector<BoundaryLayerProfile> compute_boundary_layer_profile_samples(const std::vector<int>& wall_faces,
                                                                              double u_edge,
                                                                              int max_cells_per_march) const;

    // Point-location alternative to compute_boundary_layer_profile_samples()
    // above -- identical contract and parameters as
    // NavierStokesFVMSolver::compute_boundary_layer_profile_samples_point_location().
    std::vector<BoundaryLayerProfile> compute_boundary_layer_profile_samples_point_location(
        const std::vector<int>& wall_faces, double u_edge, double max_distance, int n_samples) const;

private:
    UnstructuredMesh mesh;
    std::vector<EulerState> U;
    std::vector<double> nut;
    std::vector<RANSBoundaryCondition> bcs;
    double gamma, gas_constant, mu, prandtl, prandtl_t, cfl;
    NumericalFluxScheme flux_scheme;
    GradientCalculator gradient_calc; // Reconstructs cell gradients of u, v, T, and nut every step
    double exact_riemann_tol;
    int exact_riemann_max_iter;
    SAModelConstants sa_constants;
    EulerResidualNorms last_residual;
    double last_nut_residual = 0.0;
    std::vector<double> wall_distance; // Precomputed once at construction (WallDistance.h); one entry per cell

    // Combined inviscid-CFL + explicit-viscous-diffusion stable time step,
    // identical in structure to NavierStokesFVMSolver::compute_dt() except
    // that each face's kinematic viscosity term uses nu + nu_t (the LARGER
    // of the two adjoining cells' nu_t, mirroring how rho_min already picks
    // the more restrictive side), not just molecular nu -- nu_t can be
    // 10-1000x molecular nu in a real turbulent boundary layer, so this is
    // not a cosmetic change; see docs/archive/rans-spalart-allmaras-tracker.md
    // Phase 3's own risk callout.
    double compute_dt() const;

    // Builds the state on the far side of a boundary face for the INVISCID
    // flux only, per that face's patch's NSBoundaryType. Identical
    // methodology to NavierStokesFVMSolver::ghost_state() (duplicated, not
    // shared, per this tracker's architecture decision).
    EulerState ghost_state(const Face& face, const EulerState& U_L) const;

    // Builds the mesh.faces-sized boundary-value arrays GradientCalculator's
    // compute() needs for u, v, and T. Identical methodology to
    // NavierStokesFVMSolver::build_boundary_fields() (duplicated, not shared).
    void build_boundary_fields(const std::vector<double>& u, const std::vector<double>& v,
                                const std::vector<double>& T, std::vector<double>& boundary_u,
                                std::vector<double>& boundary_v, std::vector<double>& boundary_T) const;

    // Builds the mesh.faces-sized boundary-value array GradientCalculator's
    // compute() needs for nut: exactly 0 at a NoSlipWall (the SA model's own
    // wall condition, not case-configurable), the patch's farfield_nut at a
    // Farfield boundary, and a zero-order extrapolation (equal to the owning
    // cell's own current nut) at Outflow.
    void build_boundary_nut(const std::vector<double>& nut_field, std::vector<double>& boundary_nut) const;
};

#endif // RANSFVMSOLVER_H_INCLUDED
