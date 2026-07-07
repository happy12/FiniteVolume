// SPDX-License-Identifier: GPL-3.0-only
#ifndef NAVIERSTOKESFVMSOLVER_H_INCLUDED
#define NAVIERSTOKESFVMSOLVER_H_INCLUDED

#include <cassert>
#include <vector>

#include "EulerFVMSolver.h" // reuses EulerICMode/EulerInitialCondition/NumericalFluxScheme
#include "EulerState.h"
#include "GradientReconstruction.h"
#include "UnstructuredMesh.h"
#include "WallTraction.h"

// Boundary condition kind for the Navier-Stokes solver. Separate from
// EulerBoundaryType (EulerFVMSolver.h) because a viscous wall is genuinely
// different from Euler's slip wall: it pins BOTH velocity components to
// zero (not just the normal component) and needs an associated thermal
// condition (isothermal or adiabatic) that an inviscid wall has no use for.
enum class NSBoundaryType : int {
    NoSlipWall = 0, // u = wall_u, v = wall_v at the wall; thermal condition per is_isothermal_wall
    Farfield = 1,   // fixed prescribed state (inviscid contribution only -- see class comment)
    Outflow = 2     // zero-gradient extrapolation (inviscid contribution only)
};

// Per-patch boundary condition for the Navier-Stokes solver, indexed in
// parallel with UnstructuredMesh::patches (i.e. by Face::patch_id).
struct NSBoundaryCondition {
    NSBoundaryType type = NSBoundaryType::NoSlipWall;
    double wall_u = 0.0, wall_v = 0.0; // only meaningful when type == NoSlipWall; the wall's own velocity
                                        // (0, 0 = stationary; nonzero = a moving wall, e.g. Couette flow)
    bool is_isothermal_wall = false; // only meaningful when type == NoSlipWall; false = adiabatic (zero heat flux)
    double wall_temperature = 0.0;   // only meaningful when is_isothermal_wall
    EulerState farfield_state;       // only meaningful when type == Farfield; units as in EulerState.h
};

// Domain-wide summary of how well the current mesh resolves the local
// strain field's dissipation length scale (see
// NavierStokesFVMSolver::compute_resolution_diagnostics()). Fields are
// aggregated only over cells with non-negligible local dissipation;
// n_active reports how many cells that was, so a caller can tell "no shear
// anywhere yet" (n_active == 0) apart from "well resolved" (n_active > 0,
// ratios near 1).
struct ResolutionDiagnostics {
    double min_ratio = 0.0;  // min over active cells of h/eta (h = sqrt(cell volume), eta = local Kolmogorov length)
    double max_ratio = 0.0;
    double mean_ratio = 0.0;
    size_t n_active = 0;     // number of cells with non-negligible local dissipation rate
};

// Explicit finite-volume solver for the 2D compressible Navier-Stokes
// equations on an arbitrary unstructured polygon mesh: the compressible
// Euler equations (see EulerFVMSolver.h) plus a Newtonian viscous stress
// tensor and Fourier heat conduction.
//
//   dU/dt + div(F_inviscid(U)) - div(F_viscous(u, v, T, grad_u, grad_v, grad_T)) = 0
//
// Methodology: identical inviscid flux machinery to EulerFVMSolver (same
// NumericalFluxScheme choices, same two-OpenMP-pass structure -- this is
// deliberately a second class rather than a toggle added to EulerFVMSolver,
// per the 2026-07-05 architecture decision recorded in
// docs/navier-stokes-tracker.md), with an added viscous flux built from the
// non-orthogonality-corrected face gradients (see GradientReconstruction.h)
// of the primitive velocity components and temperature:
//
//   tau_xx = mu*(2*du/dx - (2/3)*(du/dx + dv/dy))
//   tau_yy = mu*(2*dv/dy - (2/3)*(du/dx + dv/dy))
//   tau_xy = mu*(du/dy + dv/dx)
//   q_x = -k*dT/dx, q_y = -k*dT/dy   (k = mu*cp/Pr, cp = gamma*R/(gamma-1))
//
//   F_viscous . n = [0,
//                    tau_xx*nx + tau_xy*ny,
//                    tau_xy*nx + tau_yy*ny,
//                    (u*tau_xx + v*tau_xy - q_x)*nx + (u*tau_xy + v*tau_yy - q_y)*ny]
//
// Since face_gradient() only corrects the NORMAL component of a face's
// gradient (the component that would otherwise decouple/checkerboard -- see
// docs/navier-stokes-tracker.md Phase 1), the full 2D gradient vector needed
// for the stress tensor's off-diagonal/tangential terms is assembled by
// taking the interpolated average gradient and replacing just its normal
// component with the corrected one (see corrected_face_gradient_vector() in
// GradientReconstruction.h).
//
// Temperature is recovered via the ideal-gas relation p = rho*R*T (R = the
// specific gas constant, a case-file parameter this solver introduces --
// EulerState.h's temperature() needs it as an explicit argument since Euler
// itself never needed R).
//
// Viscous flux is set to exactly zero at Farfield/Outflow boundaries -- the
// intended use is placing those boundaries far enough from any wall that
// viscous gradients there are negligible; this solver does not attempt to
// extrapolate a meaningful viscous flux through an open boundary.
class NavierStokesFVMSolver {
public:
    // Input:
    //   input_mesh          - mesh topology/geometry to solve on (copied;
    //                          its compute_geometry() is called here)
    //   boundary_conditions - per-patch BC, indexed like mesh.patches
    //   gamma               - ratio of specific heats, dimensionless (1.4 for air)
    //   gas_constant        - specific gas constant R, in mesh-consistent
    //                          units (p = rho*R*T); only used for temperature/
    //                          heat conduction, never for pressure/sound speed
    //   mu                  - dynamic viscosity, mesh-consistent units; 0 =
    //                          no viscous stress (but heat conduction, if
    //                          mu > 0 via k = mu*cp/Pr, also then vanishes,
    //                          since k scales with mu)
    //   prandtl             - Prandtl number, dimensionless (0.72 for air)
    //   cfl                 - CFL number, dimensionless, used to size each
    //                          step's dt (see compute_dt())
    //   initial_condition   - primitive-variable initial condition (see
    //                          EulerInitialCondition, EulerFVMSolver.h)
    //   flux_scheme         - which numerical flux the inviscid term uses at
    //                          every face (see NumericalFluxScheme, EulerFVMSolver.h)
    //   gradient_scheme     - which GradientCalculator scheme reconstructs
    //                          cell gradients of u/v/T every step (see GradientReconstruction.h)
    //   exact_riemann_tol, exact_riemann_max_iter - ignored unless flux_scheme == Exact
    // Output: constructs a solver with U sized to mesh.cells.size() and the
    //         initial condition above applied.
    NavierStokesFVMSolver(const UnstructuredMesh& input_mesh, const std::vector<NSBoundaryCondition>& boundary_conditions,
                           double gamma, double gas_constant, double mu, double prandtl, double cfl,
                           const EulerInitialCondition& initial_condition, NumericalFluxScheme flux_scheme,
                           GradientScheme gradient_scheme, double exact_riemann_tol, int exact_riemann_max_iter);

    // Advances the solution by one explicit step, at a step size recomputed
    // every call from compute_dt() (both the inviscid CFL limit and the
    // explicit-viscous-diffusion stability limit).
    void step();

    // Repeatedly calls step() to advance the solution.
    void run(int total_steps);

    const std::vector<EulerState>& field() const { return U; }
    double gas_gamma() const { return gamma; }
    const EulerResidualNorms& residual() const { return last_residual; }

    // Overwrites the current conserved-state field, e.g. to inject a
    // checkpoint-resumed state after construction.
    void set_field(const std::vector<EulerState>& new_U) {
        assert(new_U.size() == U.size());
        U = new_U;
        last_residual = EulerResidualNorms{};
    }

    // Overwrites the CFL number compute_dt() uses on the next step() call --
    // e.g. for residual-based CFL ramping (see CflRamp.h/main.cpp's cfl_mode
    // handling). Takes effect starting with the next step(), not retroactively.
    void set_cfl(double new_cfl) { cfl = new_cfl; }

    // Estimates, from the CURRENT flow field, how well the mesh resolves the
    // smallest (Kolmogorov) length scale of the local strain field --
    // meant as a resolution-adequacy diagnostic for treating a run as a
    // fully-resolved unsteady 2D simulation, NOT a claim about 3D
    // turbulence: this solver is 2D-only, and classical Kolmogorov scaling
    // (and the energy cascade it comes from) is a 3D result -- 2D
    // turbulence has different phenomenology (no vortex stretching, an
    // inverse energy cascade). Treat this as a pragmatic resolution
    // heuristic, not a validated 2D theory. See
    // docs/navier-stokes-tracker.md Phase 5 for the full discussion.
    //
    // Methodology: per cell, from the reconstructed velocity gradients,
    //   S_11 = du/dx, S_22 = dv/dy, S_12 = 0.5*(du/dy + dv/dx)
    //   epsilon = 2*nu*(S_11^2 + S_22^2 + 2*S_12^2)   (local dissipation rate)
    //   eta = (nu^3 / epsilon)^0.25                    (Kolmogorov length)
    //   h = sqrt(cell volume)                          (characteristic cell length)
    // h/eta noticeably greater than 1 indicates the mesh is coarser than
    // the locally estimated dissipation length scale there. Cells whose
    // local dissipation rate is below a small floor (no meaningful local
    // shear -- e.g. undisturbed far-field, or before any flow has
    // developed) are excluded from the aggregate: eta -> infinity there is
    // an absence of anything to resolve, not a resolution problem.
    //
    // Input:  none (uses the current field and mesh; recomputes u/v
    //         gradients fresh via the same GradientCalculator used by
    //         step(), since they aren't retained between steps)
    // Returns: domain-wide min/max/mean of h/eta over cells with
    //          non-negligible local dissipation, and how many cells that was
    ResolutionDiagnostics compute_resolution_diagnostics() const;

    // Wall traction diagnostics (see WallTraction.h) at every face in
    // 'wall_faces', from the CURRENT flow field -- same "on-demand
    // post-processing pass, nothing retained from step()" contract as
    // compute_resolution_diagnostics(): u/v/T, boundary fields, and
    // gradients are all recomputed fresh here. effective_viscosity is
    // uniform mu, since this solver has no turbulence closure (see
    // RANSTurbulenceSASolver for the mu + rho*nu_t analogue).
    //
    // Input:  wall_faces - indices into mesh.faces to sample; every entry
    //         should be a boundary face whose patch's NSBoundaryType is
    //         NoSlipWall (not checked here -- see WallTraction.h)
    // Returns: one WallFaceSample per entry of 'wall_faces', in the same order
    std::vector<WallFaceSample> compute_wall_traction_samples(const std::vector<int>& wall_faces) const;

    // Boundary-layer thickness diagnostics (see
    // compute_boundary_layer_profiles(), WallTraction.h) at every face in
    // 'wall_faces', from the CURRENT flow field -- same on-demand,
    // nothing-retained-from-step() contract as
    // compute_wall_traction_samples()/compute_resolution_diagnostics().
    //
    // Input:  wall_faces          - as in compute_wall_traction_samples()
    //         u_edge              - boundary-layer edge velocity magnitude
    //                               (see compute_boundary_layer_profiles())
    //         max_cells_per_march - marching cap (see compute_boundary_layer_profiles())
    // Returns: one BoundaryLayerProfile per entry of 'wall_faces', in the same order
    std::vector<BoundaryLayerProfile> compute_boundary_layer_profile_samples(const std::vector<int>& wall_faces,
                                                                              double u_edge,
                                                                              int max_cells_per_march) const;

    // Point-location alternative to compute_boundary_layer_profile_samples()
    // above (see compute_boundary_layer_profiles_point_location(),
    // WallTraction.h) -- same on-demand, nothing-retained-from-step()
    // contract, useful when the mesh doesn't have a clean wall-normal cell
    // stacking near the wall (see that function's own class comment).
    //
    // Input:  wall_faces    - as in compute_wall_traction_samples()
    //         u_edge        - as in compute_boundary_layer_profile_samples()
    //         max_distance, n_samples - see compute_boundary_layer_profiles_point_location()
    // Returns: one BoundaryLayerProfile per entry of 'wall_faces', in the same order
    std::vector<BoundaryLayerProfile> compute_boundary_layer_profile_samples_point_location(
        const std::vector<int>& wall_faces, double u_edge, double max_distance, int n_samples) const;

private:
    UnstructuredMesh mesh;
    std::vector<EulerState> U;
    std::vector<NSBoundaryCondition> bcs;
    double gamma, gas_constant, mu, prandtl, cfl;
    NumericalFluxScheme flux_scheme;
    GradientCalculator gradient_calc; // Reconstructs cell gradients of u, v, and T every step
    double exact_riemann_tol;
    int exact_riemann_max_iter;
    EulerResidualNorms last_residual;

    // Combined inviscid-CFL + explicit-viscous-diffusion stable time step
    // (Blazek-style): each face's candidate dt uses length / (Smax +
    // 2*nu/length), nu = mu/rho being the kinematic viscosity, since
    // explicit diffusion's own stability limit is dt < length^2/(2*nu); the
    // true dt is the minimum candidate over the whole mesh.
    double compute_dt() const;

    // Builds the state on the far side of a boundary face (the "ghost"
    // state used in place of a real neighbor cell) for the INVISCID flux
    // only, per that face's patch's NSBoundaryType.
    // Methodology: Farfield/Outflow are identical to EulerFVMSolver's; for
    // NoSlipWall, BOTH velocity components are mirrored about the wall's own
    // velocity (u_ghost = 2*wall_u - u_L, not just the normal component
    // mirrored about zero, unlike Euler's slip wall), so the inviscid
    // flux's own central-average term comes out to exactly (wall_u, wall_v)
    // at the wall, consistent with the no-slip condition the viscous term
    // also enforces.
    EulerState ghost_state(const Face& face, const EulerState& U_L) const;

    // Builds the mesh.faces-sized boundary-value arrays GradientCalculator's
    // compute() needs for u, v, and T. A no-slip wall pins u = wall_u,
    // v = wall_v (Dirichlet) and either pins T = wall_temperature (isothermal) or
    // extrapolates T from the owning cell (adiabatic -- the heat flux
    // itself, not this gradient stencil value, is what step() forces to
    // zero for an adiabatic wall). A Farfield/Outflow face is always
    // extrapolated (zero-order), since no viscous flux is used there.
    void build_boundary_fields(const std::vector<double>& u, const std::vector<double>& v,
                                const std::vector<double>& T, std::vector<double>& boundary_u,
                                std::vector<double>& boundary_v, std::vector<double>& boundary_T) const;
};

#endif // NAVIERSTOKESFVMSOLVER_H_INCLUDED
