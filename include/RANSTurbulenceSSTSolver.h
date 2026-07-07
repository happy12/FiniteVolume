// SPDX-License-Identifier: GPL-3.0-only
#ifndef RANSTURBULENCESSTSOLVER_H_INCLUDED
#define RANSTURBULENCESSTSOLVER_H_INCLUDED

#include <cassert>
#include <vector>

#include "EulerFVMSolver.h" // reuses EulerICMode/EulerInitialCondition/NumericalFluxScheme/EulerResidualNorms
#include "EulerState.h"
#include "GradientReconstruction.h"
#include "NavierStokesFVMSolver.h" // reuses NSBoundaryType/NSBoundaryCondition unchanged
#include "SSTKOmega.h"
#include "UnstructuredMesh.h"
#include "WallTraction.h"

// Per-patch boundary condition for the SST RANS solver, indexed in parallel
// with UnstructuredMesh::patches. Wraps NSBoundaryCondition unchanged, same
// rationale as RANSBoundaryConditionSA, plus the two extra prescribed values
// the k/omega transport equations need at a Farfield boundary. NoSlipWall's
// k/omega boundary values are NOT configurable here -- they follow the SST
// model's own wall definition (see build_boundary_k_omega() in
// RANSTurbulenceSSTSolver.cpp): k_wall = 0 exactly, omega_wall from the
// standard near-wall analytic formula. Deriving farfield_k/farfield_omega
// from a turbulence-intensity/eddy-viscosity-ratio case-file convention
// (matching SA's initial_nut precedent) is docs/sst-komega-tracker.md's
// Phase 3 scope, not this struct's -- these are raw prescribed values.
struct RANSBoundaryConditionSST {
    NSBoundaryCondition ns;
    double farfield_k = 0.0;
    double farfield_omega = 0.0;
};

// Explicit finite-volume solver for the 2D compressible RANS equations,
// closed with Menter's k-omega SST two-equation turbulence model (see
// SSTKOmega.h): NavierStokesFVMSolver's exact mean-flow equations (inviscid
// flux, viscous stress tensor, Fourier heat conduction), plus two extra
// transported scalars k and omega advected by the coupled mean flow
// velocity, with molecular viscosity/conductivity in the mean-flow viscous
// terms replaced by their turbulence-inclusive effective values:
//
//   mu_eff = mu + rho*nu_t     (nu_t = sst_eddy_viscosity(k, omega, S, Omega, F2, variant))
//   k_eff  = cp*(mu/Pr + rho*nu_t/Pr_t)
//
// Per the 2026-07-06 architecture decision confirmed in
// docs/sst-komega-tracker.md, this is deliberately a FOURTH solver class
// (after EulerFVMSolver/NavierStokesFVMSolver/RANSTurbulenceSASolver),
// duplicating NavierStokesFVMSolver's inviscid/viscous flux, ghost-state, and
// compute_dt() structure rather than adding a "turbulence_model" toggle to
// either existing RANS solver.
//
// Unlike RANSTurbulenceSASolver's nut (transported with NO rho-weighting, per
// SA's own classical form), k and omega are transported in Menter's
// compressible rho-weighted form (see SSTKOmega.h's own note on this genuine
// difference): the quantities actually conserved across a step are rho*k and
// rho*omega, exactly like EulerState conserves rho_u/rho_v (not u/v)
// directly. This class does NOT store rho*k/rho*omega as persistent members,
// though -- k/omega (PRIMITIVE, like nut) are the stored/exposed fields, and
// step() does the rho*k/rho*omega bookkeeping as local scratch values,
// dividing by the freshly-updated U[c].rho at the end of each step to recover
// primitive k/omega -- see step()'s own comment for why, and SSTKOmega.h's
// note on why k/omega need this rho-weighting that nut never did.
//
// k's transport equation (rho-weighted form, matching SSTKOmega.h's source terms exactly):
//
//   d(rho*k)/dt + div(rho*u*k) = k_production - k_destruction
//                              + div((mu + sigma_k*mu_t)*grad(k))
//
// omega's transport equation:
//
//   d(rho*omega)/dt + div(rho*u*omega) = omega_production - omega_destruction
//                                       + cross_diffusion
//                                       + div((mu + sigma_omega*mu_t)*grad(omega))
//
// sigma_k/sigma_omega are blended per-face via F1 (see SSTKOmega.h), mirroring
// how NavierStokesFVMSolver/RANSTurbulenceSASolver already face-interpolate
// their own diffusivities. k/omega's own advective flux upwinds rho*k/rho*omega
// directly using the same face-normal velocity Vn as the inviscid mean-flow
// flux, exactly mirroring RANSTurbulenceSASolver's nut_advective structure.
// k/omega's diffusive flux (unlike their advective flux) is forced to exactly
// zero at Farfield/Outflow boundaries, matching every other diffusive term's
// convention in this codebase.
class RANSTurbulenceSSTSolver {
public:
    // Input:
    //   input_mesh          - mesh topology/geometry to solve on (copied;
    //                          its compute_geometry() is called here)
    //   boundary_conditions - per-patch BC, indexed like mesh.patches
    //   gamma, gas_constant - as in NavierStokesFVMSolver
    //   mu, prandtl, prandtl_t - as in RANSTurbulenceSASolver
    //   cfl                 - CFL number, dimensionless, used to size each
    //                          step's dt (see compute_dt())
    //   initial_condition   - primitive-variable mean-flow initial condition
    //   initial_k, initial_omega - uniform initial/freestream k/omega applied
    //                          to every cell
    //   flux_scheme, gradient_scheme, exact_riemann_tol, exact_riemann_max_iter - as in RANSTurbulenceSASolver
    //   variant             - which SSTLimiterVariant this solver uses for
    //                          every eddy-viscosity/production-clip evaluation
    //                          (see SSTKOmega.h); fixed for the solver's
    //                          lifetime, not per-call
    //   sst_constants       - SST model constants (see SSTKOmega.h); defaults
    //                          to the NASA TMR set
    //   kato_launder        - if true, every step's source-term evaluation
    //                          uses the Kato-Launder production limiter (see
    //                          SSTKOmega.h's compute_sst_source_terms());
    //                          fixed for the solver's lifetime, default false
    // Output: constructs a solver with U/k/omega sized to mesh.cells.size()
    //         and the initial conditions above applied; also precomputes the
    //         general wall-distance field (WallDistance.h, feeding F1/F2 at
    //         every cell) once from the faces whose patch's NSBoundaryType is
    //         NoSlipWall -- same reasoning/precomputation-once rationale as
    //         RANSTurbulenceSASolver.
    RANSTurbulenceSSTSolver(const UnstructuredMesh& input_mesh, const std::vector<RANSBoundaryConditionSST>& boundary_conditions,
                   double gamma, double gas_constant, double mu, double prandtl, double prandtl_t, double cfl,
                   const EulerInitialCondition& initial_condition, double initial_k, double initial_omega,
                   NumericalFluxScheme flux_scheme, GradientScheme gradient_scheme, double exact_riemann_tol,
                   int exact_riemann_max_iter, SSTLimiterVariant variant,
                   const SSTModelConstants& sst_constants = SSTModelConstants{}, bool kato_launder = false);

    // Advances the solution by one explicit step, at a step size recomputed
    // every call from compute_dt() -- see compute_dt()'s own comment for why
    // this now also caps dt against the k/omega destruction terms' own
    // explicit stability limit, not just the diffusion limit RANSTurbulenceSASolver
    // and NavierStokesFVMSolver already cap against.
    void step();

    // Repeatedly calls step() to advance the solution.
    void run(int total_steps);

    const std::vector<EulerState>& field() const { return U; }
    const std::vector<double>& k_field() const { return k; }
    const std::vector<double>& omega_field() const { return omega; }
    double gas_gamma() const { return gamma; }
    const EulerResidualNorms& residual() const { return last_residual; }
    double k_residual() const { return last_k_residual; }
    double omega_residual() const { return last_omega_residual; }

    // Overwrites the current conserved-state AND k/omega fields, e.g. to
    // inject a checkpoint-resumed state after construction -- see
    // RANSTurbulenceSASolver::set_field() for the analogous precedent (k/omega
    // are PRIMITIVE here, same as that function's nut, despite step()'s
    // internal rho-weighted bookkeeping -- see this class's own comment).
    void set_field(const std::vector<EulerState>& new_U, const std::vector<double>& new_k,
                    const std::vector<double>& new_omega) {
        assert(new_U.size() == U.size());
        assert(new_k.size() == k.size());
        assert(new_omega.size() == omega.size());
        U = new_U;
        k = new_k;
        omega = new_omega;
        last_residual = EulerResidualNorms{};
        last_k_residual = 0.0;
        last_omega_residual = 0.0;
    }

    // Overwrites the CFL number compute_dt() uses on the next step() call --
    // e.g. for residual-based CFL ramping (see CflRamp.h/main.cpp's cfl_mode
    // handling). Takes effect starting with the next step(), not retroactively.
    void set_cfl(double new_cfl) { cfl = new_cfl; }

    // Wall traction diagnostics (see WallTraction.h) at every face in
    // 'wall_faces', from the CURRENT flow field -- identical on-demand,
    // nothing-retained-from-step() contract as
    // RANSTurbulenceSASolver::compute_wall_traction_samples(), except
    // effective_viscosity uses this solver's own nu_t (sst_eddy_viscosity()).
    std::vector<WallFaceSample> compute_wall_traction_samples(const std::vector<int>& wall_faces) const;

    // Boundary-layer thickness diagnostics -- identical contract to
    // RANSTurbulenceSASolver::compute_boundary_layer_profile_samples().
    std::vector<BoundaryLayerProfile> compute_boundary_layer_profile_samples(const std::vector<int>& wall_faces,
                                                                              double u_edge,
                                                                              int max_cells_per_march) const;

    // Point-location alternative -- identical contract to
    // RANSTurbulenceSASolver::compute_boundary_layer_profile_samples_point_location().
    std::vector<BoundaryLayerProfile> compute_boundary_layer_profile_samples_point_location(
        const std::vector<int>& wall_faces, double u_edge, double max_distance, int n_samples) const;

private:
    UnstructuredMesh mesh;
    std::vector<EulerState> U;
    std::vector<double> k;     // primitive (not rho-weighted); see class comment
    std::vector<double> omega; // primitive (not rho-weighted); see class comment
    std::vector<RANSBoundaryConditionSST> bcs;
    double gamma, gas_constant, mu, prandtl, prandtl_t, cfl;
    NumericalFluxScheme flux_scheme;
    GradientCalculator gradient_calc; // Reconstructs cell gradients of u, v, T, k, and omega every step
    double exact_riemann_tol;
    int exact_riemann_max_iter;
    SSTLimiterVariant variant;
    SSTModelConstants sst_constants;
    bool kato_launder;
    EulerResidualNorms last_residual;
    double last_k_residual = 0.0;
    double last_omega_residual = 0.0;
    std::vector<double> wall_distance; // Precomputed once at construction (WallDistance.h); feeds F1/F2 at every cell

    // Combined inviscid-CFL + explicit-viscous-diffusion + explicit-source
    // stable time step. The diffusion part is identical in structure to
    // RANSTurbulenceSASolver::compute_dt() (same nu_lam + nu_t_max kinematic
    // viscosity term, same face_normal_distance() length scale) -- and stays
    // sufficient for k/omega's OWN diffusion terms too without a separate
    // sigma_k/sigma_omega-scaled term, since sigma_k/sigma_omega are both
    // <= 1.0 everywhere (see SSTKOmega.h's constant table), so k/omega's
    // effective diffusivity (nu_lam + sigma*nu_t) never exceeds the
    // mean-flow's own (nu_lam + nu_t) term this formula already uses.
    //
    // The source part is new, and is NOT a diffusion-stability estimate:
    // omega's near-wall wall boundary value (60*nu/(beta1*d1^2), see
    // build_boundary_k_omega()) can be very large on a well-resolved
    // near-wall mesh, and once the near-wall cell's own omega rises toward
    // that value (via diffusive coupling to the wall), its destruction terms
    // (beta_star*omega*k, beta*omega^2 -- both per-unit-mass rates once
    // divided through by rho, see step()'s bookkeeping) can become large
    // enough that a plain diffusion-CFL dt overshoots them in a single
    // explicit step -- the same class of stiffness RANSTurbulenceSASolver's own
    // destruction term (r/g/fw chain) can hit, but reachable here at a much
    // larger magnitude since omega itself, not just a clipped ratio, is what
    // grows large. Capped per cell at cfl / max(beta_star*omega, beta*omega)
    // (beta blended via that cell's own F1), same 'cfl' parameter reused for
    // this limit as for the inviscid/diffusive ones above, rather than a
    // second independently-tunable safety factor.
    double compute_dt() const;

    // Builds the state on the far side of a boundary face for the INVISCID
    // flux only. Identical methodology to RANSTurbulenceSASolver::ghost_state()
    // (duplicated, not shared, per this project's architecture precedent).
    EulerState ghost_state(const Face& face, const EulerState& U_L) const;

    // Builds the mesh.faces-sized boundary-value arrays GradientCalculator's
    // compute() needs for u, v, and T. Identical methodology to
    // RANSTurbulenceSASolver::build_boundary_fields() (duplicated, not shared).
    void build_boundary_fields(const std::vector<double>& u, const std::vector<double>& v,
                                const std::vector<double>& T, std::vector<double>& boundary_u,
                                std::vector<double>& boundary_v, std::vector<double>& boundary_T) const;

    // Builds the mesh.faces-sized boundary-value arrays needed for k/omega:
    // PRIMITIVE boundary_k/boundary_omega (what GradientCalculator's compute()
    // needs) and CONSERVED boundary_rho_k/boundary_rho_omega (what the
    // advective flux's upwind "R" value needs, see class comment). At a
    // NoSlipWall: boundary_k = 0 exactly (the SST model's own wall condition,
    // not case-configurable); boundary_omega = 60*nu/(beta1*d1^2), d1 =
    // face_normal_distance(mesh, face) -- see docs/sst-komega-tracker.md's
    // Phase 2 section for why this reuses that existing generic helper
    // rather than the general per-cell wall_distance field (which reports
    // distance to the nearest wall from an arbitrary cell centroid, not
    // specifically this face's own first-cell-centroid distance); both
    // boundary_rho_k/boundary_rho_omega then use the OWNING cell's own
    // (interior) density, there being no separately-tracked wall density.
    // At a Farfield boundary: boundary_k/omega = the patch's prescribed
    // farfield_k/farfield_omega, and boundary_rho_k/rho_omega use the
    // prescribed farfield_state's own density. At Outflow: zero-order
    // extrapolation of k/omega AND of the already-conserved rho*k/rho*omega
    // (i.e. boundary_rho_k[i] = U[cl].rho * k[cl] exactly, not recomputed
    // from a separately-extrapolated density).
    void build_boundary_k_omega(const std::vector<double>& k_field, const std::vector<double>& omega_field,
                                  std::vector<double>& boundary_k, std::vector<double>& boundary_omega,
                                  std::vector<double>& boundary_rho_k, std::vector<double>& boundary_rho_omega) const;
};

#endif // RANSTURBULENCESSTSOLVER_H_INCLUDED
