// SPDX-License-Identifier: GPL-3.0-only
#ifndef WALLTRACTION_H_INCLUDED
#define WALLTRACTION_H_INCLUDED

#include <vector>

#include "GradientReconstruction.h"
#include "UnstructuredMesh.h"

// Selects which of compute_boundary_layer_profiles()/
// compute_boundary_layer_profiles_point_location() a caller (e.g. a
// wall_profile_file case-file key) uses to derive boundary-layer thickness
// estimates.
enum class BoundaryLayerMethod {
    Marching,      // cheap, assumes a wall-normal cell stacking near the wall (see compute_boundary_layer_profiles())
    PointLocation  // more expensive, no near-wall-topology assumption (see compute_boundary_layer_profiles_point_location())
};

// Wall diagnostics (skin friction, pressure/force/moment integration,
// boundary-layer thickness) as a separate, on-demand post-processing pass
// over only a mesh's NoSlipWall boundary faces -- see
// docs/wall-diagnostics-plan.md for the full design discussion. Kept
// solver-agnostic (no dependency on NavierStokesFVMSolver/RANSTurbulenceSASolver),
// same reasoning as WallDistance.h: it depends only on mesh geometry, a wall
// face list, and flow-field snapshots handed in by whichever solver owns
// them, not on which solver produced those fields. 'effective_viscosity' is
// per-cell so a caller can pass mu (NavierStokesFVMSolver) or mu + rho*nu_t
// (RANSTurbulenceSASolver) without this module needing to know which.

// One wall-tangential traction sample at a single NoSlipWall boundary face.
// Deliberately reference-agnostic (no rho_ref/V_ref/p_ref baked in) --
// Cf/Cp/y+ are derived downstream from these raw fields (see
// skin_friction_coefficient()/pressure_coefficient()/wall_y_plus() below),
// keeping this struct reusable regardless of which reference quantities a
// caller ultimately picks.
struct WallFaceSample {
    int face_index = -1;
    double x_mid = 0.0, y_mid = 0.0;      // face midpoint
    double p = 0.0;                        // wall pressure (see compute_wall_traction()'s "wall pressure" note)
    double rho = 0.0;                      // wall-adjacent cell's density (needed for y+)
    double effective_viscosity = 0.0;      // mu (NS) or mu + rho*nu_t (RANS) in the wall-adjacent cell
    double tau_wall = 0.0;                 // wall-tangential shear stress magnitude, signed: positive = same sense as (tx, ty)
    double tx = 0.0, ty = 0.0;             // unit tangent direction used for tau_wall's sign (face normal rotated +90 deg)
    double y_wall_normal = 0.0;            // first-cell wall-normal distance (face_normal_distance() to the owning cell's centroid)
};

// Computes a WallFaceSample for every face index in 'wall_faces' (caller-
// supplied, e.g. by filtering mesh.faces for NoSlipWall patches -- see
// WallDistance.h's identical "caller supplies the face list" precedent).
//
// Methodology: reuses face_gradient() and corrected_face_gradient_vector()
// (GradientReconstruction.h) to get the full 2D gradient vector of u and v at
// each wall face exactly as NavierStokesFVMSolver::step()/RANSTurbulenceSASolver::step()
// do for their own viscous flux, then assembles the same Newtonian stress
// tensor (tau_xx, tau_yy, tau_xy) and projects it onto the face's outward
// normal to get the traction vector, and onto the tangent direction
// (face normal rotated +90 deg) for the signed scalar tau_wall.
//
// Wall pressure: neither solver retains a separate boundary pressure state --
// a NoSlipWall's ghost state (see NavierStokesFVMSolver::ghost_state()) exists
// only for the inviscid flux, not as a reconstructed boundary value. Wall
// pressure here is therefore the owning cell's own pressure, p[face.cell_left]
// -- a first-order approximation, not a boundary-extrapolated one, consistent
// with this project's existing first-order-approximation conventions
// elsewhere (e.g. AdvectionDiffusionFVMSolver's Neumann-boundary gradient
// stencil value).
//
// Input:
//   mesh                 - mesh 'wall_faces' indexes into; must already have
//                          compute_geometry() called
//   wall_faces           - indices into mesh.faces to sample (every entry
//                          should be a boundary face, i.e. cell_right == -1)
//   u, v, p, rho          - cell-centered fields, one value per mesh.cells entry
//   effective_viscosity   - cell-centered field: mu (uniform, NS) or
//                          mu + rho*nu_t (RANS), one value per mesh.cells entry
//   boundary_u, boundary_v - mesh.faces-sized prescribed boundary values (the
//                          same arrays a caller's own
//                          build_boundary_fields()-style helper already builds
//                          for its step(); only entries at 'wall_faces'
//                          indices are read -- for a NoSlipWall these are
//                          exactly (wall_u, wall_v), the wall's own velocity)
//   grad_u, grad_v         - cell-centered reconstructed gradients (see
//                          GradientCalculator), one value per mesh.cells entry
// Returns: one WallFaceSample per entry of 'wall_faces', in the same order
std::vector<WallFaceSample> compute_wall_traction(
    const UnstructuredMesh& mesh, const std::vector<int>& wall_faces,
    const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& p, const std::vector<double>& rho,
    const std::vector<double>& effective_viscosity,
    const std::vector<double>& boundary_u, const std::vector<double>& boundary_v,
    const std::vector<Gradient2>& grad_u, const std::vector<Gradient2>& grad_v);

// Freestream/reference quantities used to nondimensionalize WallFaceSample
// into standard coefficients (Cf, Cp, Cd, Cl, Cm) and to define the point
// moments are taken about. Kept separate from WallFaceSample itself so
// compute_wall_traction() stays reference-agnostic (see its class comment) --
// a caller only needs these once it wants a *coefficient*, never for the raw
// tau_wall/y+ diagnostics.
struct WallReferenceQuantities {
    double rho_ref = 1.0;                                   // freestream density, for dynamic pressure
    double velocity_x_ref = 1.0, velocity_y_ref = 0.0;       // freestream velocity; drag direction = this vector normalized
    double p_ref = 0.0;                                      // freestream static pressure, for Cp
    double length_ref = 1.0;                                 // length scale for Cd/Cl/Cm (e.g. chord); 1.0 = per-unit-span coefficient
    double moment_reference_x = 0.0, moment_reference_y = 0.0; // point Cm is taken about
};

// Cf = tau_wall / (0.5 * rho_ref * V_ref^2)
double skin_friction_coefficient(const WallFaceSample& sample, const WallReferenceQuantities& ref);

// Cp = (p - p_ref) / (0.5 * rho_ref * V_ref^2)
double pressure_coefficient(const WallFaceSample& sample, const WallReferenceQuantities& ref);

// y+ = rho * u_tau * y_wall_normal / effective_viscosity, u_tau = sqrt(|tau_wall| / rho).
// Needs no reference quantities at all -- every term is already on 'sample'.
double wall_y_plus(const WallFaceSample& sample);

// Per-patch (and, at patch_id == -1, domain-total-across-every-sample) force
// and moment integration, nondimensionalized against 'ref'. An airfoil patch
// and a tunnel-wall patch are never silently summed together unless the
// caller passes samples from both under the same call and also wants the
// patch_id == -1 total -- see compute_wall_forces()'s own contract.
struct WallForceReport {
    int patch_id = -1;                                    // -1 = domain total across every sample passed in
    double friction_drag = 0.0, pressure_drag = 0.0, total_drag = 0.0;
    double cd_friction = 0.0, cd_pressure = 0.0, cd_total = 0.0;
    double lift = 0.0, cl = 0.0;
    double moment = 0.0, cm = 0.0;
};

// Sums, over 'samples' grouped by their face's patch_id (via
// mesh.faces[sample.face_index].patch_id):
//   friction_force = sum( face.area * tau_wall * (tx, ty) )
//   pressure_force = sum( face.area * p * (nx, ny) )   -- outward normal; force ON the wall from the fluid
//   moment         = sum( r x (friction_force_i + pressure_force_i) ), r from ref.moment_reference_{x,y}
// then projects the force onto ref's freestream direction for drag and its
// perpendicular for lift.
//
// Input:  mesh, samples - as computed by compute_wall_traction()
//         ref            - see WallReferenceQuantities
// Returns: one WallForceReport per distinct patch_id present in 'samples',
//          plus one final entry with patch_id == -1 summing every sample
//          (the domain total). Empty if 'samples' is empty.
std::vector<WallForceReport> compute_wall_forces(const UnstructuredMesh& mesh,
                                                   const std::vector<WallFaceSample>& samples,
                                                   const WallReferenceQuantities& ref);

// One wall-normal velocity profile starting from a single NoSlipWall
// boundary face, and the boundary-layer thickness estimates derived from it.
// Produced by either compute_boundary_layer_profiles() (cell-marching) or
// compute_boundary_layer_profiles_point_location() (point-location) -- see
// each function's own class comment for its sampling methodology.
struct BoundaryLayerProfile {
    int wall_face_index = -1;
    double x_mid = 0.0, y_mid = 0.0;
    double delta_99 = 0.0;
    double displacement_thickness = 0.0;
    double momentum_thickness = 0.0;
    double shape_factor = 0.0;             // H = displacement_thickness / momentum_thickness
    int n_cells_marched = 0;               // diagnostic: number of samples taken beyond the wall itself (cells
                                            // marched through, or point-location samples along the profile) before
                                            // stopping -- did it hit its cap/plateau at 0.99*u_edge first?
};

// For every face in 'wall_faces', marches cell-to-cell away from the wall
// along the locally most-wall-normal-aligned neighbor, sampling the
// tangential velocity component at each visited cell centroid, then derives
// delta_99/displacement_thickness/momentum_thickness/shape_factor from the
// resulting (distance, u_tangential) profile.
//
// Methodology (see docs/wall-diagnostics-plan.md's "Boundary-layer thickness"
// section for the full discussion of why this approach was chosen over a
// point-location/ray-casting alternative):
//   1. Start at the wall face's owning cell. The wall face's own INWARD
//      normal (into the fluid, i.e. the negation of Face::nx/ny, which points
//      outward from the fluid cell into the wall) is the fixed marching
//      direction, and the face's tangent (rotated +90 deg) is the fixed
//      sampling direction for the whole march -- neither is re-derived per cell.
//   2. From the current cell, among its OTHER faces (excluding the one just
//      crossed), cross into whichever neighbor's shared face has an outward
//      normal most nearly parallel to the marching direction (largest dot
//      product, and strictly positive) -- the same "most-aligned" idea
//      face_normal_distance() uses for a single step, repeated cell-to-cell.
//      A boundary face (no neighbor) is never chosen, and a non-positive
//      best alignment stops the march rather than accepting a lateral/
//      backward step (this matters at the top of a column: once the
//      directly-ahead face is a boundary, the only remaining interior faces
//      can be purely lateral, and picking the "least bad" of those would
//      silently walk sideways instead of correctly stopping). Stop marching
//      if no interior face remains, no candidate has positive alignment, or
//      once max_cells_per_march is reached.
//   3. Accumulate wall-normal distance via face_normal_distance() at each
//      step; sample u_tangential = |u . tangent| (fabs -- a "top"-type wall's
//      outward-normal-derived tangent points opposite a "bottom"-type wall's
//      even when the physical flow moves the same direction past both, and
//      the ratios below assume u/u_edge approaches +1, not -1) at each
//      visited cell centroid.
//   4. Stop early once a sample's u_tangential first reaches 0.99*u_edge
//      (nothing past the boundary layer is needed for delta_99/delta*/theta).
//   5. delta_99: linear interpolation of the profile where u/u_edge first
//      crosses 0.99. displacement_thickness/momentum_thickness: trapezoidal
//      integration of (1 - u/u_edge) and (u/u_edge)*(1 - u/u_edge) from
//      distance 0 (using the wall's own prescribed tangential velocity,
//      boundary_u/boundary_v, as the distance-0 sample) out to the last
//      sampled point. shape_factor = displacement_thickness / momentum_thickness.
//
// Explicit, disclosed limitation, DEMONSTRATED not hypothetical: this
// assumes a locally wall-normal-ish structured cell stacking near the wall
// (true for boundary-layer-clustered meshes, reproducing
// --verify-flat-plate's hardcoded column-walk result to under 1%) --
// `--verify-bl-marching-unstructured` shows a systematic ~20-25% delta_99
// error on a genuinely unstructured (checkerboard-triangulated) near-wall
// mesh, root-caused to face_normal_distance() projecting onto each crossed
// face's own normal rather than the fixed march direction once a step
// crosses a diagonal face (see docs/wall-diagnostics-plan.md's "Phase 2
// stress-test finding"). compute_boundary_layer_profiles_point_location()
// below is the point-location alternative built specifically to not share
// this failure mode -- prefer it on a mesh without a clean wall-normal cell
// stacking, at the cost of its own O(cells) per-sample-point search.
//
// Input:
//   mesh                 - mesh 'wall_faces' indexes into; must already have
//                          compute_geometry() called
//   wall_faces           - indices into mesh.faces to march from (every entry
//                          should be a boundary face, i.e. cell_right == -1)
//   u, v                  - cell-centered velocity fields, one per mesh.cells entry
//   boundary_u, boundary_v - mesh.faces-sized prescribed boundary values (see
//                          compute_wall_traction()'s identical parameter);
//                          only entries at 'wall_faces' indices are read
//   u_edge                - boundary-layer edge velocity magnitude (e.g. the
//                          freestream speed) that delta_99/delta*/theta are
//                          normalized against
//   max_cells_per_march   - marching cap; a march that reaches this without
//                          its sample plateauing at 0.99*u_edge stops anyway
//                          (n_cells_marched reports this on the result)
// Returns: one BoundaryLayerProfile per entry of 'wall_faces', in the same order
std::vector<BoundaryLayerProfile> compute_boundary_layer_profiles(
    const UnstructuredMesh& mesh, const std::vector<int>& wall_faces,
    const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& boundary_u, const std::vector<double>& boundary_v,
    double u_edge, int max_cells_per_march);

// Point-location alternative to compute_boundary_layer_profiles() above,
// built specifically to not share its structured-cell-stacking assumption
// (see that function's "Explicit, disclosed limitation" note and
// docs/wall-diagnostics-plan.md's Phase 2 stress-test finding, which is what
// motivated this). Produces the same BoundaryLayerProfile via the identical
// downstream delta_99/thickness-integral math (see derive_boundary_layer_profile()
// in WallTraction.cpp) -- the two functions differ only in how the raw
// (distance, u_tangential) samples are obtained.
//
// Methodology: for every face in 'wall_faces', walks straight out from the
// face midpoint along its fixed inward normal at 'n_samples' evenly-spaced
// offsets up to 'max_distance', literally locating (brute-force,
// point-in-polygon over every cell -- no spatial acceleration structure,
// same "simple first, document scaling limits" precedent as WallDistance.h)
// which mesh cell contains each sample point and reading that cell's
// (piecewise-constant, per this project's FVM convention) velocity there.
// This is exact regardless of the mesh's near-wall topology -- unlike
// cell-marching, there is no neighbor-alignment choice to get wrong -- at
// the cost of an O(cells) search per sample point (versus marching's O(1)-ish
// per step), so this is meaningfully more expensive on a large mesh; prefer
// compute_boundary_layer_profiles() when the mesh does have a clean
// wall-normal cell stacking (the common case for a boundary-layer-clustered
// mesh) and reach for this one when it doesn't.
//
// Stops early (same 0.99*u_edge plateau rule as the marching version) or
// once a sample point falls outside the mesh domain (locate_cell() finds no
// containing cell) or 'n_samples' is reached, whichever comes first.
//
// Input:
//   mesh                 - mesh 'wall_faces' indexes into; must already have
//                          compute_geometry() called
//   wall_faces           - indices into mesh.faces to sample from (every
//                          entry should be a boundary face, i.e. cell_right == -1)
//   u, v                  - cell-centered velocity fields, one per mesh.cells entry
//   boundary_u, boundary_v - mesh.faces-sized prescribed boundary values (see
//                          compute_wall_traction()'s identical parameter);
//                          only entries at 'wall_faces' indices are read
//   u_edge                - boundary-layer edge velocity magnitude (e.g. the
//                          freestream speed) that delta_99/delta*/theta are
//                          normalized against
//   max_distance          - furthest wall-normal distance to sample; should
//                          comfortably exceed the expected boundary-layer
//                          thickness
//   n_samples             - number of evenly-spaced sample points between 0
//                          and max_distance (n_cells_marched on the result
//                          reports how many were actually taken before
//                          stopping early)
// Returns: one BoundaryLayerProfile per entry of 'wall_faces', in the same order
std::vector<BoundaryLayerProfile> compute_boundary_layer_profiles_point_location(
    const UnstructuredMesh& mesh, const std::vector<int>& wall_faces,
    const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& boundary_u, const std::vector<double>& boundary_v,
    double u_edge, double max_distance, int n_samples);

// Per-wall-mesh-node reduction of compute_wall_traction()/
// compute_boundary_layer_profiles()'s per-face samples, so a plotted surface
// distribution lines up with actual mesh-vertex coordinates (for arc length
// or x-coordinate) rather than a set of face-midpoint samples offset from
// them. An interior wall node (shared by exactly two wall faces on the same
// patch) gets the plain average of its two adjacent faces' values; an
// endpoint node of an open wall patch (one adjacent wall face) takes that
// face's value directly; a node shared by two DIFFERENT wall patches (a
// corner) is reported once per patch it belongs to, not silently merged
// (mirroring the force/moment totals being kept per patch).
struct WallNodeSample {
    int node_index = -1;
    int patch_id = -1;
    double x = 0.0, y = 0.0;
    double tau_wall = 0.0, Cf = 0.0, Cp = 0.0, y_plus = 0.0;
    double delta_99 = 0.0, displacement_thickness = 0.0, momentum_thickness = 0.0, shape_factor = 0.0;
};

// Input:
//   mesh          - mesh the samples/profiles were computed on
//   traction       - compute_wall_traction()'s output
//   ref            - reference quantities for Cf/Cp (see WallReferenceQuantities)
//   profiles       - compute_boundary_layer_profiles()'s output, matched to
//                    'traction' by face_index; pass an empty vector if
//                    boundary-layer thickness wasn't computed (delta_99/
//                    displacement_thickness/momentum_thickness/shape_factor
//                    are left at 0 on every resulting WallNodeSample)
// Returns: one WallNodeSample per (wall node, patch) pair present among
//          'traction' (see the corner-node convention above); order is
//          otherwise unspecified
std::vector<WallNodeSample> average_wall_samples_to_nodes(const UnstructuredMesh& mesh,
                                                            const std::vector<WallFaceSample>& traction,
                                                            const WallReferenceQuantities& ref,
                                                            const std::vector<BoundaryLayerProfile>& profiles);

#endif // WALLTRACTION_H_INCLUDED
