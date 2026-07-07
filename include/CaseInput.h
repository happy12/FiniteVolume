// SPDX-License-Identifier: GPL-3.0-only
#ifndef CASEINPUT_H_INCLUDED
#define CASEINPUT_H_INCLUDED

#include <string>
#include <vector>

#include "UnstructuredMesh.h"
#include "CflRamp.h"
#include "EulerFVMSolver.h"
#include "NavierStokesFVMSolver.h"
#include "RANSTurbulenceSASolver.h"
#include "RANSTurbulenceSSTSolver.h"
#include "SSTKOmega.h"
#include "GradientReconstruction.h"

// A boundary condition assignment read from the case file, matched by patch
// name against the boundary patches parsed from the mesh file. Used when
// 'equation' is Diffusion.
struct BoundaryConditionSpec {
    std::string patch_name;
    BoundaryType type;
    double value;
};

// A boundary condition assignment read from the case file, matched by patch
// name against the boundary patches parsed from the mesh file. Used when
// 'equation' is Euler.
struct EulerBoundaryConditionSpec {
    std::string patch_name;
    EulerBoundaryType type;
    // Prescribed primitive state (density, velocity components, pressure, in
    // mesh-consistent units -- see EulerState.h); only meaningful when
    // type == Farfield.
    double rho = 0.0, u = 0.0, v = 0.0, p = 0.0;
};

// A boundary condition assignment read from the case file, matched by patch
// name against the boundary patches parsed from the mesh file. Used when
// 'equation' is NavierStokes. Distinct from EulerBoundaryConditionSpec even
// though Farfield/Outflow mean the same thing physically, since a no-slip
// wall needs an associated thermal condition an inviscid wall never does
// (see NSBoundaryType/NSBoundaryCondition in NavierStokesFVMSolver.h).
struct NSBoundaryConditionSpec {
    std::string patch_name;
    NSBoundaryType type;
    double wall_u = 0.0, wall_v = 0.0; // NoSlipWall only; the wall's own velocity (0,0 = stationary)
    bool is_isothermal_wall = false; // NoSlipWall only; false = adiabatic (zero heat flux)
    double wall_temperature = 0.0;   // NoSlipWall + is_isothermal_wall only
    // Prescribed primitive state (density, velocity components, pressure, in
    // mesh-consistent units -- see EulerState.h); only meaningful when
    // type == Farfield.
    double rho = 0.0, u = 0.0, v = 0.0, p = 0.0;
};

// A boundary condition assignment read from the case file, matched by patch
// name against the boundary patches parsed from the mesh file. Used when
// 'equation' is RANS. Distinct from NSBoundaryConditionSpec (even though a
// no-slip wall/farfield/outflow mean exactly the same thing for the mean-flow
// equations here as they do for Navier-Stokes) purely so a case file's
// ransSA_* boundary keywords land in their own vector -- see
// RANSBoundaryConditionSA in RANSTurbulenceSASolver.h -- and to carry the one extra
// value RANS's nut transport equation needs beyond that: a prescribed
// freestream nut for a Farfield boundary.
struct RANSBoundaryConditionSpecSA {
    std::string patch_name;
    NSBoundaryType type;
    double wall_u = 0.0, wall_v = 0.0; // NoSlipWall only; the wall's own velocity (0,0 = stationary)
    bool is_isothermal_wall = false; // NoSlipWall only; false = adiabatic (zero heat flux)
    double wall_temperature = 0.0;   // NoSlipWall + is_isothermal_wall only
    // Prescribed primitive state (density, velocity components, pressure, in
    // mesh-consistent units -- see EulerState.h); only meaningful when
    // type == Farfield.
    double rho = 0.0, u = 0.0, v = 0.0, p = 0.0;
    double farfield_nut = 0.0; // only meaningful when type == Farfield; see RANSBoundaryConditionSA::farfield_nut
};

// A boundary condition assignment read from the case file, matched by patch
// name against the boundary patches parsed from the mesh file. Used when
// 'equation' is RANS_SST. Same rationale as RANSBoundaryConditionSpecSA
// (its own vector, so ransSST_* keywords don't collide with ransSA_*'s), and
// carries the two extra values the k/omega transport equations need at a
// Farfield boundary -- see RANSBoundaryConditionSST in RANSTurbulenceSSTSolver.h.
struct RANSBoundaryConditionSpecSST {
    std::string patch_name;
    NSBoundaryType type;
    double wall_u = 0.0, wall_v = 0.0; // NoSlipWall only; the wall's own velocity (0,0 = stationary)
    bool is_isothermal_wall = false; // NoSlipWall only; false = adiabatic (zero heat flux)
    double wall_temperature = 0.0;   // NoSlipWall + is_isothermal_wall only
    // Prescribed primitive state (density, velocity components, pressure, in
    // mesh-consistent units -- see EulerState.h); only meaningful when
    // type == Farfield.
    double rho = 0.0, u = 0.0, v = 0.0, p = 0.0;
    double farfield_k = 0.0, farfield_omega = 0.0; // only meaningful when type == Farfield; see RANSBoundaryConditionSST
};

// Which physics 'equation' selects.
enum class EquationSet {
    Diffusion,
    Euler,
    AdvectionDiffusion,
    NavierStokes,
    RANS_SA,
    RANS_SST
};

// Parses a simple key=value case file describing solver run parameters,
// e.g. for the scalar diffusion equation:
//   mesh_file = mesh/square.msh
//   output_file = output/result.vtk
//   alpha = 0.01
//   dt = 0.001
//   nsteps = 500
//   initial_value = 1.0
//   initial_radius = 0.25
//
//   boundary wall dirichlet 0.0
//   boundary inlet dirichlet 1.0
//
// or for the advection-diffusion equation of a passive scalar
// (equation = advection_diffusion) -- same alpha/dt/initial_value/
// initial_radius/boundary keys as diffusion above, plus:
//   u_adv = 1.0                    (uniform advection velocity components,
//   v_adv = 0.0                     mesh length units / time unit; default 0,0)
//   gradient_scheme = least-squares (or 'green-gauss'; selects
//                                    GradientCalculator's scheme -- see
//                                    GradientReconstruction.h; default least-squares)
//
// or for the compressible Euler equations (equation = euler):
//   equation = euler
//   gamma = 1.4
//   cfl = 0.5
//   cfl_mode = fixed               (or 'ramp'; default 'fixed', today's behavior
//                                    unchanged. Also applies to navier_stokes/
//                                    rans_sa/rans_sst below -- ignored for
//                                    diffusion/advection_diffusion, which use
//                                    dt/alpha directly. See docs/adaptive-cfl-
//                                    ramp-plan.md and CflRamp.h/next_cfl() for
//                                    the ramp rule. In 'ramp' mode, 'cfl' above
//                                    is the ceiling the ramp grows *toward*,
//                                    not a fixed value used from step 1.)
//   cfl_min = 0.05                 (ramp only: starting/floor CFL)
//   cfl_ramp_growth = 1.5          (ramp only: multiplier applied when the
//                                    monitored residual's windowed trend is decreasing)
//   cfl_ramp_shrink = 0.5          (ramp only: multiplier applied when the
//                                    trend is mildly rising)
//   cfl_ramp_window = 15           (ramp only: number of steps of residual
//                                    history the trend is computed over)
//   cfl_ramp_divergence_threshold = 2.0  (ramp only: log10(residual) rise over
//                                    the window that forces a hard reset to cfl_min)
//   flux_scheme = rusanov      (or 'hllc'/'exact'; omit for the default 'rusanov')
//   exact_riemann_tol = 1e-6       (exact solver's Newton-Raphson tolerance; not meant to be tuned)
//   exact_riemann_max_iter = 20    (exact solver's Newton-Raphson iteration cap; not meant to be tuned)
//   nsteps = 500
//
//   euler_init = freestream 1.0 0.0 0.0 1.0
//   (or) euler_init = tworegion 1.0 0.0 0.0 1.0  0.125 0.0 0.0 0.1  0.0
//
//   boundary wall1 wall
//   boundary inlet farfield 1.0 2.0 0.0 1.0
//   boundary outlet outflow
//
// or for the compressible Navier-Stokes equations (equation = navier_stokes)
// -- same gamma/cfl/flux_scheme/exact_riemann_*/nsteps keys as Euler above,
// plus:
//   mu = 0.01                      (dynamic viscosity, mesh-consistent units; default 0)
//   prandtl = 0.72                 (Prandtl number, dimensionless; default 0.72, air's value)
//   gas_constant = 1.0              (specific gas constant R in p = rho*R*T,
//                                    mesh-consistent units; only used for
//                                    temperature/heat conduction; default 1.0)
//   gradient_scheme = least-squares (or 'green-gauss'; same key as
//                                    advection-diffusion's, selects
//                                    GradientCalculator's scheme)
//
//   ns_init = freestream 1.0 0.0 0.0 1.0      (same grammar as euler_init, own storage)
//   (or) ns_init = tworegion 1.0 0.0 0.0 1.0  0.125 0.0 0.0 0.1  0.0
//
//   boundary wall1 ns_wall                       (no-slip, stationary, adiabatic)
//   boundary wall2 ns_wall_isothermal 1.0         (no-slip, stationary, fixed wall temperature)
//   boundary wall3 ns_wall_moving 0.5 0.0         (no-slip, moving at (u,v)=(0.5,0), adiabatic -- e.g. Couette flow)
//   boundary wall4 ns_wall_moving_isothermal 0.5 0.0 1.0  (moving + fixed wall temperature)
//   boundary inlet ns_farfield 1.0 2.0 0.0 1.0
//   boundary outlet ns_outflow
//
//   resolution_report_file = output/resolution.csv   (Navier-Stokes only;
//                                    omit to disable; see
//                                    NavierStokesFVMSolver::compute_resolution_diagnostics()
//                                    -- a resolution-adequacy heuristic for
//                                    treating a run as fully-resolved 2D
//                                    unsteady, NOT a 3D DNS claim)
//   resolution_report_interval = 10  (row written every N steps, default 1)
//
//   wall_forces_file = output/wall_forces.csv    (Navier-Stokes / RANS only;
//                                    omit to disable; see WallTraction.h --
//                                    per-(step, patch) friction/pressure/total
//                                    drag, lift, and moment, plus their
//                                    coefficients, over every NoSlipWall patch)
//   wall_forces_interval = 10        (row block written every N steps, default 1)
//   wall_profile_file = output/wall_profile.csv  (Navier-Stokes / RANS only;
//                                    omit to disable; see WallTraction.h --
//                                    one row per wall mesh node: Cf, Cp, y+,
//                                    and boundary-layer thickness estimates)
//   wall_profile_interval = 10       (snapshot written every N steps; default
//                                    0 = write once at the run's natural end,
//                                    like write_interval = 0's VTK convention)
//
//   reference_density = 1.0          (freestream rho for dynamic pressure;
//                                    omit to default to the first ns_farfield
//                                    patch's rho -- only resolved/validated
//                                    when wall_forces_file or
//                                    wall_profile_file is set)
//   reference_velocity_x = 1.0       (freestream velocity components; drag
//   reference_velocity_y = 0.0        direction = this vector normalized;
//                                    omit either to default from the first
//                                    ns_farfield patch, same rule as above)
//   reference_pressure = 0.0         (freestream p for Cp; omit to default
//                                    from the first ns_farfield patch, same
//                                    rule as above)
//   reference_length = 1.0           (length scale for Cd/Cl/Cm, e.g. chord;
//                                    1.0 degrades gracefully to a
//                                    per-unit-span coefficient)
//   moment_reference_x = 0.0         (point Cm is taken about, e.g. an
//   moment_reference_y = 0.0          airfoil's quarter-chord)
//   boundary_layer_method = marching (or 'point-location'; selects which of
//                                    compute_boundary_layer_profiles()/
//                                    compute_boundary_layer_profiles_point_location()
//                                    (WallTraction.h) wall_profile_file's
//                                    thickness columns use -- point-location
//                                    is more expensive but has no near-wall
//                                    cell-stacking assumption; default 'marching')
//   boundary_layer_max_distance = 0.05 (point-location only; furthest
//                                    wall-normal distance to sample; omit or
//                                    <= 0 for an automatic default -- the
//                                    mesh's own bounding-box diagonal)
//   boundary_layer_n_samples = 200   (point-location only; number of
//                                    evenly-spaced sample points, default 200)
//
// or for the RANS (Spalart-Allmaras) equations (equation = rans_sa) -- same
// gamma/cfl/flux_scheme/exact_riemann_*/mu/prandtl/gas_constant/
// gradient_scheme/nsteps keys as Navier-Stokes above (mu/prandtl/gas_constant
// govern the mean-flow molecular terms exactly as they do for Navier-Stokes),
// plus:
//   prandtl_t = 0.9                 (turbulent Prandtl number, dimensionless; default 0.9, air's value)
//   initial_nut = 3e-5               (uniform initial nut applied to every cell,
//                                    mesh-consistent units; default 0 -- SA is
//                                    documented as sensitive to this choice)
//   sa_cb1 = 0.1355                 (SA-noft2 model constants -- see
//   sa_cb2 = 0.622                   SpalartAllmaras.h; each independently
//   sa_sigma = 0.6667                optional, defaults to the standard set;
//   sa_kappa = 0.41                  exposed for experimentation, not meant
//   sa_cw2 = 0.3                     to be tuned in normal use)
//   sa_cw3 = 2.0
//   sa_cv1 = 7.1
//   sa_cv2 = 0.7
//   sa_cv3 = 0.9
//
//   ransSA_init = freestream 1.0 0.0 0.0 1.0     (same grammar as ns_init, own storage)
//   (or) ransSA_init = tworegion 1.0 0.0 0.0 1.0  0.125 0.0 0.0 0.1  0.0
//
//   boundary wall1 ransSA_wall                       (no-slip, stationary, adiabatic)
//   boundary wall2 ransSA_wall_isothermal 1.0         (no-slip, stationary, fixed wall temperature)
//   boundary wall3 ransSA_wall_moving 0.5 0.0         (no-slip, moving at (u,v)=(0.5,0), adiabatic)
//   boundary wall4 ransSA_wall_moving_isothermal 0.5 0.0 1.0  (moving + fixed wall temperature)
//   boundary inlet ransSA_farfield 1.0 2.0 0.0 1.0 3e-5   (rho u v p farfield_nut)
//   boundary outlet ransSA_outflow
//
// wall_forces_file/wall_profile_file/reference_*/boundary_layer_* (see
// Navier-Stokes above) all apply to RANS too, using mu + rho*nu_t as the
// effective viscosity instead of just mu.
//
// or for the RANS (k-omega SST) equations (equation = rans_sst) -- same
// gamma/cfl/flux_scheme/exact_riemann_*/mu/prandtl/gas_constant/
// gradient_scheme/nsteps/prandtl_t keys as RANS (Spalart-Allmaras) above,
// plus:
//   initial_k = 1e-4                (uniform initial k applied to every cell,
//                                    mesh-consistent units; default 0)
//   initial_omega = 10.0             (uniform initial omega applied to every
//                                    cell, mesh-consistent units; default 0)
//   sst_turbulence_intensity = 0.05 (Tu, dimensionless; used to DERIVE
//   sst_eddy_viscosity_ratio = 10.0   initial_k/initial_omega from the
//                                    ransSST_init freestream velocity and
//                                    mu/rho -- k = 1.5*(Tu*|U|)^2,
//                                    omega = k/(ratio*nu) -- ONLY when
//                                    initial_k/initial_omega are not set
//                                    explicitly above; omitted entirely,
//                                    both default to 0, same as SA's
//                                    initial_nut having no derivation)
//   sst_limiter = vorticity          (or 'strain_rate'; selects
//                                    SSTLimiterVariant -- see SSTKOmega.h;
//                                    default 'vorticity', the original 1994 SST)
//   sst_kato_launder = false         (Kato-Launder production limiter, true/false; see
//                                    SSTKOmega.h's compute_sst_source_terms(); default false)
//   sst_beta_star = 0.09             (SST model constants -- see SSTKOmega.h;
//   sst_kappa = 0.41                  each independently optional, defaults
//   sst_a1 = 0.31                     to the NASA TMR set; exposed for
//   sst_sigma_k1 = 0.85               experimentation, not meant to be tuned
//   sst_sigma_k2 = 1.0                in normal use)
//   sst_sigma_omega1 = 0.5
//   sst_sigma_omega2 = 0.856
//   sst_beta1 = 0.075
//   sst_beta2 = 0.0828
//
//   ransSST_init = freestream 1.0 0.0 0.0 1.0     (same grammar as ransSA_init, own storage)
//   (or) ransSST_init = tworegion 1.0 0.0 0.0 1.0  0.125 0.0 0.0 0.1  0.0
//
//   boundary wall1 ransSST_wall                       (no-slip, stationary, adiabatic)
//   boundary wall2 ransSST_wall_isothermal 1.0         (no-slip, stationary, fixed wall temperature)
//   boundary wall3 ransSST_wall_moving 0.5 0.0         (no-slip, moving at (u,v)=(0.5,0), adiabatic)
//   boundary wall4 ransSST_wall_moving_isothermal 0.5 0.0 1.0  (moving + fixed wall temperature)
//   boundary inlet ransSST_farfield 1.0 2.0 0.0 1.0 1e-4 10.0   (rho u v p farfield_k farfield_omega,
//                                    explicit per-patch values -- NOT derived
//                                    from sst_turbulence_intensity/
//                                    sst_eddy_viscosity_ratio, unlike
//                                    initial_k/initial_omega above)
//   boundary outlet ransSST_outflow
//
// wall_forces_file/wall_profile_file/reference_*/boundary_layer_* (see
// Navier-Stokes above) all apply to RANS (k-omega SST) too, using
// mu + rho*nu_t as the effective viscosity instead of just mu, same as SA.
//
// residual_tolerance_k/residual_tolerance_omega (see the stopping-criteria
// keys below) are this equation set's analogue of RANS (Spalart-Allmaras)'s
// residual_tolerance_nut.
//
// All equation sets also accept optional monitoring keys:
//   residual_file = output/residual.csv   (omit to disable residual tracking)
//   residual_interval = 10                (rows written every N steps, default 1)
//   write_interval = 100                  (numbered VTK snapshots every N steps, default 0 = off)
//
// ...and an optional performance key:
//   num_threads = 4                       (OpenMP thread count; omit/0 = OpenMP's own default,
//                                           which respects the OMP_NUM_THREADS environment variable)
//
// ...and optional stopping-criteria / checkpoint keys (nsteps is always the
// absolute target step count -- raising it, or loosening/tightening a
// residual tolerance, and rerunning the same command is how a stopped run is
// continued further):
//   residual_tolerance = 1e-6             (diffusion / advection-diffusion only; omit = no residual-based stop)
//
//   residual_tolerance_rho = 1e-6         (Euler / Navier-Stokes / RANS only; each independently optional --
//   residual_tolerance_rho_u = 1e-6        "converged" requires every tolerance the user set
//   residual_tolerance_rho_v = 1e-6        to be satisfied simultaneously)
//   residual_tolerance_E = 1e-6
//   residual_tolerance_nut = 1e-6         (RANS (Spalart-Allmaras) only; checked alongside the four above)
//   residual_tolerance_k = 1e-6           (RANS (k-omega SST) only; checked alongside the four above)
//   residual_tolerance_omega = 1e-6       (RANS (k-omega SST) only; checked alongside the four above)
//
//   checkpoint_file = output/checkpoint.bin  (omit to disable checkpointing; written when
//                                              the run stops for any reason except divergence;
//                                              auto-resumed from if this file already exists)
//
// ...and optional output keys (shared by all equation sets):
//   scratch_dir = output                  (base directory for relative output_file/checkpoint_file/
//                                           residual_file paths; omit to resolve relative to the
//                                           working directory as before; absolute paths and mesh_file
//                                           are never affected)
//   output_precision = 6                  (significant digits written for every output double --
//                                           VTK results/snapshots, residual_file CSV, and
//                                           FvMeshWriter; valid range 1-17)
//
// A run also always stops immediately (regardless of the above) if its
// residual goes NaN/Inf, which is treated as an unrecoverable divergence: a
// final VTK snapshot is still written for inspection, but no checkpoint is,
// and the program exits with a non-zero status.
class CaseInput {
public:
    // Loads and parses a case file. See the class comment above for the
    // expected format.
    // Input:   filename - path to the case file
    // Output:  *this is populated with the parsed run parameters (mesh_file,
    //          output_file, nsteps, equation, and whichever of the
    //          diffusion-only or Euler-only fields apply; fields not present
    //          in the file keep their default value)
    // Returns: true on success (file opened, mesh_file/output_file set,
    //          every boundary/equation/euler_init/ns_init keyword recognized,
    //          every numeric "key = value" parsed as a number, and any value
    //          with an obvious valid domain -- nsteps/dt/cfl/mu/gamma/
    //          output_precision -- in range); false otherwise (a descriptive
    //          message naming the file, line, and key is printed to stderr).
    bool load(const std::string& filename);

    std::string mesh_file;                 // Path to the mesh file to solve on: Gmsh ASCII (.msh) or
                                            // this project's own FVMESH format (.fvmesh, see
                                            // docs/fvmesh-format.md), dispatched by extension
    std::string output_file;                // Path to write the result .vtk file to
    int nsteps = 500;                       // Number of time steps to advance
    EquationSet equation = EquationSet::Diffusion; // Which physics to solve

    // --- Monitoring (shared by both equation sets) ---
    std::string residual_file;              // Path to write a per-step residual history CSV to; empty = disabled
    int residual_interval = 1;              // Write a residual_file row every N steps (only used if residual_file is set)
    int write_interval = 0;                 // Write a numbered VTK snapshot every N steps; 0 = disabled (only the final output_file write)

    // --- Performance (shared by both equation sets) ---
    int num_threads = 0;                    // OpenMP thread count, dimensionless; 0 = OpenMP's own default

    // --- Stopping criteria (shared by both equation sets; a negative value means disabled,
    //     since a real tolerance is always positive) ---
    double residual_tolerance = -1.0;        // diffusion only
    double residual_tolerance_rho = -1.0;    // Euler only, independent per conserved variable
    double residual_tolerance_rho_u = -1.0;
    double residual_tolerance_rho_v = -1.0;
    double residual_tolerance_E = -1.0;
    double residual_tolerance_nut = -1.0;    // RANS (Spalart-Allmaras) only, checked alongside rho/rho_u/rho_v/E above
    double residual_tolerance_k = -1.0;      // RANS (k-omega SST) only, checked alongside rho/rho_u/rho_v/E above
    double residual_tolerance_omega = -1.0;  // RANS (k-omega SST) only, checked alongside rho/rho_u/rho_v/E above

    // --- Checkpoint/restart (shared by both equation sets) ---
    std::string checkpoint_file;             // Path to save/resume solver state; empty = disabled

    // --- Output (shared by both equation sets) ---
    std::string scratch_dir;                 // Base directory for relative output_file/checkpoint_file/
                                              // residual_file paths; empty = resolve relative to the
                                              // working directory (unchanged behavior). Absolute paths
                                              // and mesh_file are never affected by this.
    int output_precision = 6;                // Significant digits written for every output double
                                              // (VTK results/snapshots, residual_file CSV, FvMeshWriter);
                                              // valid range 1-17 (17 = std::numeric_limits<double>::
                                              // max_digits10, the ceiling beyond which extra digits are
                                              // meaningless for a double)

    // --- Diffusion / Advection-Diffusion parameters (shared, since a scalar
    //     Dirichlet/Neumann BC means the same thing for both) ---
    double alpha = 0.01;                    // Diffusion coefficient, in (mesh length units)^2 / (time units)
    double dt = 0.001;                      // Time step size, in time units
    double initial_value = 1.0;             // phi value inside the initial condition disc
    double initial_radius = 0.25;           // Radius (from the origin) of the initial condition disc, in mesh length units
    std::vector<BoundaryConditionSpec> boundary_conditions; // Per-patch boundary condition assignments

    // --- Advection-Diffusion-only parameters ---
    double u_adv = 0.0, v_adv = 0.0;        // Uniform advection velocity components, in mesh length units / time unit
    GradientScheme gradient_scheme = GradientScheme::LeastSquares; // GradientCalculator scheme used to reconstruct cell gradients each step

    // --- Euler-only parameters ---
    double gamma = 1.4;                     // Ratio of specific heats, dimensionless (1.4 for air)
    double cfl = 0.5;                       // CFL number, dimensionless, controlling the adaptive time step;
                                              // in cfl_mode == Ramp, the ceiling the ramp grows toward
    CflMode cfl_mode = CflMode::Fixed;      // Fixed (unchanged behavior) or Ramp -- see CflRamp.h; also
                                              // applies to navier_stokes/rans_sa/rans_sst below
    CflRampParams cfl_ramp;                 // Ramp-only parameters (cfl_min/cfl_ramp_growth/cfl_ramp_shrink/
                                              // cfl_ramp_window/cfl_ramp_divergence_threshold); ignored when
                                              // cfl_mode == Fixed
    NumericalFluxScheme flux_scheme = NumericalFluxScheme::Rusanov; // Numerical flux used at every face (see EulerFVMSolver.h)
    // Exact Riemann solver's Newton-Raphson relative pressure tolerance and
    // iteration cap (see ExactRiemannFlux.h); ignored unless flux_scheme is
    // 'exact'. Exposed here for transparency, not because they are meant to
    // be tuned in normal use -- the defaults below are the values this
    // solver always used before they became configurable.
    double exact_riemann_tol = 1e-6;
    int exact_riemann_max_iter = 20;
    EulerInitialCondition euler_ic;         // Initial condition, in primitive variables (see EulerFVMSolver.h)
    std::vector<EulerBoundaryConditionSpec> euler_boundary_conditions; // Per-patch boundary condition assignments

    // --- Navier-Stokes-only parameters ---
    double mu = 0.0;                        // Dynamic viscosity, mesh-consistent units; 0 = no viscous stress
    double prandtl = 0.72;                  // Prandtl number, dimensionless (0.72 for air)
    double gas_constant = 1.0;              // Specific gas constant R in p = rho*R*T, mesh-consistent units
    EulerInitialCondition ns_ic;             // Initial condition, in primitive variables (own storage, same grammar as euler_ic)
    std::vector<NSBoundaryConditionSpec> ns_boundary_conditions; // Per-patch boundary condition assignments

    // --- RANS (Spalart-Allmaras)-only parameters --- (reuses gamma/cfl/
    //     flux_scheme/gradient_scheme/exact_riemann_*/mu/prandtl/gas_constant above)
    double prandtl_t = 0.9;                 // Turbulent Prandtl number, dimensionless (~0.9 for air); no Navier-Stokes analog
    double initial_nut = 0.0;               // Uniform initial nut applied to every cell; SA is documented as sensitive to this choice
    EulerInitialCondition ransSA_ic;           // Mean-flow initial condition, in primitive variables (own storage, same grammar as ns_ic)
    std::vector<RANSBoundaryConditionSpecSA> ransSA_boundary_conditions; // Per-patch boundary condition assignments
    // SA-noft2 model constants (see SpalartAllmaras.h); default to the
    // standard set. Exposed for experimentation, not because they are meant
    // to be tuned in normal use.
    SAModelConstants sa_constants;

    // --- RANS (k-omega SST)-only parameters --- (reuses gamma/cfl/
    //     flux_scheme/gradient_scheme/exact_riemann_*/mu/prandtl/gas_constant/
    //     prandtl_t above)
    double initial_k = 0.0, initial_omega = 0.0; // Uniform initial k/omega applied to every cell; see class comment
    bool initial_k_set = false, initial_omega_set = false; // Recorded so sst_turbulence_intensity/
                                              // sst_eddy_viscosity_ratio only derive a value the case
                                              // file didn't set explicitly -- same *_set convention as
                                              // reference_density_set below.
    double sst_turbulence_intensity = 0.0;   // Tu, dimensionless; only used to derive initial_k/initial_omega
    double sst_eddy_viscosity_ratio = 0.0;   // nu_t/nu, dimensionless; only used to derive initial_k/initial_omega
    SSTLimiterVariant sst_limiter_variant = SSTLimiterVariant::Vorticity; // see SSTKOmega.h
    bool sst_kato_launder = false;           // see SSTKOmega.h's compute_sst_source_terms()
    EulerInitialCondition ransSST_ic;          // Mean-flow initial condition, in primitive variables (own storage, same grammar as ransSA_ic)
    std::vector<RANSBoundaryConditionSpecSST> ransSST_boundary_conditions; // Per-patch boundary condition assignments
    // SST model constants (see SSTKOmega.h); default to the NASA TMR set.
    // Exposed for experimentation, not because they are meant to be tuned in
    // normal use.
    SSTModelConstants sst_constants;

    // Resolution-adequacy diagnostic (see NavierStokesFVMSolver::compute_resolution_diagnostics()):
    // NOT a claim about 3D DNS -- this solver is 2D-only, and classical
    // Kolmogorov scaling is a 3D result. A pragmatic resolution heuristic
    // for a fully-resolved unsteady 2D run, nothing more.
    std::string resolution_report_file;      // Path to write a resolution-diagnostic CSV to; empty = disabled
    int resolution_report_interval = 1;      // Write a resolution_report_file row every N steps

    // --- Wall diagnostics (Navier-Stokes / RANS only; see WallTraction.h and
    //     docs/wall-diagnostics-plan.md) ---
    std::string wall_forces_file;             // Path to write a per-(step, patch) force/moment CSV to; empty = disabled
    int wall_forces_interval = 1;             // Write a wall_forces_file row block every N steps
    std::string wall_profile_file;            // Path to write a per-wall-node Cf/Cp/y+/boundary-layer-thickness CSV to; empty = disabled
    int wall_profile_interval = 0;            // Write a wall_profile_file snapshot every N steps; 0 = write once at the run's natural end

    // Reference/freestream quantities used to nondimensionalize Cf/Cp/Cd/Cl/Cm
    // (see WallReferenceQuantities, WallTraction.h). Each *_set flag records
    // whether the case file explicitly set the corresponding value; when
    // false, CaseInput::load() defaults it from the first ns_farfield (or, for
    // equation == RANS_SA, ransSA_farfield) patch's prescribed state -- resolved
    // (and validated as an error if no such patch exists to default from)
    // only when wall_forces_file or wall_profile_file is actually set, so an
    // unrelated case file is never broken by these keys' absence.
    double reference_density = 0.0;
    bool reference_density_set = false;
    double reference_velocity_x = 0.0;
    bool reference_velocity_x_set = false;
    double reference_velocity_y = 0.0;
    bool reference_velocity_y_set = false;
    double reference_pressure = 0.0;
    bool reference_pressure_set = false;
    double reference_length = 1.0;            // Length scale for Cd/Cl/Cm (e.g. chord); 1.0 = per-unit-span coefficient
    double moment_reference_x = 0.0, moment_reference_y = 0.0; // Point Cm is taken about

    // Which of compute_boundary_layer_profiles()/
    // compute_boundary_layer_profiles_point_location() (WallTraction.h)
    // wall_profile_file's boundary-layer thickness columns use. Point-location
    // is more expensive (O(cells) per sample point) but has no near-wall
    // cell-stacking assumption -- see docs/wall-diagnostics-plan.md's Phase 2
    // stress-test finding for why marching can fail on some meshes.
    BoundaryLayerMethod boundary_layer_method = BoundaryLayerMethod::Marching;
    double boundary_layer_max_distance = -1.0; // point-location only; <= 0 = auto (mesh bounding-box diagonal)
    int boundary_layer_n_samples = 200;        // point-location only
};

#endif // CASEINPUT_H_INCLUDED
