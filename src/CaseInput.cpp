// SPDX-License-Identifier: GPL-3.0-only
#include "CaseInput.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>

namespace {

// Strips leading/trailing whitespace (spaces, tabs, CR, LF) from a string.
// Input:  s - the string to trim
// Returns: s with surrounding whitespace removed, or "" if s is all whitespace
std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

} // namespace

// See CaseInput.h for the file format. Each line is either:
//   - a "boundary <patch_name> <type_word> <value...>" entry, or
//   - a "key = value" run parameter, or
//   - blank/a comment (from '#' to end of line), which is ignored.
//
// Input:
//   filename - path to the case file to parse
// Output:
//   *this    - populated from the file's keys (see CaseInput.h's field list);
//              fields not present in the file keep their default value
// Returns:
//   true if the file was opened, both mesh_file and output_file were set, and
//   every boundary/equation/euler_init keyword was recognized; false
//   otherwise (a descriptive message is printed to stderr for an
//   unrecognized keyword).
bool CaseInput::load(const std::string& filename) {
    std::ifstream in(filename);
    if (!in.is_open()) {
        return false;
    }

    std::string line;
    while (std::getline(in, line)) {
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }
        line = trim(line);
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string first;
        iss >> first;

        // "boundary <patch_name> <type_word> <value...>". type_word is
        // "dirichlet"/"neumann" (diffusion/advection-diffusion, one trailing
        // value; see BoundaryType in UnstructuredMesh.h), "wall"/"farfield"/
        // "outflow" (Euler; see EulerBoundaryType in EulerFVMSolver.h --
        // "farfield" takes 4 trailing values "rho u v p", the others take
        // none), or "ns_wall"/"ns_wall_isothermal"/"ns_wall_moving"/
        // "ns_wall_moving_isothermal"/"ns_farfield"/"ns_outflow"
        // (Navier-Stokes; see NSBoundaryType in NavierStokesFVMSolver.h --
        // prefixed since "farfield"/"outflow" mean the same thing for Euler
        // and Navier-Stokes but need to land in different spec vectors;
        // "ns_wall_isothermal" takes 1 trailing value "wall_temperature",
        // "ns_wall_moving" takes 2 "wall_u wall_v", "ns_wall_moving_isothermal"
        // takes 3 "wall_u wall_v wall_temperature", "ns_farfield" takes 4
        // "rho u v p", "ns_wall"/"ns_outflow" take none).
        if (first == "boundary") {
            std::string patch_name, type_word;
            iss >> patch_name >> type_word;

            if (type_word == "dirichlet" || type_word == "neumann") {
                BoundaryConditionSpec bc;
                bc.patch_name = patch_name;
                bc.type = (type_word == "dirichlet") ? BoundaryType::Dirichlet : BoundaryType::Neumann;
                iss >> bc.value;
                boundary_conditions.push_back(bc);
            } else if (type_word == "wall" || type_word == "farfield" || type_word == "outflow") {
                EulerBoundaryConditionSpec bc;
                bc.patch_name = patch_name;
                if (type_word == "wall") bc.type = EulerBoundaryType::Wall;
                else if (type_word == "outflow") bc.type = EulerBoundaryType::Outflow;
                else {
                    bc.type = EulerBoundaryType::Farfield;
                    iss >> bc.rho >> bc.u >> bc.v >> bc.p;
                }
                euler_boundary_conditions.push_back(bc);
            } else if (type_word == "ns_wall" || type_word == "ns_wall_isothermal" ||
                       type_word == "ns_wall_moving" || type_word == "ns_wall_moving_isothermal" ||
                       type_word == "ns_farfield" || type_word == "ns_outflow") {
                NSBoundaryConditionSpec bc;
                bc.patch_name = patch_name;
                if (type_word == "ns_wall") {
                    bc.type = NSBoundaryType::NoSlipWall;
                } else if (type_word == "ns_wall_isothermal") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    bc.is_isothermal_wall = true;
                    iss >> bc.wall_temperature;
                } else if (type_word == "ns_wall_moving") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    iss >> bc.wall_u >> bc.wall_v;
                } else if (type_word == "ns_wall_moving_isothermal") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    bc.is_isothermal_wall = true;
                    iss >> bc.wall_u >> bc.wall_v >> bc.wall_temperature;
                } else if (type_word == "ns_outflow") {
                    bc.type = NSBoundaryType::Outflow;
                } else {
                    bc.type = NSBoundaryType::Farfield;
                    iss >> bc.rho >> bc.u >> bc.v >> bc.p;
                }
                ns_boundary_conditions.push_back(bc);
            } else if (type_word == "rans_wall" || type_word == "rans_wall_isothermal" ||
                       type_word == "rans_wall_moving" || type_word == "rans_wall_moving_isothermal" ||
                       type_word == "rans_farfield" || type_word == "rans_outflow") {
                RANSBoundaryConditionSpec bc;
                bc.patch_name = patch_name;
                if (type_word == "rans_wall") {
                    bc.type = NSBoundaryType::NoSlipWall;
                } else if (type_word == "rans_wall_isothermal") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    bc.is_isothermal_wall = true;
                    iss >> bc.wall_temperature;
                } else if (type_word == "rans_wall_moving") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    iss >> bc.wall_u >> bc.wall_v;
                } else if (type_word == "rans_wall_moving_isothermal") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    bc.is_isothermal_wall = true;
                    iss >> bc.wall_u >> bc.wall_v >> bc.wall_temperature;
                } else if (type_word == "rans_outflow") {
                    bc.type = NSBoundaryType::Outflow;
                } else {
                    bc.type = NSBoundaryType::Farfield;
                    iss >> bc.rho >> bc.u >> bc.v >> bc.p >> bc.farfield_nut;
                }
                rans_boundary_conditions.push_back(bc);
            } else {
                std::cerr << "Error in case file '" << filename << "': unknown boundary type '"
                           << type_word << "' for patch '" << patch_name
                           << "' (expected 'dirichlet', 'neumann', 'wall', 'farfield', 'outflow', 'ns_wall', "
                           << "'ns_wall_isothermal', 'ns_wall_moving', 'ns_wall_moving_isothermal', "
                           << "'ns_farfield', 'ns_outflow', 'rans_wall', 'rans_wall_isothermal', "
                           << "'rans_wall_moving', 'rans_wall_moving_isothermal', 'rans_farfield' or 'rans_outflow')\n";
                return false;
            }
            continue;
        }

        // Otherwise expect a "key = value" run parameter.
        size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;

        std::string key = trim(line.substr(0, eq_pos));
        std::string value = trim(line.substr(eq_pos + 1));

        if (key == "mesh_file") mesh_file = value;                    // Path to the .msh or .fvmesh file to load
        else if (key == "output_file") output_file = value;           // Path to write the result .vtk file to
        else if (key == "alpha") alpha = std::stod(value);            // Diffusion coefficient, in (mesh length units)^2 / (time units)
        else if (key == "dt") dt = std::stod(value);                  // Time step size, in time units (must satisfy the explicit-scheme stability limit)
        else if (key == "nsteps") nsteps = std::stoi(value);          // Number of time steps to advance
        else if (key == "initial_value") initial_value = std::stod(value);   // phi value inside the initial condition disc
        else if (key == "initial_radius") initial_radius = std::stod(value); // Radius (from the origin) of the initial condition disc, in mesh length units
        else if (key == "u_adv") u_adv = std::stod(value);            // Uniform advection velocity x-component (advection-diffusion)
        else if (key == "v_adv") v_adv = std::stod(value);            // Uniform advection velocity y-component (advection-diffusion)
        else if (key == "gamma") gamma = std::stod(value);            // Ratio of specific heats, dimensionless (Euler / Navier-Stokes)
        else if (key == "cfl") cfl = std::stod(value);                // CFL number, dimensionless, controlling the adaptive time step (Euler / Navier-Stokes)
        else if (key == "mu") mu = std::stod(value);                  // Dynamic viscosity, mesh-consistent units (Navier-Stokes)
        else if (key == "prandtl") prandtl = std::stod(value);        // Prandtl number, dimensionless (Navier-Stokes)
        else if (key == "gas_constant") gas_constant = std::stod(value); // Specific gas constant R in p = rho*R*T (Navier-Stokes / RANS)
        else if (key == "prandtl_t") prandtl_t = std::stod(value);        // Turbulent Prandtl number, dimensionless (RANS)
        else if (key == "initial_nut") initial_nut = std::stod(value);    // Uniform initial nut applied to every cell (RANS)
        else if (key == "sa_cb1") sa_constants.cb1 = std::stod(value);     // SA-noft2 model constants (RANS) -- see SpalartAllmaras.h
        else if (key == "sa_cb2") sa_constants.cb2 = std::stod(value);
        else if (key == "sa_sigma") sa_constants.sigma = std::stod(value);
        else if (key == "sa_kappa") sa_constants.kappa = std::stod(value);
        else if (key == "sa_cw2") sa_constants.cw2 = std::stod(value);
        else if (key == "sa_cw3") sa_constants.cw3 = std::stod(value);
        else if (key == "sa_cv1") sa_constants.cv1 = std::stod(value);
        else if (key == "sa_cv2") sa_constants.cv2 = std::stod(value);
        else if (key == "sa_cv3") sa_constants.cv3 = std::stod(value);
        else if (key == "resolution_report_file") resolution_report_file = value; // Path to resolution-diagnostic CSV; empty = disabled (Navier-Stokes)
        else if (key == "resolution_report_interval") resolution_report_interval = std::stoi(value); // Write a row every N steps
        else if (key == "wall_forces_file") wall_forces_file = value; // Path to per-(step, patch) force/moment CSV; empty = disabled (Navier-Stokes)
        else if (key == "wall_forces_interval") wall_forces_interval = std::stoi(value); // Write a row block every N steps
        else if (key == "wall_profile_file") wall_profile_file = value; // Path to per-wall-node Cf/Cp/y+/BL-thickness CSV; empty = disabled (Navier-Stokes)
        else if (key == "wall_profile_interval") wall_profile_interval = std::stoi(value); // 0 = write once at the run's natural end
        else if (key == "reference_density") { reference_density = std::stod(value); reference_density_set = true; }
        else if (key == "reference_velocity_x") { reference_velocity_x = std::stod(value); reference_velocity_x_set = true; }
        else if (key == "reference_velocity_y") { reference_velocity_y = std::stod(value); reference_velocity_y_set = true; }
        else if (key == "reference_pressure") { reference_pressure = std::stod(value); reference_pressure_set = true; }
        else if (key == "reference_length") reference_length = std::stod(value); // Length scale for Cd/Cl/Cm; 1.0 = per-unit-span
        else if (key == "moment_reference_x") moment_reference_x = std::stod(value); // Point Cm is taken about
        else if (key == "moment_reference_y") moment_reference_y = std::stod(value);
        else if (key == "boundary_layer_max_distance") boundary_layer_max_distance = std::stod(value); // point-location only; <= 0 = auto
        else if (key == "boundary_layer_n_samples") boundary_layer_n_samples = std::stoi(value);        // point-location only
        else if (key == "residual_file") residual_file = value;              // Path to the residual history CSV; empty = disabled
        else if (key == "residual_interval") residual_interval = std::stoi(value); // Write a residual row every N steps
        else if (key == "write_interval") write_interval = std::stoi(value); // Write a numbered VTK snapshot every N steps; 0 = disabled
        else if (key == "num_threads") num_threads = std::stoi(value);       // OpenMP thread count, dimensionless; 0 = OpenMP's own default
        else if (key == "residual_tolerance") residual_tolerance = std::stod(value);           // diffusion only; < 0 = disabled
        else if (key == "residual_tolerance_rho") residual_tolerance_rho = std::stod(value);   // Euler only; < 0 = disabled
        else if (key == "residual_tolerance_rho_u") residual_tolerance_rho_u = std::stod(value);
        else if (key == "residual_tolerance_rho_v") residual_tolerance_rho_v = std::stod(value);
        else if (key == "residual_tolerance_E") residual_tolerance_E = std::stod(value);
        else if (key == "residual_tolerance_nut") residual_tolerance_nut = std::stod(value); // RANS only; < 0 = disabled
        else if (key == "checkpoint_file") checkpoint_file = value;          // Path to save/resume solver state; empty = disabled
        // Exact Riemann solver's Newton-Raphson tolerance/iteration cap (see
        // ExactRiemannFlux.h); exposed for transparency, not meant to be
        // tuned in normal use.
        else if (key == "exact_riemann_tol") exact_riemann_tol = std::stod(value);
        else if (key == "exact_riemann_max_iter") exact_riemann_max_iter = std::stoi(value);
        else if (key == "scratch_dir") scratch_dir = value;                  // Base dir for relative output_file/checkpoint_file/residual_file; empty = disabled
        // "output_precision = <1-17>" -- significant digits written for every
        // output double (VTK results/snapshots, residual_file CSV, FvMeshWriter).
        else if (key == "output_precision") {
            output_precision = std::stoi(value);
            if (output_precision < 1 || output_precision > 17) {
                std::cerr << "Error in case file '" << filename << "': output_precision must be between 1 and 17 (got "
                           << output_precision << ")\n";
                return false;
            }
        }
        // "equation = diffusion|euler|advection_diffusion|navier_stokes|rans" --
        // selects which physics run_diffusion()/run_euler()/
        // run_advection_diffusion()/run_navier_stokes()/run_rans() (see
        // main.cpp) is used to solve.
        else if (key == "equation") {
            if (value == "diffusion") equation = EquationSet::Diffusion;
            else if (value == "euler") equation = EquationSet::Euler;
            else if (value == "advection_diffusion") equation = EquationSet::AdvectionDiffusion;
            else if (value == "navier_stokes") equation = EquationSet::NavierStokes;
            else if (value == "rans") equation = EquationSet::RANS;
            else {
                std::cerr << "Error in case file '" << filename << "': unknown equation '"
                           << value << "' (expected 'diffusion', 'euler', 'advection_diffusion', "
                              "'navier_stokes' or 'rans')\n";
                return false;
            }
        // "gradient_scheme = least-squares|green-gauss" -- selects
        // GradientCalculator's scheme (see GradientReconstruction.h),
        // currently used by the advection-diffusion equation set.
        } else if (key == "gradient_scheme") {
            if (value == "least-squares") gradient_scheme = GradientScheme::LeastSquares;
            else if (value == "green-gauss") gradient_scheme = GradientScheme::GreenGauss;
            else {
                std::cerr << "Error in case file '" << filename << "': unknown gradient_scheme '"
                           << value << "' (expected 'least-squares' or 'green-gauss')\n";
                return false;
            }
        // "boundary_layer_method = marching|point-location" -- selects which
        // of compute_boundary_layer_profiles()/
        // compute_boundary_layer_profiles_point_location() (WallTraction.h)
        // wall_profile_file's thickness columns use (Navier-Stokes only).
        } else if (key == "boundary_layer_method") {
            if (value == "marching") boundary_layer_method = BoundaryLayerMethod::Marching;
            else if (value == "point-location") boundary_layer_method = BoundaryLayerMethod::PointLocation;
            else {
                std::cerr << "Error in case file '" << filename << "': unknown boundary_layer_method '"
                           << value << "' (expected 'marching' or 'point-location')\n";
                return false;
            }
        // "euler_init = <mode> <values...>" -- mode is either "freestream"
        // (4 values: rho u v p) or "tworegion" (9 values: rho_l u_l v_l p_l
        // rho_r u_r v_r p_r x0). Values are primitive variables in
        // mesh-consistent units (density, velocity components, pressure; x0
        // is a mesh length-unit coordinate). See EulerInitialCondition in
        // EulerFVMSolver.h.
        // "flux_scheme = rusanov|hllc|exact" -- selects the numerical flux
        // used at every face by EulerFVMSolver::step() (see
        // NumericalFluxScheme in EulerFVMSolver.h).
        } else if (key == "flux_scheme") {
            if (value == "rusanov") flux_scheme = NumericalFluxScheme::Rusanov;
            else if (value == "hllc") flux_scheme = NumericalFluxScheme::HLLC;
            else if (value == "exact") flux_scheme = NumericalFluxScheme::Exact;
            else {
                std::cerr << "Error in case file '" << filename << "': unknown flux_scheme '"
                           << value << "' (expected 'rusanov', 'hllc' or 'exact')\n";
                return false;
            }
        } else if (key == "euler_init") {
            std::istringstream vs(value);
            std::string mode_word;
            vs >> mode_word;
            if (mode_word == "freestream") {
                euler_ic.mode = EulerICMode::Freestream;
                vs >> euler_ic.rho >> euler_ic.u >> euler_ic.v >> euler_ic.p;
            } else if (mode_word == "tworegion") {
                euler_ic.mode = EulerICMode::TwoRegion;
                vs >> euler_ic.rho_l >> euler_ic.u_l >> euler_ic.v_l >> euler_ic.p_l
                   >> euler_ic.rho_r >> euler_ic.u_r >> euler_ic.v_r >> euler_ic.p_r
                   >> euler_ic.x0;
            } else {
                std::cerr << "Error in case file '" << filename << "': unknown euler_init mode '"
                           << mode_word << "' (expected 'freestream' or 'tworegion')\n";
                return false;
            }
        // "ns_init = <mode> <values...>" -- same grammar as euler_init, own
        // storage (ns_ic), since a Navier-Stokes case file shouldn't need to
        // know that its initial condition happens to reuse Euler's format.
        } else if (key == "ns_init") {
            std::istringstream vs(value);
            std::string mode_word;
            vs >> mode_word;
            if (mode_word == "freestream") {
                ns_ic.mode = EulerICMode::Freestream;
                vs >> ns_ic.rho >> ns_ic.u >> ns_ic.v >> ns_ic.p;
            } else if (mode_word == "tworegion") {
                ns_ic.mode = EulerICMode::TwoRegion;
                vs >> ns_ic.rho_l >> ns_ic.u_l >> ns_ic.v_l >> ns_ic.p_l
                   >> ns_ic.rho_r >> ns_ic.u_r >> ns_ic.v_r >> ns_ic.p_r
                   >> ns_ic.x0;
            } else {
                std::cerr << "Error in case file '" << filename << "': unknown ns_init mode '"
                           << mode_word << "' (expected 'freestream' or 'tworegion')\n";
                return false;
            }
        // "rans_init = <mode> <values...>" -- same grammar as ns_init, own
        // storage (rans_ic).
        } else if (key == "rans_init") {
            std::istringstream vs(value);
            std::string mode_word;
            vs >> mode_word;
            if (mode_word == "freestream") {
                rans_ic.mode = EulerICMode::Freestream;
                vs >> rans_ic.rho >> rans_ic.u >> rans_ic.v >> rans_ic.p;
            } else if (mode_word == "tworegion") {
                rans_ic.mode = EulerICMode::TwoRegion;
                vs >> rans_ic.rho_l >> rans_ic.u_l >> rans_ic.v_l >> rans_ic.p_l
                   >> rans_ic.rho_r >> rans_ic.u_r >> rans_ic.v_r >> rans_ic.p_r
                   >> rans_ic.x0;
            } else {
                std::cerr << "Error in case file '" << filename << "': unknown rans_init mode '"
                           << mode_word << "' (expected 'freestream' or 'tworegion')\n";
                return false;
            }
        }
    }

    // Resolve reference_density/velocity_x/velocity_y/pressure's "auto"
    // default from the first ns_farfield (or, for equation == RANS,
    // rans_farfield) patch found, only when the wall diagnostics they feed
    // (wall_forces_file/wall_profile_file) are actually configured -- an
    // unrelated Navier-Stokes/RANS case file (e.g. one with no farfield patch
    // at all, like a Couette-flow setup) is never broken by these keys being
    // unset. See WallReferenceQuantities (WallTraction.h) and
    // docs/wall-diagnostics-plan.md.
    if ((equation == EquationSet::NavierStokes || equation == EquationSet::RANS) &&
        (!wall_forces_file.empty() || !wall_profile_file.empty())) {
        double farfield_rho = 0.0, farfield_u = 0.0, farfield_v = 0.0, farfield_p = 0.0;
        bool have_farfield = false;
        int farfield_count = 0;
        const char* farfield_keyword = (equation == EquationSet::RANS) ? "rans_farfield" : "ns_farfield";

        if (equation == EquationSet::RANS) {
            for (const auto& bc : rans_boundary_conditions) {
                if (bc.type == NSBoundaryType::Farfield) {
                    if (!have_farfield) {
                        farfield_rho = bc.rho; farfield_u = bc.u; farfield_v = bc.v; farfield_p = bc.p;
                        have_farfield = true;
                    }
                    ++farfield_count;
                }
            }
        } else {
            for (const auto& bc : ns_boundary_conditions) {
                if (bc.type == NSBoundaryType::Farfield) {
                    if (!have_farfield) {
                        farfield_rho = bc.rho; farfield_u = bc.u; farfield_v = bc.v; farfield_p = bc.p;
                        have_farfield = true;
                    }
                    ++farfield_count;
                }
            }
        }

        if (farfield_count > 1) {
            std::cerr << "Warning: multiple " << farfield_keyword << " patches found in '" << filename
                       << "'; using the first (in patch declaration order) as the freestream for auto "
                          "reference quantities\n";
        }
        if (!reference_density_set || !reference_velocity_x_set || !reference_velocity_y_set ||
            !reference_pressure_set) {
            if (!have_farfield) {
                std::cerr << "Error in case file '" << filename << "': reference_density/reference_velocity_x/"
                              "reference_velocity_y/reference_pressure must be set explicitly when no "
                           << farfield_keyword
                           << " patch exists to default them from (needed by wall_forces_file/wall_profile_file)\n";
                return false;
            }
            if (!reference_density_set) reference_density = farfield_rho;
            if (!reference_velocity_x_set) reference_velocity_x = farfield_u;
            if (!reference_velocity_y_set) reference_velocity_y = farfield_v;
            if (!reference_pressure_set) reference_pressure = farfield_p;
        }
    }

    return !mesh_file.empty() && !output_file.empty();
}
