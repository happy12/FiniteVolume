// SPDX-License-Identifier: GPL-3.0-only
#include "CaseInput.h"

#include <cmath>
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
//   true if the file was opened, both mesh_file and output_file were set,
//   every boundary/equation/euler_init keyword was recognized, every numeric
//   "key = value" parsed as a number, and any value with an obvious valid
//   domain (nsteps/dt/cfl/mu/gamma) was in range; false otherwise (a
//   descriptive message naming the file/line/key is printed to stderr).
bool CaseInput::load(const std::string& filename) {
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Error: could not open case file '" << filename << "'\n";
        return false;
    }

    // line_number is 1-based and only meaningful for the per-line errors
    // reported from inside the loop below; it is not incorporated into the
    // whole-file semantic checks (e.g. reference-quantity resolution) after
    // the loop, since those aren't tied to one line.
    int line_number = 0;

    // std::stod/std::stoi throw on a malformed value, which would otherwise
    // propagate uncaught out of load() and crash the program on a typo'd
    // case-file number. These wrap that in a line-numbered error instead.
    auto parse_double = [&](const std::string& key, const std::string& raw, double& out) -> bool {
        try {
            size_t pos = 0;
            double parsed = std::stod(raw, &pos);
            if (pos != raw.size()) throw std::invalid_argument("trailing characters");
            out = parsed;
            return true;
        } catch (const std::exception&) {
            std::cerr << "Error in case file '" << filename << "' at line " << line_number
                       << ": key '" << key << "' expects a number, got '" << raw << "'\n";
            return false;
        }
    };
    auto parse_int = [&](const std::string& key, const std::string& raw, int& out) -> bool {
        try {
            size_t pos = 0;
            int parsed = std::stoi(raw, &pos);
            if (pos != raw.size()) throw std::invalid_argument("trailing characters");
            out = parsed;
            return true;
        } catch (const std::exception&) {
            std::cerr << "Error in case file '" << filename << "' at line " << line_number
                       << ": key '" << key << "' expects an integer, got '" << raw << "'\n";
            return false;
        }
    };

    std::string line;
    while (std::getline(in, line)) {
        ++line_number;
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
            } else if (type_word == "ransSA_wall" || type_word == "ransSA_wall_isothermal" ||
                       type_word == "ransSA_wall_moving" || type_word == "ransSA_wall_moving_isothermal" ||
                       type_word == "ransSA_farfield" || type_word == "ransSA_outflow") {
                RANSBoundaryConditionSpecSA bc;
                bc.patch_name = patch_name;
                if (type_word == "ransSA_wall") {
                    bc.type = NSBoundaryType::NoSlipWall;
                } else if (type_word == "ransSA_wall_isothermal") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    bc.is_isothermal_wall = true;
                    iss >> bc.wall_temperature;
                } else if (type_word == "ransSA_wall_moving") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    iss >> bc.wall_u >> bc.wall_v;
                } else if (type_word == "ransSA_wall_moving_isothermal") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    bc.is_isothermal_wall = true;
                    iss >> bc.wall_u >> bc.wall_v >> bc.wall_temperature;
                } else if (type_word == "ransSA_outflow") {
                    bc.type = NSBoundaryType::Outflow;
                } else {
                    bc.type = NSBoundaryType::Farfield;
                    iss >> bc.rho >> bc.u >> bc.v >> bc.p >> bc.farfield_nut;
                }
                ransSA_boundary_conditions.push_back(bc);
            } else if (type_word == "ransSST_wall" || type_word == "ransSST_wall_isothermal" ||
                       type_word == "ransSST_wall_moving" || type_word == "ransSST_wall_moving_isothermal" ||
                       type_word == "ransSST_farfield" || type_word == "ransSST_outflow") {
                RANSBoundaryConditionSpecSST bc;
                bc.patch_name = patch_name;
                if (type_word == "ransSST_wall") {
                    bc.type = NSBoundaryType::NoSlipWall;
                } else if (type_word == "ransSST_wall_isothermal") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    bc.is_isothermal_wall = true;
                    iss >> bc.wall_temperature;
                } else if (type_word == "ransSST_wall_moving") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    iss >> bc.wall_u >> bc.wall_v;
                } else if (type_word == "ransSST_wall_moving_isothermal") {
                    bc.type = NSBoundaryType::NoSlipWall;
                    bc.is_isothermal_wall = true;
                    iss >> bc.wall_u >> bc.wall_v >> bc.wall_temperature;
                } else if (type_word == "ransSST_outflow") {
                    bc.type = NSBoundaryType::Outflow;
                } else {
                    bc.type = NSBoundaryType::Farfield;
                    iss >> bc.rho >> bc.u >> bc.v >> bc.p >> bc.farfield_k >> bc.farfield_omega;
                }
                ransSST_boundary_conditions.push_back(bc);
            } else {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": unknown boundary type '"
                           << type_word << "' for patch '" << patch_name
                           << "' (expected 'dirichlet', 'neumann', 'wall', 'farfield', 'outflow', 'ns_wall', "
                           << "'ns_wall_isothermal', 'ns_wall_moving', 'ns_wall_moving_isothermal', "
                           << "'ns_farfield', 'ns_outflow', 'ransSA_wall', 'ransSA_wall_isothermal', "
                           << "'ransSA_wall_moving', 'ransSA_wall_moving_isothermal', 'ransSA_farfield', "
                           << "'ransSA_outflow', 'ransSST_wall', 'ransSST_wall_isothermal', 'ransSST_wall_moving', "
                           << "'ransSST_wall_moving_isothermal', 'ransSST_farfield' or 'ransSST_outflow')\n";
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
        else if (key == "alpha") {                                    // Diffusion coefficient, in (mesh length units)^2 / (time units)
            if (!parse_double(key, value, alpha)) return false;
        }
        else if (key == "dt") {                                       // Time step size, in time units (must satisfy the explicit-scheme stability limit)
            if (!parse_double(key, value, dt)) return false;
            if (dt <= 0.0) {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": dt must be positive (got " << dt << ")\n";
                return false;
            }
        }
        else if (key == "nsteps") {                                   // Number of time steps to advance
            if (!parse_int(key, value, nsteps)) return false;
            if (nsteps <= 0) {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": nsteps must be positive (got " << nsteps << ")\n";
                return false;
            }
        }
        else if (key == "initial_value") { if (!parse_double(key, value, initial_value)) return false; }   // phi value inside the initial condition disc
        else if (key == "initial_radius") { if (!parse_double(key, value, initial_radius)) return false; } // Radius (from the origin) of the initial condition disc, in mesh length units
        else if (key == "u_adv") { if (!parse_double(key, value, u_adv)) return false; }            // Uniform advection velocity x-component (advection-diffusion)
        else if (key == "v_adv") { if (!parse_double(key, value, v_adv)) return false; }            // Uniform advection velocity y-component (advection-diffusion)
        else if (key == "gamma") {                                    // Ratio of specific heats, dimensionless (Euler / Navier-Stokes)
            if (!parse_double(key, value, gamma)) return false;
            if (gamma <= 1.0) {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": gamma must be greater than 1 (got " << gamma << ")\n";
                return false;
            }
        }
        else if (key == "cfl") {                                      // CFL number, dimensionless, controlling the adaptive time step (Euler / Navier-Stokes)
            if (!parse_double(key, value, cfl)) return false;
            if (cfl <= 0.0) {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": cfl must be positive (got " << cfl << ")\n";
                return false;
            }
        }
        else if (key == "cfl_min") { if (!parse_double(key, value, cfl_ramp.cfl_min)) return false; }                              // ramp only: starting/floor CFL
        else if (key == "cfl_ramp_growth") { if (!parse_double(key, value, cfl_ramp.growth_factor)) return false; }                // ramp only: multiplier when residual trend is decreasing
        else if (key == "cfl_ramp_shrink") { if (!parse_double(key, value, cfl_ramp.shrink_factor)) return false; }                // ramp only: multiplier when residual trend is mildly rising
        else if (key == "cfl_ramp_window") { if (!parse_int(key, value, cfl_ramp.window)) return false; }                          // ramp only: steps of residual history the trend is computed over
        else if (key == "cfl_ramp_divergence_threshold") { if (!parse_double(key, value, cfl_ramp.divergence_threshold)) return false; } // ramp only: log10(residual) rise over the window that forces a hard reset to cfl_min
        else if (key == "mu") {                                       // Dynamic viscosity, mesh-consistent units (Navier-Stokes)
            if (!parse_double(key, value, mu)) return false;
            if (mu < 0.0) {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": mu must be non-negative (got " << mu << ")\n";
                return false;
            }
        }
        else if (key == "prandtl") { if (!parse_double(key, value, prandtl)) return false; }        // Prandtl number, dimensionless (Navier-Stokes)
        else if (key == "gas_constant") { if (!parse_double(key, value, gas_constant)) return false; } // Specific gas constant R in p = rho*R*T (Navier-Stokes / RANS)
        else if (key == "prandtl_t") { if (!parse_double(key, value, prandtl_t)) return false; }        // Turbulent Prandtl number, dimensionless (RANS)
        else if (key == "initial_nut") { if (!parse_double(key, value, initial_nut)) return false; }    // Uniform initial nut applied to every cell (RANS)
        else if (key == "sa_cb1") { if (!parse_double(key, value, sa_constants.cb1)) return false; }     // SA-noft2 model constants (RANS) -- see SpalartAllmaras.h
        else if (key == "sa_cb2") { if (!parse_double(key, value, sa_constants.cb2)) return false; }
        else if (key == "sa_sigma") { if (!parse_double(key, value, sa_constants.sigma)) return false; }
        else if (key == "sa_kappa") { if (!parse_double(key, value, sa_constants.kappa)) return false; }
        else if (key == "sa_cw2") { if (!parse_double(key, value, sa_constants.cw2)) return false; }
        else if (key == "sa_cw3") { if (!parse_double(key, value, sa_constants.cw3)) return false; }
        else if (key == "sa_cv1") { if (!parse_double(key, value, sa_constants.cv1)) return false; }
        else if (key == "sa_cv2") { if (!parse_double(key, value, sa_constants.cv2)) return false; }
        else if (key == "sa_cv3") { if (!parse_double(key, value, sa_constants.cv3)) return false; }
        else if (key == "initial_k") { if (!parse_double(key, value, initial_k)) return false; initial_k_set = true; }         // Uniform initial k applied to every cell (RANS k-omega SST)
        else if (key == "initial_omega") { if (!parse_double(key, value, initial_omega)) return false; initial_omega_set = true; } // Uniform initial omega applied to every cell (RANS k-omega SST)
        else if (key == "sst_turbulence_intensity") { if (!parse_double(key, value, sst_turbulence_intensity)) return false; } // Tu, dimensionless; derives initial_k/initial_omega if unset (RANS k-omega SST)
        else if (key == "sst_eddy_viscosity_ratio") { if (!parse_double(key, value, sst_eddy_viscosity_ratio)) return false; }  // nu_t/nu, dimensionless; derives initial_k/initial_omega if unset (RANS k-omega SST)
        else if (key == "sst_kato_launder") sst_kato_launder = (value == "true" || value == "1"); // Kato-Launder production limiter (RANS k-omega SST) -- see SSTKOmega.h
        else if (key == "sst_beta_star") { if (!parse_double(key, value, sst_constants.beta_star)) return false; }  // SST model constants (RANS k-omega SST) -- see SSTKOmega.h
        else if (key == "sst_kappa") { if (!parse_double(key, value, sst_constants.kappa)) return false; }
        else if (key == "sst_a1") { if (!parse_double(key, value, sst_constants.a1)) return false; }
        else if (key == "sst_sigma_k1") { if (!parse_double(key, value, sst_constants.sigma_k1)) return false; }
        else if (key == "sst_sigma_k2") { if (!parse_double(key, value, sst_constants.sigma_k2)) return false; }
        else if (key == "sst_sigma_omega1") { if (!parse_double(key, value, sst_constants.sigma_omega1)) return false; }
        else if (key == "sst_sigma_omega2") { if (!parse_double(key, value, sst_constants.sigma_omega2)) return false; }
        else if (key == "sst_beta1") { if (!parse_double(key, value, sst_constants.beta1)) return false; }
        else if (key == "sst_beta2") { if (!parse_double(key, value, sst_constants.beta2)) return false; }
        else if (key == "resolution_report_file") resolution_report_file = value; // Path to resolution-diagnostic CSV; empty = disabled (Navier-Stokes)
        else if (key == "resolution_report_interval") { if (!parse_int(key, value, resolution_report_interval)) return false; } // Write a row every N steps
        else if (key == "wall_forces_file") wall_forces_file = value; // Path to per-(step, patch) force/moment CSV; empty = disabled (Navier-Stokes)
        else if (key == "wall_forces_interval") { if (!parse_int(key, value, wall_forces_interval)) return false; } // Write a row block every N steps
        else if (key == "wall_profile_file") wall_profile_file = value; // Path to per-wall-node Cf/Cp/y+/BL-thickness CSV; empty = disabled (Navier-Stokes)
        else if (key == "wall_profile_interval") { if (!parse_int(key, value, wall_profile_interval)) return false; } // 0 = write once at the run's natural end
        else if (key == "reference_density") { if (!parse_double(key, value, reference_density)) return false; reference_density_set = true; }
        else if (key == "reference_velocity_x") { if (!parse_double(key, value, reference_velocity_x)) return false; reference_velocity_x_set = true; }
        else if (key == "reference_velocity_y") { if (!parse_double(key, value, reference_velocity_y)) return false; reference_velocity_y_set = true; }
        else if (key == "reference_pressure") { if (!parse_double(key, value, reference_pressure)) return false; reference_pressure_set = true; }
        else if (key == "reference_length") { if (!parse_double(key, value, reference_length)) return false; } // Length scale for Cd/Cl/Cm; 1.0 = per-unit-span
        else if (key == "moment_reference_x") { if (!parse_double(key, value, moment_reference_x)) return false; } // Point Cm is taken about
        else if (key == "moment_reference_y") { if (!parse_double(key, value, moment_reference_y)) return false; }
        else if (key == "boundary_layer_max_distance") { if (!parse_double(key, value, boundary_layer_max_distance)) return false; } // point-location only; <= 0 = auto
        else if (key == "boundary_layer_n_samples") { if (!parse_int(key, value, boundary_layer_n_samples)) return false; }        // point-location only
        else if (key == "residual_file") residual_file = value;              // Path to the residual history CSV; empty = disabled
        else if (key == "residual_interval") { if (!parse_int(key, value, residual_interval)) return false; } // Write a residual row every N steps
        else if (key == "write_interval") { if (!parse_int(key, value, write_interval)) return false; } // Write a numbered VTK snapshot every N steps; 0 = disabled
        else if (key == "num_threads") { if (!parse_int(key, value, num_threads)) return false; }       // OpenMP thread count, dimensionless; 0 = OpenMP's own default
        else if (key == "residual_tolerance") { if (!parse_double(key, value, residual_tolerance)) return false; }           // diffusion only; < 0 = disabled
        else if (key == "residual_tolerance_rho") { if (!parse_double(key, value, residual_tolerance_rho)) return false; }   // Euler only; < 0 = disabled
        else if (key == "residual_tolerance_rho_u") { if (!parse_double(key, value, residual_tolerance_rho_u)) return false; }
        else if (key == "residual_tolerance_rho_v") { if (!parse_double(key, value, residual_tolerance_rho_v)) return false; }
        else if (key == "residual_tolerance_E") { if (!parse_double(key, value, residual_tolerance_E)) return false; }
        else if (key == "residual_tolerance_nut") { if (!parse_double(key, value, residual_tolerance_nut)) return false; } // RANS (Spalart-Allmaras) only; < 0 = disabled
        else if (key == "residual_tolerance_k") { if (!parse_double(key, value, residual_tolerance_k)) return false; } // RANS (k-omega SST) only; < 0 = disabled
        else if (key == "residual_tolerance_omega") { if (!parse_double(key, value, residual_tolerance_omega)) return false; } // RANS (k-omega SST) only; < 0 = disabled
        else if (key == "checkpoint_file") checkpoint_file = value;          // Path to save/resume solver state; empty = disabled
        // Exact Riemann solver's Newton-Raphson tolerance/iteration cap (see
        // ExactRiemannFlux.h); exposed for transparency, not meant to be
        // tuned in normal use.
        else if (key == "exact_riemann_tol") { if (!parse_double(key, value, exact_riemann_tol)) return false; }
        else if (key == "exact_riemann_max_iter") { if (!parse_int(key, value, exact_riemann_max_iter)) return false; }
        else if (key == "scratch_dir") scratch_dir = value;                  // Base dir for relative output_file/checkpoint_file/residual_file; empty = disabled
        // "output_precision = <1-17>" -- significant digits written for every
        // output double (VTK results/snapshots, residual_file CSV, FvMeshWriter).
        else if (key == "output_precision") {
            if (!parse_int(key, value, output_precision)) return false;
            if (output_precision < 1 || output_precision > 17) {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": output_precision must be between 1 and 17 (got "
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
            else if (value == "rans_sa") equation = EquationSet::RANS_SA;
            else if (value == "rans_sst") equation = EquationSet::RANS_SST;
            else {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": unknown equation '"
                           << value << "' (expected 'diffusion', 'euler', 'advection_diffusion', "
                              "'navier_stokes', 'rans_sa' or 'rans_sst')\n";
                return false;
            }
        // "cfl_mode = fixed|ramp" -- selects CflMode (see CflRamp.h), used by
        // the four cfl-driven equation sets (euler/navier_stokes/rans_sa/rans_sst).
        } else if (key == "cfl_mode") {
            if (value == "fixed") cfl_mode = CflMode::Fixed;
            else if (value == "ramp") cfl_mode = CflMode::Ramp;
            else {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": unknown cfl_mode '"
                           << value << "' (expected 'fixed' or 'ramp')\n";
                return false;
            }
        // "sst_limiter = vorticity|strain_rate" -- selects SSTLimiterVariant
        // (see SSTKOmega.h), used by the RANS (k-omega SST) equation set.
        } else if (key == "sst_limiter") {
            if (value == "vorticity") sst_limiter_variant = SSTLimiterVariant::Vorticity;
            else if (value == "strain_rate") sst_limiter_variant = SSTLimiterVariant::StrainRate;
            else {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": unknown sst_limiter '"
                           << value << "' (expected 'vorticity' or 'strain_rate')\n";
                return false;
            }
        // "gradient_scheme = least-squares|green-gauss" -- selects
        // GradientCalculator's scheme (see GradientReconstruction.h),
        // currently used by the advection-diffusion equation set.
        } else if (key == "gradient_scheme") {
            if (value == "least-squares") gradient_scheme = GradientScheme::LeastSquares;
            else if (value == "green-gauss") gradient_scheme = GradientScheme::GreenGauss;
            else {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": unknown gradient_scheme '"
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
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": unknown boundary_layer_method '"
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
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": unknown flux_scheme '"
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
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": unknown euler_init mode '"
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
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": unknown ns_init mode '"
                           << mode_word << "' (expected 'freestream' or 'tworegion')\n";
                return false;
            }
        // "ransSA_init = <mode> <values...>" -- same grammar as ns_init, own
        // storage (ransSA_ic).
        } else if (key == "ransSA_init") {
            std::istringstream vs(value);
            std::string mode_word;
            vs >> mode_word;
            if (mode_word == "freestream") {
                ransSA_ic.mode = EulerICMode::Freestream;
                vs >> ransSA_ic.rho >> ransSA_ic.u >> ransSA_ic.v >> ransSA_ic.p;
            } else if (mode_word == "tworegion") {
                ransSA_ic.mode = EulerICMode::TwoRegion;
                vs >> ransSA_ic.rho_l >> ransSA_ic.u_l >> ransSA_ic.v_l >> ransSA_ic.p_l
                   >> ransSA_ic.rho_r >> ransSA_ic.u_r >> ransSA_ic.v_r >> ransSA_ic.p_r
                   >> ransSA_ic.x0;
            } else {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": unknown ransSA_init mode '"
                           << mode_word << "' (expected 'freestream' or 'tworegion')\n";
                return false;
            }
        // "ransSST_init = <mode> <values...>" -- same grammar as ransSA_init,
        // own storage (ransSST_ic).
        } else if (key == "ransSST_init") {
            std::istringstream vs(value);
            std::string mode_word;
            vs >> mode_word;
            if (mode_word == "freestream") {
                ransSST_ic.mode = EulerICMode::Freestream;
                vs >> ransSST_ic.rho >> ransSST_ic.u >> ransSST_ic.v >> ransSST_ic.p;
            } else if (mode_word == "tworegion") {
                ransSST_ic.mode = EulerICMode::TwoRegion;
                vs >> ransSST_ic.rho_l >> ransSST_ic.u_l >> ransSST_ic.v_l >> ransSST_ic.p_l
                   >> ransSST_ic.rho_r >> ransSST_ic.u_r >> ransSST_ic.v_r >> ransSST_ic.p_r
                   >> ransSST_ic.x0;
            } else {
                std::cerr << "Error in case file '" << filename << "' at line " << line_number
                           << ": unknown ransSST_init mode '"
                           << mode_word << "' (expected 'freestream' or 'tworegion')\n";
                return false;
            }
        }
    }

    // Resolve reference_density/velocity_x/velocity_y/pressure's "auto"
    // default from the first ns_farfield (or, for equation == RANS_SA,
    // ransSA_farfield; equation == RANS_SST, ransSST_farfield) patch found,
    // only when the wall diagnostics they feed (wall_forces_file/
    // wall_profile_file) are actually configured -- an unrelated
    // Navier-Stokes/RANS case file (e.g. one with no farfield patch at all,
    // like a Couette-flow setup) is never broken by these keys being unset.
    // See WallReferenceQuantities (WallTraction.h) and
    // docs/wall-diagnostics-plan.md.
    if ((equation == EquationSet::NavierStokes || equation == EquationSet::RANS_SA ||
         equation == EquationSet::RANS_SST) &&
        (!wall_forces_file.empty() || !wall_profile_file.empty())) {
        double farfield_rho = 0.0, farfield_u = 0.0, farfield_v = 0.0, farfield_p = 0.0;
        bool have_farfield = false;
        int farfield_count = 0;
        const char* farfield_keyword = (equation == EquationSet::RANS_SA)
                                            ? "ransSA_farfield"
                                            : (equation == EquationSet::RANS_SST) ? "ransSST_farfield" : "ns_farfield";

        if (equation == EquationSet::RANS_SA) {
            for (const auto& bc : ransSA_boundary_conditions) {
                if (bc.type == NSBoundaryType::Farfield) {
                    if (!have_farfield) {
                        farfield_rho = bc.rho; farfield_u = bc.u; farfield_v = bc.v; farfield_p = bc.p;
                        have_farfield = true;
                    }
                    ++farfield_count;
                }
            }
        } else if (equation == EquationSet::RANS_SST) {
            for (const auto& bc : ransSST_boundary_conditions) {
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

    // Derive initial_k/initial_omega from sst_turbulence_intensity (Tu) and
    // sst_eddy_viscosity_ratio, per CaseInput.h's class comment -- ONLY when
    // the case file didn't set either explicitly, and only from a
    // ransSST_init freestream (a TwoRegion initial condition has no single
    // velocity magnitude to derive Tu-based k/omega from, so this
    // derivation is skipped there; initial_k/initial_omega then stay at
    // their explicit-or-default values, same as if sst_turbulence_intensity/
    // sst_eddy_viscosity_ratio had been omitted entirely).
    if (equation == EquationSet::RANS_SST && ransSST_ic.mode == EulerICMode::Freestream &&
        sst_turbulence_intensity > 0.0 && sst_eddy_viscosity_ratio > 0.0) {
        double speed = std::sqrt(ransSST_ic.u * ransSST_ic.u + ransSST_ic.v * ransSST_ic.v);
        double derived_k = 1.5 * (sst_turbulence_intensity * speed) * (sst_turbulence_intensity * speed);
        double nu_freestream = (ransSST_ic.rho > 0.0) ? mu / ransSST_ic.rho : 0.0;
        double derived_omega = (nu_freestream > 0.0) ? derived_k / (sst_eddy_viscosity_ratio * nu_freestream) : 0.0;
        if (!initial_k_set) initial_k = derived_k;
        if (!initial_omega_set) initial_omega = derived_omega;
    }

    if (mesh_file.empty()) {
        std::cerr << "Error in case file '" << filename << "': required key 'mesh_file' was not set\n";
        return false;
    }
    if (output_file.empty()) {
        std::cerr << "Error in case file '" << filename << "': required key 'output_file' was not set\n";
        return false;
    }
    return true;
}
