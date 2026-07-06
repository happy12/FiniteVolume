// SPDX-License-Identifier: GPL-3.0-only
#include <algorithm>
#include <cmath>
#include <csignal>
#include <functional>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include "UnstructuredFVMSolver.h"
#include "EulerFVMSolver.h"
#include "AdvectionDiffusionFVMSolver.h"
#include "NavierStokesFVMSolver.h"
#include "MeshReader.h"
#include "CaseInput.h"
#include "Checkpoint.h"
#include "VtkWriter.h"
#include "GradientReconstruction.h"
#include "WallDistance.h"
#include "WallTraction.h"
#include "SpalartAllmaras.h"
#include "RANSFVMSolver.h"
#include "Version.h"

namespace {

// Set by handle_sigint() when the user presses Ctrl+C, checked once per loop
// iteration in run_diffusion()/run_euler() so a requested interrupt breaks
// the loop the same way any other "natural" stop does (sharing the same
// checkpoint-write + final-output code path). sig_atomic_t is the only type
// the C++ standard guarantees is safe to write from a signal handler; the
// handler does nothing else (no I/O, no allocation).
volatile sig_atomic_t g_interrupt_requested = 0;
void handle_sigint(int) { g_interrupt_requested = 1; }

// Resolves a case-file output path (output_file/checkpoint_file/
// residual_file) against scratch_dir, so relative paths can be redirected to
// one place without editing every key individually. Absolute paths and an
// unset scratch_dir/path are passed through unchanged; mesh_file is never
// resolved this way (it's an input, not output).
// Input:  path        - the path as read from the case file (may be empty,
//                        meaning that output is disabled -- left as-is)
//         scratch_dir - case_input.scratch_dir; empty means disabled
// Returns: 'path' unchanged, or 'path' rebased under scratch_dir
std::string resolve_output_path(const std::string& path, const std::string& scratch_dir) {
    if (path.empty() || scratch_dir.empty()) {
        return path;
    }
    std::filesystem::path p(path);
    if (p.is_absolute()) {
        return path;
    }
    return (std::filesystem::path(scratch_dir) / p).string();
}

// Creates 'path's parent directory (recursively) if it doesn't already
// exist, so a long-running solve doesn't fail at its very last step just
// because a destination folder was never created. Input: path - a file path
// (may be relative/absolute; empty means "disabled", skipped; no directory
// component means the working directory, already guaranteed to exist).
// Output: the directory is created if missing. Failures are ignored here --
// the subsequent ofstream::open() at the actual write site will fail (and
// report) the same way it always has if creation didn't succeed.
void ensure_parent_directory(const std::string& path) {
    if (path.empty()) {
        return;
    }
    std::filesystem::path parent = std::filesystem::path(path).parent_path();
    if (!parent.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(parent, ec);
    }
}

// Builds a per-step snapshot filename by inserting a zero-padded step number
// before output_file's extension (or at the end, if it has none).
// Input:  output_file - e.g. "output/result.vtk"
//         step_index  - current step number
// Returns: e.g. "output/result_000100.vtk" for step_index == 100
std::string numbered_filename(const std::string& output_file, int step_index) {
    size_t slash = output_file.find_last_of("/\\");
    size_t dot = output_file.find_last_of('.');
    // Only treat a '.' as the extension separator if it comes after the last
    // path separator (so a directory name containing '.' isn't mistaken for one).
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) {
        dot = output_file.size();
    }

    std::ostringstream oss;
    oss << output_file.substr(0, dot) << "_" << std::setw(6) << std::setfill('0') << step_index
        << output_file.substr(dot);
    return oss.str();
}

// Returns the diagonal length of 'mesh's node bounding box -- used as an
// automatic default for boundary_layer_max_distance (CaseInput.h) when the
// case file doesn't set one: comfortably larger than any boundary-layer
// thickness the mesh could contain, and compute_boundary_layer_profiles_point_location()
// stops as soon as a sample point falls outside the mesh domain, so an
// over-generous distance costs a few extra out-of-domain samples, not
// meaningfully more work.
double mesh_bounding_box_diagonal(const UnstructuredMesh& mesh) {
    double min_x = mesh.nodes[0].x, max_x = mesh.nodes[0].x;
    double min_y = mesh.nodes[0].y, max_y = mesh.nodes[0].y;
    for (const Node& n : mesh.nodes) {
        min_x = std::min(min_x, n.x); max_x = std::max(max_x, n.x);
        min_y = std::min(min_y, n.y); max_y = std::max(max_y, n.y);
    }
    return std::hypot(max_x - min_x, max_y - min_y);
}

// Collects the indices of every boundary face whose patch's NSBoundaryType is
// NoSlipWall, for WallTraction.h's compute_wall_traction_samples() /
// compute_boundary_layer_profiles() -- identical purpose to
// RANSFVMSolver.cpp's own private collect_wall_faces(), duplicated here since
// that one isn't reachable from main.cpp.
std::vector<int> collect_wall_faces(const UnstructuredMesh& mesh, const std::vector<NSBoundaryCondition>& bcs) {
    std::vector<int> wall_faces;
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        if (face.cell_right == -1 && bcs.at(face.patch_id).type == NSBoundaryType::NoSlipWall) {
            wall_faces.push_back((int)i);
        }
    }
    return wall_faces;
}

// RANS overload: identical purpose, indexed by RANSBoundaryCondition::ns.type instead.
std::vector<int> collect_wall_faces(const UnstructuredMesh& mesh, const std::vector<RANSBoundaryCondition>& bcs) {
    std::vector<int> wall_faces;
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        if (face.cell_right == -1 && bcs.at(face.patch_id).ns.type == NSBoundaryType::NoSlipWall) {
            wall_faces.push_back((int)i);
        }
    }
    return wall_faces;
}

// Builds WallReferenceQuantities (WallTraction.h) from a CaseInput whose
// reference_* fields have already been resolved (auto-fallback from the
// first ns_farfield patch, or a load()-time error if that wasn't possible --
// see CaseInput::load()).
WallReferenceQuantities build_wall_reference_quantities(const CaseInput& case_input) {
    WallReferenceQuantities ref;
    ref.rho_ref = case_input.reference_density;
    ref.velocity_x_ref = case_input.reference_velocity_x;
    ref.velocity_y_ref = case_input.reference_velocity_y;
    ref.p_ref = case_input.reference_pressure;
    ref.length_ref = case_input.reference_length;
    ref.moment_reference_x = case_input.moment_reference_x;
    ref.moment_reference_y = case_input.moment_reference_y;
    return ref;
}

// Writes one wall_forces_file row per WallForceReport (one per NoSlipWall
// patch, plus a final "TOTAL" domain-total row -- see compute_wall_forces()).
void write_wall_forces_rows(std::ofstream& out, long long step, const UnstructuredMesh& mesh,
                              const std::vector<WallForceReport>& reports) {
    for (const auto& r : reports) {
        std::string patch_name = (r.patch_id == -1) ? "TOTAL" : mesh.patches[r.patch_id].name;
        out << step << "," << patch_name << "," << r.friction_drag << "," << r.pressure_drag << ","
            << r.total_drag << "," << r.cd_friction << "," << r.cd_pressure << "," << r.cd_total << ","
            << r.lift << "," << r.cl << "," << r.moment << "," << r.cm << "\n";
    }
}

// Writes a wall_profile_file snapshot (opened fresh each call, unlike
// wall_forces_file's append-forever stream, since this is a per-location
// distribution at a single instant, not a time history).
bool write_wall_profile_snapshot(const std::string& path, int precision, const UnstructuredMesh& mesh,
                                   const std::vector<WallNodeSample>& nodes) {
    std::ofstream out(path);
    if (!out.is_open()) return false;
    out << std::setprecision(precision);
    out << "x,y,patch_name,tau_wall,Cf,Cp,y_plus,delta_99,displacement_thickness,momentum_thickness,shape_factor\n";
    for (const auto& n : nodes) {
        out << n.x << "," << n.y << "," << mesh.patches[n.patch_id].name << "," << n.tau_wall << "," << n.Cf << ","
            << n.Cp << "," << n.y_plus << "," << n.delta_99 << "," << n.displacement_thickness << ","
            << n.momentum_thickness << "," << n.shape_factor << "\n";
    }
    return true;
}

// Builds a deterministic n_x x n_y structured quad mesh on the unit square,
// mapping each reference-grid node (X, Y) = (i/n_x, j/n_y) through
// 'displace' to get its physical position. Shared topology-building code
// for build_skewed_verification_mesh() and build_sheared_verification_mesh()
// below, which differ only in that mapping (and, for the latter, in wanting
// n_x != n_y -- see its own comment). Its 4 sides are tagged as named
// boundary patches ("bottom"/"left"/"top"/"right") so a caller can assign
// different BCs per side; callers that don't care about per-side BCs (e.g.
// run_verify_gradient(), which only cares whether a face is a boundary face
// at all, via cell_right == -1) can ignore them.
//
// Faces are built directly from the known structured connectivity (no
// generic edge-matching needed, since the grid topology is known up front);
// cell volumes/centroids are computed via the standard polygon shoelace
// formula, matching how MeshReader computes them at parse time.
// mesh.compute_geometry() must still be called afterwards to fill in face
// midpoints/areas/normals.
//
// Input:  n_x, n_y - number of cells along each reference axis (n_x*n_y cells total)
//         displace - maps a reference-grid node (X, Y), X in [0,1], Y in
//                     [0,1], to its physical (x, y) position; must keep
//                     every cell a simple (non-self-intersecting) polygon
// Output: mesh is replaced with the generated nodes/cells/faces/patches
void build_structured_verification_mesh(UnstructuredMesh& mesh, int n_x, int n_y,
                                          const std::function<Node(double, double)>& displace) {
    double hx = 1.0 / n_x, hy = 1.0 / n_y;
    auto node_index = [n_x](int i, int j) { return j * (n_x + 1) + i; };

    mesh.patches = {{"bottom"}, {"left"}, {"top"}, {"right"}};
    const int PATCH_BOTTOM = 0, PATCH_LEFT = 1, PATCH_TOP = 2, PATCH_RIGHT = 3;

    mesh.nodes.resize((n_x + 1) * (n_y + 1));
    for (int j = 0; j <= n_y; ++j) {
        for (int i = 0; i <= n_x; ++i) {
            mesh.nodes[node_index(i, j)] = displace(i * hx, j * hy);
        }
    }

    mesh.cells.resize((size_t)n_x * n_y);
    for (int j = 0; j < n_y; ++j) {
        for (int i = 0; i < n_x; ++i) {
            Cell& cell = mesh.cells[j * n_x + i];
            cell.node_ids = {node_index(i, j), node_index(i + 1, j),
                              node_index(i + 1, j + 1), node_index(i, j + 1)};

            // Polygon shoelace formula for signed area and centroid (CCW
            // winding, as constructed above, gives a positive signed area).
            double a2 = 0.0, cx = 0.0, cy = 0.0;
            for (size_t k = 0; k < cell.node_ids.size(); ++k) {
                const Node& p1 = mesh.nodes[cell.node_ids[k]];
                const Node& p2 = mesh.nodes[cell.node_ids[(k + 1) % cell.node_ids.size()]];
                double cross = p1.x * p2.y - p2.x * p1.y;
                a2 += cross;
                cx += (p1.x + p2.x) * cross;
                cy += (p1.y + p2.y) * cross;
            }
            cell.volume = std::abs(a2) * 0.5;
            cell.x_centroid = cx / (3.0 * a2);
            cell.y_centroid = cy / (3.0 * a2);
        }
    }

    // Each cell contributes its own "bottom" and "left" edge (shared with
    // the cell below/to the left, or a domain boundary if there is none),
    // plus a "top"/"right" edge only for the top row/right column, where no
    // cell above/to the right exists to have already claimed that edge as
    // its own bottom/left. Node order follows each cell's CCW node_ids
    // traversal, so UnstructuredMesh::compute_geometry()'s outward-normal
    // convention comes out correct for free.
    auto add_face = [&mesh](int n1, int n2, int cell_left, int cell_right, int patch_id) {
        Face face;
        face.node1 = n1;
        face.node2 = n2;
        face.cell_left = cell_left;
        face.cell_right = cell_right;
        face.patch_id = (cell_right == -1) ? patch_id : -1;
        mesh.faces.push_back(face);
        int face_idx = (int)mesh.faces.size() - 1;
        mesh.cells[cell_left].faces.push_back(face_idx);
        if (cell_right != -1) mesh.cells[cell_right].faces.push_back(face_idx);
    };

    for (int j = 0; j < n_y; ++j) {
        for (int i = 0; i < n_x; ++i) {
            int c = j * n_x + i;
            add_face(node_index(i, j), node_index(i + 1, j), c, j > 0 ? c - n_x : -1, PATCH_BOTTOM);
            add_face(node_index(i, j + 1), node_index(i, j), c, i > 0 ? c - 1 : -1, PATCH_LEFT);
            if (j == n_y - 1) add_face(node_index(i + 1, j + 1), node_index(i, j + 1), c, -1, PATCH_TOP);
            if (i == n_x - 1) add_face(node_index(i + 1, j), node_index(i + 1, j + 1), c, -1, PATCH_RIGHT);
        }
    }
}

// A smooth, non-affine "wavy" perturbation of a Cartesian grid (the outer
// boundary is left exactly straight, since the perturbation vanishes at
// X/Y = 0 or 1), used as a self-contained non-orthogonal test mesh by
// run_verify_gradient() and run_verify_advection_diffusion().
//
// A merely SHEARED (affine) Cartesian grid would not actually stress a cell
// gradient scheme: an affine map preserves collinearity, so if every face
// midpoint lies on the line joining its two cell centroids in the
// (undistorted) reference grid -- true for any Cartesian mesh -- it still
// does after the affine map, and Green-Gauss stays exact "by accident". This
// smooth but non-affine wobble breaks that collinearity for real, which is
// what actually exercises Green-Gauss's known non-orthogonality error (see
// GradientReconstruction.h) and gives Phase 1's face-gradient correction
// something genuine to fix. It is NOT used for run_verify_couette() (see
// build_sheared_verification_mesh()) precisely because it varies with x:
// Couette flow's analytic solution assumes x-translation invariance, and a
// mesh whose own distortion varies with x breaks that invariance in the
// discretization itself, contaminating the comparison with a spurious
// x-varying truncation-error "forcing" (found the hard way -- see
// docs/navier-stokes-tracker.md Phase 4).
//
// Input:  n      - number of cells per side (n*n cells total)
//         wobble - perturbation amplitude, in mesh length units (0 = exact
//                  Cartesian grid; must stay small enough relative to 1/n
//                  that cells remain convex/non-degenerate)
// Output: mesh is replaced with the generated nodes/cells/faces/patches
void build_skewed_verification_mesh(UnstructuredMesh& mesh, int n, double wobble) {
    const double pi = 3.14159265358979323846;
    build_structured_verification_mesh(mesh, n, n, [wobble, pi](double X, double Y) {
        double perturb = wobble * std::sin(pi * X) * std::sin(pi * Y);
        return Node{X + perturb, Y + perturb};
    });
}

// A pure shear of a Cartesian grid (x = X + shear*Y, y = Y), used as a
// non-orthogonal test mesh specifically where x-translation invariance must
// be preserved (run_verify_couette()): a shear's Jacobian is exactly 1
// everywhere regardless of 'shear' (so cells never become degenerate/
// inverted, unlike an unbounded wobble amplitude), and every cell at a given
// row is an identical parallelogram just translated in x, so the
// discretization doesn't distinguish between x-positions the way
// build_skewed_verification_mesh()'s wavy perturbation does. It is still
// genuinely non-orthogonal in Jasak's sense: a horizontal face's normal
// stays exactly vertical after shearing, but the line joining the two
// cells' centroids picks up a horizontal offset (since shear displaces a
// cell's centroid by an amount that depends on the cell's own y-position),
// so the face normal and the centroid-connecting line are no longer
// parallel -- exactly the condition face_gradient()'s over-relaxed
// correction (Phase 1) exists to handle.
//
// n_x defaults to 1 in run_verify_couette() for a reason beyond "cheaper":
// NSBoundaryType::Outflow is a zero-gradient/do-nothing condition, not a
// true periodic one, so it does NOT actually force the two ends of an x
// column to match. Nothing in that BC prevents (or damps) a slow, smoothly
// linear-in-x drift mode from developing across an x-extended domain --
// found the hard way (see docs/navier-stokes-tracker.md Phase 4) as a
// persistent, non-decaying secondary flow that survived 10x more time
// steps unchanged. A single cell of x-extent removes the problem
// structurally: there is no room for an x-varying mode to exist at all,
// while the shear still makes every internal y-face non-orthogonal.
//
// Input:  n_x, n_y - number of cells along each axis (n_x*n_y cells total)
//         shear    - horizontal shift per unit y (0 = exact Cartesian grid)
// Output: mesh is replaced with the generated nodes/cells/faces/patches
void build_sheared_verification_mesh(UnstructuredMesh& mesh, int n_x, int n_y, double shear) {
    build_structured_verification_mesh(mesh, n_x, n_y, [shear](double X, double Y) {
        return Node{X + shear * Y, Y};
    });
}

// Builds a flat-plate boundary-layer mesh: n_x uniform streamwise cells over
// [0, L], wall-normal cells over [0, H] geometrically stretched from the
// wall (first_cell_height, then each successive cell first_cell_height
// times growth_ratio taller) -- needed to put the first cell's centroid at
// a wall-resolved y+ ~ 1 for Phase 4's flat-plate turbulent boundary layer
// validation; Couette flow's uniform mesh
// (docs/navier-stokes-tracker.md Phase 4) never needed wall-normal
// resolution at all. Reuses build_structured_verification_mesh() (same
// generic topology builder as every other mesh in this file) with a
// displace() that looks up each reference node's precomputed stretched y.
//
// Grows cells geometrically until their cumulative height reaches or
// exceeds H, then rescales every cell's height by the (small, since growth
// stops as soon as H is reached) correction factor needed to land on
// exactly H -- simpler than solving the geometric series for an exact
// growth_ratio/n_y combination, at the cost of the actual first cell height
// coming out slightly different from the literal 'first_cell_height' input.
//
// Input:  n_x               - number of streamwise cells (uniform)
//         L, H              - streamwise/wall-normal domain extents, mesh length units
//         first_cell_height - target height of the wall-adjacent cell, before rescaling
//         growth_ratio      - each successive wall-normal cell's height
//                              divided by the previous one's
// Output: mesh is replaced with the generated nodes/cells/faces/patches;
//         patches are ("bottom"=wall, "left"=inlet, "top"=farfield,
//         "right"=outlet), same ordering/naming as
//         build_structured_verification_mesh()
// Returns: the number of wall-normal cells generated (n_y), so the caller
//          can locate any cell by (i, j) as j*n_x + i
int build_flat_plate_mesh(UnstructuredMesh& mesh, int n_x, double L, double H, double first_cell_height,
                            double growth_ratio) {
    std::vector<double> y_phys = {0.0};
    double h = first_cell_height;
    while (y_phys.back() < H) {
        y_phys.push_back(y_phys.back() + h);
        h *= growth_ratio;
    }
    double rescale = H / y_phys.back();
    for (double& y : y_phys) y *= rescale;

    int n_y = (int)y_phys.size() - 1;
    build_structured_verification_mesh(mesh, n_x, n_y, [&y_phys, n_y, L](double X, double Y) {
        int j = (int)std::lround(Y * n_y);
        return Node{X * L, y_phys[j]};
    });
    return n_y;
}

// Builds an n_x x n_y rectangular domain triangulated by splitting each quad
// into 2 triangles via a checkerboard-alternating diagonal (a "criss-cross"
// pattern) -- deliberately NOT column-aligned, unlike every other mesh in
// this file, to stress-test compute_boundary_layer_profiles()'s
// (WallTraction.h) "most-aligned-neighbor" marching heuristic against a
// genuinely unstructured near-wall triangulation (see
// docs/wall-diagnostics-plan.md's Validation step 4 and that function's own
// disclosed structured-cell-stacking assumption).
//
// Face/adjacency uses the same node-pair edge-matching idea as
// MeshReader.cpp's build_cells_faces_patches() (reimplemented here, not
// shared -- this builds an in-memory verification mesh directly rather than
// parsing a file, the same "own face-building, not the shared parser"
// precedent as build_structured_verification_mesh() elsewhere in this file).
// Unlike this file's structured quad meshes, a triangulated mesh's
// shared-edge topology depends on each quad's checkerboard parity and isn't
// knowable directly from (i, j) alone, which is why this needs the generic
// matcher rather than build_structured_verification_mesh()'s direct
// index-based add_face().
//
// Input:  n_x, n_y - number of quad cells per side before triangulation
//                     (2*n_x*n_y triangles total)
//         L, H      - domain extents, mesh length units
// Output: mesh is replaced with the generated nodes/cells/faces/patches;
//         patches are ("bottom"=wall, "left", "top", "right"), same naming
//         as build_structured_verification_mesh(), assigned directly from
//         node position since this builder knows the domain's exact
//         rectangular extent (no external tagged-edge list, unlike
//         MeshReader's Gmsh-driven version)
void build_criss_cross_triangulated_mesh(UnstructuredMesh& mesh, int n_x, int n_y, double L, double H) {
    auto node_index = [n_x](int i, int j) { return j * (n_x + 1) + i; };

    mesh.patches = {{"bottom"}, {"left"}, {"top"}, {"right"}};
    const int PATCH_BOTTOM = 0, PATCH_LEFT = 1, PATCH_TOP = 2, PATCH_RIGHT = 3;

    mesh.nodes.resize((n_x + 1) * (n_y + 1));
    for (int j = 0; j <= n_y; ++j) {
        for (int i = 0; i <= n_x; ++i) {
            mesh.nodes[node_index(i, j)] = Node{i * L / n_x, j * H / n_y};
        }
    }

    // Split each quad into 2 CCW triangles via a diagonal that alternates
    // with (i + j) parity, so no two quads in a row or column share the same
    // diagonal orientation.
    std::vector<std::vector<int>> cell_node_ids;
    cell_node_ids.reserve((size_t)2 * n_x * n_y);
    for (int j = 0; j < n_y; ++j) {
        for (int i = 0; i < n_x; ++i) {
            int bl = node_index(i, j), br = node_index(i + 1, j);
            int tr = node_index(i + 1, j + 1), tl = node_index(i, j + 1);
            if ((i + j) % 2 == 0) {
                cell_node_ids.push_back({bl, br, tr});
                cell_node_ids.push_back({bl, tr, tl});
            } else {
                cell_node_ids.push_back({bl, br, tl});
                cell_node_ids.push_back({br, tr, tl});
            }
        }
    }

    mesh.cells.resize(cell_node_ids.size());
    std::map<std::pair<int, int>, int> edge_to_face;
    for (size_t ci = 0; ci < cell_node_ids.size(); ++ci) {
        const std::vector<int>& ids = cell_node_ids[ci];
        mesh.cells[ci].node_ids = ids;

        double a2 = 0.0, cx = 0.0, cy = 0.0;
        for (size_t k = 0; k < ids.size(); ++k) {
            const Node& p1 = mesh.nodes[ids[k]];
            const Node& p2 = mesh.nodes[ids[(k + 1) % ids.size()]];
            double cross = p1.x * p2.y - p2.x * p1.y;
            a2 += cross;
            cx += (p1.x + p2.x) * cross;
            cy += (p1.y + p2.y) * cross;
        }
        mesh.cells[ci].volume = std::abs(a2) * 0.5;
        mesh.cells[ci].x_centroid = cx / (3.0 * a2);
        mesh.cells[ci].y_centroid = cy / (3.0 * a2);

        int n = (int)ids.size();
        for (int e = 0; e < n; ++e) {
            int a = ids[e], b = ids[(e + 1) % n];
            auto key = std::make_pair(std::min(a, b), std::max(a, b));
            auto it = edge_to_face.find(key);
            if (it == edge_to_face.end()) {
                Face face;
                face.node1 = a;
                face.node2 = b;
                face.cell_left = (int)ci;
                face.cell_right = -1;
                int face_index = (int)mesh.faces.size();
                mesh.faces.push_back(face);
                mesh.cells[ci].faces.push_back(face_index);
                edge_to_face[key] = face_index;
            } else {
                int face_index = it->second;
                mesh.faces[face_index].cell_right = (int)ci;
                mesh.cells[ci].faces.push_back(face_index);
            }
        }
    }

    const double eps = 1e-9;
    for (Face& face : mesh.faces) {
        if (face.cell_right != -1) continue;
        double y1 = mesh.nodes[face.node1].y, y2 = mesh.nodes[face.node2].y;
        double x1 = mesh.nodes[face.node1].x, x2 = mesh.nodes[face.node2].x;
        if (std::abs(y1) < eps && std::abs(y2) < eps) face.patch_id = PATCH_BOTTOM;
        else if (std::abs(y1 - H) < eps && std::abs(y2 - H) < eps) face.patch_id = PATCH_TOP;
        else if (std::abs(x1) < eps && std::abs(x2) < eps) face.patch_id = PATCH_LEFT;
        else face.patch_id = PATCH_RIGHT;
    }
}

// Runs Phase 4 of the RANS (Spalart-Allmaras) plan's own verification gate,
// the hardest and last-committed phase per this tracker's plan -- originally
// aimed at a flat-plate turbulent boundary layer compared against the
// log-law velocity profile u+ = ln(y+)/kappa + B. Two real attempts at that
// (Re_L = 1e4 and 1e5) both showed nut decaying to a negligible fraction of
// nu with no positive trend -- not a bug, a genuine physical scaling limit:
// SA needs enough local shear relative to wall destruction to sustain
// itself, and reaching that regime for a flat plate needs a Reynolds number
// near real transition values (~5e5-1e6), which is not tractable for this
// project's explicit, uniform-mesh, 2D compressible time-stepping within a
// reasonable run (independently estimated and confirmed with Mathieu; see
// docs/archive/rans-spalart-allmaras-tracker.md Phase 4's notes for the full story).
//
// Given that, the verification target was deliberately relaxed: at the
// (sub-transition) Re_L this project CAN run tractably, SA correctly
// predicts no sustained turbulence, so the mean-flow equations should
// produce an ordinary LAMINAR flat-plate boundary layer -- comparable
// against the classical Pohlhausen quartic approximation,
// u/U = 2*eta - 2*eta^3 + eta^4 (eta = y/delta_99), rather than the
// turbulent log-law originally targeted. delta_99 is extracted from the
// simulation's OWN exit-station profile (interpolated y where u first
// reaches 0.99*U), not from an independent analytic estimate, so this
// isolates the profile's SHAPE (which Blasius/Pohlhausen theory says is
// self-similar regardless of the exact delta) from any small mismatch in
// the boundary-layer growth rate itself. This is still a genuine,
// quantitative comparison against a known reference profile -- not merely
// "did it look boundary-layer-shaped" -- just not the turbulent one Phase 4
// originally set out to reach.
//
// A further finding while choosing H (domain height): delta_99 tracks H
// itself (checked directly -- H = 0.1, 0.2, 0.6 gave delta_99/H ~ 0.92,
// 0.93, 0.73 respectively), not a fixed value independent of H the way a
// genuinely isolated flat-plate boundary layer under an untouched freestream
// should. At this tractable (low) Re_L, molecular viscosity is large enough
// that the viscous-affected region reaches the domain's top Farfield
// boundary well before x = L, so this setup is closer to a developing,
// shallow channel/duct entrance flow than a classical isolated flat plate
// with a clean freestream above it -- the Pohlhausen comparison is still
// meaningful (it approximates a laminar viscous shear layer's self-similar
// SHAPE, not specifically an unbounded external boundary layer), but this
// is a real caveat on what was actually validated, not the textbook flat
// plate the phase's title implies. Making H large enough to avoid this
// entirely would need a proportionally larger Re_L/mesh to keep the
// boundary layer thin relative to it, running back into this phase's
// original compute-cost wall.
//
// Setup: a flat plate of length L, stationary adiabatic no-slip wall along
// the bottom, a Farfield inlet (freestream state, freestream nut = 3*nu) on
// the left, a Farfield "undisturbed freestream" condition (same state) along
// the top (domain height H chosen well above the expected boundary-layer
// thickness at x = L), and a zero-gradient Outflow on the right. SA-noft2
// (no trip term) means the whole plate is treated as turbulent from the
// leading edge -- there is deliberately no laminar-to-turbulent transition
// zone to resolve; nut is free to develop, it just doesn't sustain itself at
// this Re_L, per the finding above.
//
// Input:  none
// Output: prints the exit-station wall-shear/nu_t diagnostics (confirming
//         nu_t stays negligible, i.e. effectively laminar), delta_99, and
//         the Pohlhausen-profile comparison
// Returns: true if at least 3 cells land inside the boundary layer (y <
//          delta_99) and their L2 error against the Pohlhausen profile is
//          below 15%; false otherwise (including if the run diverges)
bool run_verify_flat_plate_boundary_layer() {
    const double L = 1.0, H = 0.2;
    const int n_x = 40;
    const double first_cell_height = 2.949e-3; // sized for y+ ~ 1 at x=L; not load-bearing for the
                                                 // laminar-profile comparison below, but kept from the
                                                 // original turbulent-log-law mesh design since it's
                                                 // already validated and doesn't hurt accuracy either way
    const double growth_ratio = 1.15;

    const double gamma = 1.4, gas_constant = 1.0, prandtl = 0.72, prandtl_t = 0.9, cfl = 0.3;
    const double rho_inf = 1.0, p_inf = 1.0, u_inf = 0.2;
    const double Re_L = 1.0e4;
    const double mu = rho_inf * u_inf * L / Re_L;
    const double nu = mu / rho_inf;
    const double initial_nut = 3.0 * nu;
    const int nsteps = 20000;

    UnstructuredMesh mesh;
    int n_y = build_flat_plate_mesh(mesh, n_x, L, H, first_cell_height, growth_ratio);
    mesh.compute_geometry();

    EulerState freestream = from_primitive(rho_inf, u_inf, 0.0, p_inf, gamma);

    std::vector<RANSBoundaryCondition> bcs(mesh.patches.size());
    bcs[0].ns.type = NSBoundaryType::NoSlipWall; // bottom: the plate
    bcs[1].ns.type = NSBoundaryType::Farfield;   // left: inlet
    bcs[1].ns.farfield_state = freestream;
    bcs[1].farfield_nut = initial_nut;
    bcs[2].ns.type = NSBoundaryType::Farfield;   // top: undisturbed freestream
    bcs[2].ns.farfield_state = freestream;
    bcs[2].farfield_nut = initial_nut;
    bcs[3].ns.type = NSBoundaryType::Outflow;    // right: outlet

    EulerInitialCondition ic;
    ic.mode = EulerICMode::Freestream;
    ic.rho = rho_inf; ic.u = u_inf; ic.v = 0.0; ic.p = p_inf;

    RANSFVMSolver solver(mesh, bcs, gamma, gas_constant, mu, prandtl, prandtl_t, cfl, ic, initial_nut,
                          NumericalFluxScheme::Rusanov, GradientScheme::LeastSquares, 1e-6, 20);

    std::cout << "Flat-plate boundary-layer verification: " << n_x << "x" << n_y << " stretched mesh ("
               << mesh.cells.size() << " cells), Re_L = " << Re_L << ", first cell height = "
               << (mesh.cells[0].y_centroid * 2.0) << "\n";

    int last_step = 0;
    for (; last_step < nsteps; ++last_step) {
        solver.step();
        const EulerResidualNorms& r = solver.residual();
        double nr = solver.nut_residual();
        if (std::isnan(r.rho_u) || std::isinf(r.rho_u) || std::isnan(nr) || std::isinf(nr)) {
            std::cout << "Diverged (NaN/Inf residual) at step " << last_step + 1 << "\n";
            return false;
        }
    }
    std::cout << "Simulation completed across " << last_step << " steps.\n";

    const std::vector<EulerState>& field = solver.field();
    const std::vector<double>& nut_field = solver.nut_field();

    // Wall shear stress at the exit station, from the near-wall cell's own
    // velocity via a simple two-point estimate (u = 0 exactly at the wall,
    // y = 0): du/dy ~ u_cell / y_cell. This is the solver's own computed
    // value, not the a priori Schlichting estimate used to size the mesh.
    int wall_cell = (n_x - 1); // (i=n_x-1, j=0)
    double y_wall_cell = mesh.cells[wall_cell].y_centroid;
    double u_wall_cell = field[wall_cell].rho_u / field[wall_cell].rho;
    double rho_wall_cell = field[wall_cell].rho;
    double nu_t_wall_cell = sa_eddy_viscosity(nut_field[wall_cell], mu / rho_wall_cell);
    double mu_eff_wall = mu + rho_wall_cell * nu_t_wall_cell;
    double tau_wall = mu_eff_wall * u_wall_cell / y_wall_cell;
    double u_tau = std::sqrt(std::fabs(tau_wall) / rho_wall_cell);

    std::cout << "  exit-station wall cell: u = " << u_wall_cell << ", y = " << y_wall_cell
               << ", nu_t = " << nu_t_wall_cell << " (nu = " << nu << ", ratio = " << (nu_t_wall_cell / nu)
               << "), tau_wall = " << tau_wall << ", u_tau = " << u_tau << "\n";

    // delta_99 and the Pohlhausen-profile comparison -- see class comment
    // for why this replaced the originally-targeted turbulent log-law.
    std::vector<double> y_profile(n_y), u_profile(n_y);
    for (int j = 0; j < n_y; ++j) {
        int c = j * n_x + (n_x - 1);
        y_profile[j] = mesh.cells[c].y_centroid;
        u_profile[j] = field[c].rho_u / field[c].rho;
    }

    double delta99 = y_profile.back();
    for (int j = 0; j + 1 < n_y; ++j) {
        if (u_profile[j] < 0.99 * u_inf && u_profile[j + 1] >= 0.99 * u_inf) {
            double t = (0.99 * u_inf - u_profile[j]) / (u_profile[j + 1] - u_profile[j]);
            delta99 = y_profile[j] + t * (y_profile[j + 1] - y_profile[j]);
            break;
        }
    }

    int n_in_bl = 0;
    double sum_sq_error = 0.0, max_error = 0.0;
    for (int j = 0; j < n_y; ++j) {
        double eta = y_profile[j] / delta99;
        if (eta > 1.0) break; // outside the boundary layer -- freestream, nothing to compare
        double u_over_U = u_profile[j] / u_inf;
        double pohlhausen = 2.0 * eta - 2.0 * eta * eta * eta + eta * eta * eta * eta;
        double err = u_over_U - pohlhausen;
        sum_sq_error += err * err;
        max_error = std::max(max_error, std::abs(err));
        ++n_in_bl;
    }

    if (n_in_bl < 3) {
        std::cout << "FAIL: fewer than 3 cells landed inside the boundary layer (y < delta_99); "
                     "cannot meaningfully compare against the Pohlhausen profile.\n";
        return false;
    }

    double l2_error = std::sqrt(sum_sq_error / n_in_bl);
    std::cout << "  delta_99 (interpolated) = " << delta99 << ", " << n_in_bl
               << " cells inside the boundary layer: max |u/U - Pohlhausen| = " << max_error
               << ", L2 = " << l2_error << "\n";

    bool passed = l2_error < 0.15;
    std::cout << (passed ? "PASS" : "FAIL")
               << ": at this sub-transition Re_L, nut must stay negligible and the mean-flow profile "
                  "must match the laminar Pohlhausen approximation (not a turbulent log-law).\n";

    // Phase 2 of docs/wall-diagnostics-plan.md: rerun the exit station
    // through the GENERAL, mesh-topology-walking compute_boundary_layer_profiles()
    // (WallTraction.h) rather than the hand-indexed column walk above, and
    // confirm it reproduces the same delta_99 (regression), then add two
    // checks the hardcoded computation above never had: Blasius laminar skin
    // friction Cf(x) = 0.664/sqrt(Re_x) at the exit station, and Pohlhausen's
    // own closed-form displacement/momentum-thickness ratios
    // (delta*/delta_99 = 3/10, theta/delta_99 = 37/315, from integrating the
    // same quartic profile compared against above) for the delta_99 this
    // run's own profile established.
    std::vector<double> u_field(mesh.cells.size()), v_field(mesh.cells.size());
    std::vector<double> p_field(mesh.cells.size()), rho_field(mesh.cells.size());
    std::vector<double> effective_viscosity(mesh.cells.size());
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        u_field[c] = field[c].rho_u / field[c].rho;
        v_field[c] = field[c].rho_v / field[c].rho;
        p_field[c] = pressure(field[c], gamma);
        rho_field[c] = field[c].rho;
        double nu_t_c = sa_eddy_viscosity(nut_field[c], mu / field[c].rho);
        effective_viscosity[c] = mu + field[c].rho * nu_t_c;
    }

    // Mirrors NavierStokesFVMSolver::build_boundary_fields()'s Dirichlet
    // (NoSlipWall)/zero-order-extrapolation (Farfield/Outflow) convention;
    // RANSFVMSolver's own build_boundary_fields() is private and unreachable
    // from here. NoSlipWall entries are left at 0.0 -- this plate is
    // stationary and adiabatic.
    std::vector<double> boundary_u(mesh.faces.size(), 0.0), boundary_v(mesh.faces.size(), 0.0);
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& f = mesh.faces[i];
        if (f.cell_right != -1) continue;
        if (bcs.at(f.patch_id).ns.type != NSBoundaryType::NoSlipWall) {
            boundary_u[i] = u_field[f.cell_left];
            boundary_v[i] = v_field[f.cell_left];
        }
    }

    int exit_wall_face = -1;
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& f = mesh.faces[i];
        if (f.cell_right == -1 && f.cell_left == wall_cell && f.patch_id == 0) {
            exit_wall_face = (int)i;
            break;
        }
    }

    GradientCalculator grad_calc(mesh, GradientScheme::LeastSquares);
    std::vector<Gradient2> grad_u = grad_calc.compute(mesh, u_field, boundary_u);
    std::vector<Gradient2> grad_v = grad_calc.compute(mesh, v_field, boundary_v);

    std::vector<WallFaceSample> exit_traction = compute_wall_traction(
        mesh, {exit_wall_face}, u_field, v_field, p_field, rho_field, effective_viscosity, boundary_u,
        boundary_v, grad_u, grad_v);
    std::vector<BoundaryLayerProfile> exit_profile = compute_boundary_layer_profiles(
        mesh, {exit_wall_face}, u_field, v_field, boundary_u, boundary_v, u_inf, n_y);

    double delta99_general = exit_profile[0].delta_99;
    double delta99_regression_error = std::abs(delta99_general - delta99) / delta99;

    WallReferenceQuantities blasius_ref;
    blasius_ref.rho_ref = rho_inf;
    blasius_ref.velocity_x_ref = u_inf;
    blasius_ref.velocity_y_ref = 0.0;
    double cf_computed = skin_friction_coefficient(exit_traction[0], blasius_ref);
    double cf_blasius = 0.664 / std::sqrt(Re_L);
    double cf_rel_error = std::abs(std::fabs(cf_computed) - cf_blasius) / cf_blasius;

    double disp_ratio = exit_profile[0].displacement_thickness / delta99_general;
    double mom_ratio = exit_profile[0].momentum_thickness / delta99_general;
    const double disp_ratio_exact = 0.3;
    const double mom_ratio_exact = 37.0 / 315.0;
    double disp_ratio_error = std::abs(disp_ratio - disp_ratio_exact);
    double mom_ratio_error = std::abs(mom_ratio - mom_ratio_exact);

    std::cout << "  General compute_boundary_layer_profiles() at the exit station:\n"
               << "    delta_99 = " << delta99_general << " (hardcoded computation above: " << delta99
               << ", relative diff = " << delta99_regression_error << ")\n"
               << "    Cf = " << cf_computed << " (Blasius Cf = 0.664/sqrt(Re_L) = " << cf_blasius
               << ", relative error = " << cf_rel_error << ")\n"
               << "    displacement_thickness/delta_99 = " << disp_ratio << " (Pohlhausen exact = "
               << disp_ratio_exact << ", |diff| = " << disp_ratio_error << ")\n"
               << "    momentum_thickness/delta_99 = " << mom_ratio << " (Pohlhausen exact = " << mom_ratio_exact
               << ", |diff| = " << mom_ratio_error << ")\n"
               << "    n_cells_marched = " << exit_profile[0].n_cells_marched << " (of " << n_y << " available)\n";

    // Cf's tolerance is deliberately looser than delta_99's regression check
    // or the Couette-flow exact check (--verify-wall-forces): this class
    // comment's own earlier note on H tracking delta_99 already establishes
    // that this setup is closer to a developing, shallow duct-entrance flow
    // than a true isolated Blasius flat plate (the domain top boundary is
    // close enough to interact with the boundary layer) -- 40% captures that
    // known confinement effect without being a meaningless "always passes" bound.
    bool general_bl_passed =
        delta99_regression_error < 0.01 && cf_rel_error < 0.4 && disp_ratio_error < 0.05 && mom_ratio_error < 0.05;
    std::cout << (general_bl_passed ? "PASS" : "FAIL")
               << ": the general compute_boundary_layer_profiles() must reproduce the hardcoded delta_99 above "
                  "and match the Blasius/Pohlhausen closed-form Cf/thickness-ratio values.\n";

    return passed && general_bl_passed;
}

// Runs docs/wall-diagnostics-plan.md's Validation step 4: a stress test of
// compute_boundary_layer_profiles()'s "most-aligned-neighbor" marching
// heuristic against a genuinely unstructured (checkerboard-triangulated, not
// column-aligned) near-wall mesh, to OBSERVE and document its actual
// behavior rather than assume it degrades gracefully -- that function's own
// class comment already discloses this as untested; this is what actually
// tests it.
//
// Setup: a manufactured, exact planar-Couette-like velocity field u = U*y/H,
// v = 0 prescribed directly onto build_criss_cross_triangulated_mesh()'s
// cell centroids -- no solver run, isolating the marching GEOMETRY from any
// solver-convergence question (the same "manufactured field" precedent as
// run_verify_sa_source_terms()'s Test 2). This field's exact delta_99 (where
// u/U first reaches 0.99) is exactly 0.99*H by construction, regardless of
// how the mesh happens to be triangulated -- a march that tracks true
// wall-normal distance must still land on this exact answer.
//
// FINDING (this test fails as written -- see PASS/FAIL message below): the
// failure mode is NOT lateral wandering (the "most-aligned-neighbor" choice
// itself consistently favors the cell stacked above, even across the
// diagonal cuts) -- it's a systematic ~20-25% delta_99 error, essentially
// identical across every wall face despite their differing local diagonal
// orientation, which points at cumulative distance accounting rather than
// cell selection. compute_boundary_layer_profiles() accumulates
// face_normal_distance() at each step (per docs/wall-diagnostics-plan.md's
// own "Boundary-layer thickness" methodology), which projects onto THAT
// FACE's own normal -- correct when every crossed face's normal already
// equals the march's fixed global direction (true for a structured quad
// stack), but not when a step crosses a diagonal face whose normal is
// tilted away from vertical, where it silently sums a distance measured
// along the wrong axis instead of the true wall-normal one. This is exactly
// the "structured cell stacking" assumption docs/wall-diagnostics-plan.md
// already discloses -- now demonstrated with a concrete mechanism and
// number, not left hypothetical.
//
// Input:  none
// Output: prints the max/mean relative delta_99 error over every wall face,
//         and how many faces' marches hit max_cells_per_march without ever
//         crossing 0.99*U
// Returns: true if the max relative delta_99 error stays under 10%; false
//          otherwise -- either result is a legitimate, evidence-based
//          finding for docs/wall-diagnostics-plan.md's Phase 3 to act on,
//          not a foregone conclusion
bool run_verify_bl_marching_unstructured() {
    const int n_x = 12, n_y = 10;
    const double L = 1.0, H = 0.3;
    const double U = 1.0;

    UnstructuredMesh mesh;
    build_criss_cross_triangulated_mesh(mesh, n_x, n_y, L, H);
    mesh.compute_geometry();

    std::vector<double> u_field(mesh.cells.size()), v_field(mesh.cells.size(), 0.0);
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        u_field[c] = U * mesh.cells[c].y_centroid / H;
    }
    std::vector<double> boundary_u(mesh.faces.size(), 0.0), boundary_v(mesh.faces.size(), 0.0);

    std::vector<int> wall_faces;
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        if (mesh.faces[i].patch_id == 0) wall_faces.push_back((int)i); // "bottom" patch
    }

    const int max_cells_per_march = 3 * n_y; // generous enough that every march below genuinely crosses 0.99*U
    std::vector<BoundaryLayerProfile> profiles = compute_boundary_layer_profiles(
        mesh, wall_faces, u_field, v_field, boundary_u, boundary_v, U, max_cells_per_march);

    double delta99_exact = 0.99 * H;
    double max_rel_error = 0.0, sum_rel_error = 0.0;
    int n_hit_cap = 0;
    for (const auto& p : profiles) {
        double rel_error = std::abs(p.delta_99 - delta99_exact) / delta99_exact;
        max_rel_error = std::max(max_rel_error, rel_error);
        sum_rel_error += rel_error;
        if (p.n_cells_marched >= max_cells_per_march) ++n_hit_cap;
    }
    double mean_rel_error = sum_rel_error / profiles.size();

    std::cout << "Boundary-layer marching unstructured-mesh stress test: " << n_x << "x" << n_y
               << " criss-cross-triangulated mesh (" << mesh.cells.size() << " triangles), " << wall_faces.size()
               << " wall faces, manufactured field u = U*y/H (exact delta_99 = " << delta99_exact << ")\n"
               << "  max relative delta_99 error  = " << max_rel_error << "\n"
               << "  mean relative delta_99 error = " << mean_rel_error << "\n"
               << "  marches that hit max_cells_per_march without crossing 0.99*U: " << n_hit_cap << " / "
               << profiles.size() << "\n";

    bool passed = max_rel_error < 0.10;
    std::cout << (passed ? "PASS" : "FAIL")
               << ": on this checkerboard-triangulated near-wall mesh, "
               << (passed ? "the marching heuristic stayed accurate enough to reproduce the exact delta_99 within "
                            "10%."
                          : "delta_99 misses the exact value by more than 10%, from a systematic "
                            "distance-accounting bias when the march crosses a diagonal face whose normal isn't "
                            "aligned with the fixed march direction (see class comment) -- the point-location "
                            "alternative docs/wall-diagnostics-plan.md describes is now a concrete follow-up, "
                            "not a hypothetical one.")
               << "\n";
    return passed;
}

// Runs the same stress test as run_verify_bl_marching_unstructured() above
// (identical checkerboard-triangulated mesh and manufactured
// planar-Couette-like field, whose exact delta_99 = 0.99*H is known
// regardless of triangulation), but through
// compute_boundary_layer_profiles_point_location() (WallTraction.h) instead
// of the cell-marching version -- confirming the point-location alternative
// doesn't share cell-marching's demonstrated failure mode on this exact mesh
// (see docs/wall-diagnostics-plan.md's "Phase 2 stress-test finding" and
// Phase 3), since it never makes a neighbor-alignment choice at all.
//
// Input:  none
// Output: prints the max/mean relative delta_99 error over every wall face,
//         and how many faces' samples ran off the mesh domain or hit
//         n_samples without ever crossing 0.99*U
// Returns: true if the max relative delta_99 error stays under 1% -- much
//          tighter than cell-marching's own 10% bar here, since
//          point-location has no topology-dependent approximation at all on
//          this mesh, only sampling-resolution error from a finite n_samples
bool run_verify_bl_point_location() {
    const int n_x = 12, n_y = 10;
    const double L = 1.0, H = 0.3;
    const double U = 1.0;

    UnstructuredMesh mesh;
    build_criss_cross_triangulated_mesh(mesh, n_x, n_y, L, H);
    mesh.compute_geometry();

    std::vector<double> u_field(mesh.cells.size()), v_field(mesh.cells.size(), 0.0);
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        u_field[c] = U * mesh.cells[c].y_centroid / H;
    }
    std::vector<double> boundary_u(mesh.faces.size(), 0.0), boundary_v(mesh.faces.size(), 0.0);

    std::vector<int> wall_faces;
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        if (mesh.faces[i].patch_id == 0) wall_faces.push_back((int)i); // "bottom" patch
    }

    const double max_distance = 1.05 * H; // comfortably past the exact delta_99 = 0.99*H
    const int n_samples = 400;
    std::vector<BoundaryLayerProfile> profiles = compute_boundary_layer_profiles_point_location(
        mesh, wall_faces, u_field, v_field, boundary_u, boundary_v, U, max_distance, n_samples);

    double delta99_exact = 0.99 * H;
    double max_rel_error = 0.0, sum_rel_error = 0.0;
    int n_hit_cap = 0;
    for (const auto& p : profiles) {
        double rel_error = std::abs(p.delta_99 - delta99_exact) / delta99_exact;
        max_rel_error = std::max(max_rel_error, rel_error);
        sum_rel_error += rel_error;
        if (p.n_cells_marched >= n_samples) ++n_hit_cap;
    }
    double mean_rel_error = sum_rel_error / profiles.size();

    std::cout << "Boundary-layer point-location unstructured-mesh stress test: " << n_x << "x" << n_y
               << " criss-cross-triangulated mesh (" << mesh.cells.size() << " triangles), " << wall_faces.size()
               << " wall faces, manufactured field u = U*y/H (exact delta_99 = " << delta99_exact << ")\n"
               << "  max relative delta_99 error  = " << max_rel_error << "\n"
               << "  mean relative delta_99 error = " << mean_rel_error << "\n"
               << "  samples that ran off the mesh or hit n_samples without crossing 0.99*U: " << n_hit_cap << " / "
               << profiles.size() << "\n";

    bool passed = max_rel_error < 0.01;
    std::cout << (passed ? "PASS" : "FAIL")
               << ": point-location must reproduce the exact delta_99 on this checkerboard-triangulated mesh to "
                  "within 1%, since it has no neighbor-topology dependence at all (only finite-sample-resolution "
                  "error).\n";
    return passed;
}

// Runs Phases 0-1 of the gradient-reconstruction verification, on a known
// linear scalar field phi(x,y) = a + b*x + c*y over a deliberately
// non-orthogonal ("wavy") mesh:
//   Phase 0 - checks that each GradientScheme reconstructs the exact
//     analytic cell gradient (b, c). Least-Squares must recover it to near
//     machine precision on any mesh (by construction of its normal
//     equations); Green-Gauss is only exact when face midpoints lie on the
//     line joining their two cell centroids, which is false on this mesh,
//     so its reported error is expected and informative, not a failure.
//   Phase 1 - using the (near-exact) Least-Squares cell gradients, checks
//     that face_gradient()'s corrected normal derivative also reproduces
//     the exact analytic value at every face, and reports what today's
//     UnstructuredFVMSolver-style naive two-point difference (no
//     non-orthogonality correction) gives on the same faces for comparison.
//
// Input:  none
// Output: prints a per-scheme/per-check max-error report to stdout
// Returns: true if both Least-Squares's cell-gradient error and the
//          corrected face gradient's error are below a tight tolerance
//          (1e-9); false otherwise, which would indicate an actual bug (not
//          a mesh-skewness artifact -- Green-Gauss's and the naive
//          two-point form's errors are expected and don't gate this)
bool run_verify_gradient() {
    const int n = 8;
    const double wobble = 0.15;
    const double a = 1.7, b = 2.3, c = -1.1;

    UnstructuredMesh mesh;
    build_skewed_verification_mesh(mesh, n, wobble);
    mesh.compute_geometry();

    auto analytic = [&](double x, double y) { return a + b * x + c * y; };

    std::vector<double> phi(mesh.cells.size());
    for (size_t i = 0; i < mesh.cells.size(); ++i) {
        phi[i] = analytic(mesh.cells[i].x_centroid, mesh.cells[i].y_centroid);
    }

    std::vector<double> boundary_phi(mesh.faces.size(), 0.0);
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        if (mesh.faces[i].cell_right == -1) {
            boundary_phi[i] = analytic(mesh.faces[i].x_mid, mesh.faces[i].y_mid);
        }
    }

    std::cout << "Gradient verification: " << n << "x" << n << " skewed mesh ("
               << mesh.cells.size() << " cells), phi = " << a << " + " << b << "*x + " << c << "*y\n"
               << "Exact gradient: (" << b << ", " << c << ")\n";

    double lsq_max_error = 0.0;
    std::vector<Gradient2> lsq_grad; // kept for the face-gradient check below
    for (GradientScheme scheme : {GradientScheme::GreenGauss, GradientScheme::LeastSquares}) {
        GradientCalculator calc(mesh, scheme);
        std::vector<Gradient2> grad = calc.compute(mesh, phi, boundary_phi);

        double max_error = 0.0;
        for (const auto& g : grad) {
            max_error = std::max(max_error, std::max(std::abs(g.dphidx - b), std::abs(g.dphidy - c)));
        }

        const char* name = (scheme == GradientScheme::GreenGauss) ? "Green-Gauss" : "Least-Squares";
        std::cout << "  " << name << ": max gradient error = " << max_error << "\n";

        if (scheme == GradientScheme::LeastSquares) {
            lsq_max_error = max_error;
            lsq_grad = std::move(grad);
        }
    }
    bool cell_gradient_passed = lsq_max_error < 1e-9;
    std::cout << (cell_gradient_passed ? "PASS" : "FAIL")
               << ": Least-Squares must reproduce a linear field's gradient near machine precision.\n";

    // Phase 1: corrected face-normal derivative, using the (near-exact)
    // Least-Squares cell gradients as face_gradient()'s input. Compared
    // against both the exact analytic normal derivative and today's
    // UnstructuredFVMSolver-style naive two-point difference (no
    // non-orthogonality correction at all), to show the naive form's real
    // error on this non-orthogonal mesh, not just assert one exists.
    double corrected_max_error = 0.0, naive_max_error = 0.0;
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        const Face& face = mesh.faces[i];
        int cl = face.cell_left, cr = face.cell_right;
        bool is_boundary = (cr == -1);

        double phi_R = is_boundary ? boundary_phi[i] : phi[cr];
        Gradient2 grad_R = is_boundary ? Gradient2{} : lsq_grad[cr];

        FaceGradient fg = face_gradient(mesh, face, phi[cl], phi_R, lsq_grad[cl], grad_R, boundary_phi[i]);

        double x_R = is_boundary ? face.x_mid : mesh.cells[cr].x_centroid;
        double y_R = is_boundary ? face.y_mid : mesh.cells[cr].y_centroid;
        double dist = std::hypot(x_R - mesh.cells[cl].x_centroid, y_R - mesh.cells[cl].y_centroid);
        double naive_dphidn = (phi_R - phi[cl]) / dist;

        double exact_dphidn = b * face.nx + c * face.ny;
        corrected_max_error = std::max(corrected_max_error, std::abs(fg.dphidn - exact_dphidn));
        naive_max_error = std::max(naive_max_error, std::abs(naive_dphidn - exact_dphidn));
    }

    std::cout << "  Corrected face gradient: max normal-derivative error = " << corrected_max_error << "\n"
               << "  Naive two-point (no correction): max normal-derivative error = " << naive_max_error << "\n";
    bool face_gradient_passed = corrected_max_error < 1e-9;
    std::cout << (face_gradient_passed ? "PASS" : "FAIL")
               << ": corrected face gradient must reproduce a linear field's normal derivative near machine precision.\n";

    return cell_gradient_passed && face_gradient_passed;
}

// Runs Phase 2's advection-diffusion verification: steady 1D advection with
// diffusion of a passive scalar between two Dirichlet boundaries, on the
// same wavy (non-orthogonal) mesh used by run_verify_gradient(), with the
// mean flow direction NOT aligned with the mesh's internal faces (the
// interior is warped even though the outer boundary is a straight
// rectangle -- see build_skewed_verification_mesh()), so the corrected
// diffusive face gradient from Phase 1 is genuinely exercised rather than
// happening to align with an orthogonal mesh.
//
// Setup: u_adv = (U, 0), alpha constant; Dirichlet phi=0 on "left" (x=0),
// Dirichlet phi=1 on "right" (x=1); zero-flux Neumann on "bottom"/"top" (the
// exact solution is independent of y, and those boundary faces are exactly
// horizontal since the perturbation vanishes on the outer boundary, so
// u_adv has exactly zero normal component there -- no inflow/outflow to get
// wrong). The known steady-state solution of u*dphi/dx = alpha*d2phi/dx2
// with phi(0)=0, phi(1)=1 is the classic 1D advection-diffusion profile
// phi(x) = (exp(Pe*x) - 1) / (exp(Pe) - 1), Pe = U/alpha.
//
// Input:  none
// Output: prints the steady-state max/L2 error against the exact profile,
//         plus the final residual norm (as a steady-state sanity check)
// Returns: true if the L2 relative error is below a tolerance appropriate
//          for a first-order-upwind scheme at this resolution (10%); false
//          otherwise
bool run_verify_advection_diffusion() {
    const int n = 16;
    const double wobble = 0.15;
    const double U = 1.0, alpha = 0.1;
    const double Pe = U / alpha;
    const double dt = 0.0005;
    const int nsteps = 40000;

    UnstructuredMesh mesh;
    build_skewed_verification_mesh(mesh, n, wobble);
    mesh.compute_geometry();

    mesh.patches[1].type = BoundaryType::Dirichlet; mesh.patches[1].value = 0.0; // left
    mesh.patches[3].type = BoundaryType::Dirichlet; mesh.patches[3].value = 1.0; // right
    mesh.patches[0].type = BoundaryType::Neumann; mesh.patches[0].value = 0.0;   // bottom
    mesh.patches[2].type = BoundaryType::Neumann; mesh.patches[2].value = 0.0;   // top

    AdvectionDiffusionFVMSolver solver(mesh, alpha, U, 0.0, dt, GradientScheme::LeastSquares, 0.0, 0.0);
    solver.run(nsteps);

    auto exact = [&](double x) { return (std::exp(Pe * x) - 1.0) / (std::exp(Pe) - 1.0); };

    double max_error = 0.0, sum_sq_error = 0.0;
    const std::vector<double>& phi = solver.field();
    for (size_t i = 0; i < mesh.cells.size(); ++i) {
        double err = phi[i] - exact(mesh.cells[i].x_centroid);
        max_error = std::max(max_error, std::abs(err));
        sum_sq_error += err * err;
    }
    double l2_relative_error = std::sqrt(sum_sq_error / mesh.cells.size());

    std::cout << "Advection-diffusion verification: " << n << "x" << n << " wavy mesh (" << mesh.cells.size()
               << " cells), Pe = " << Pe << ", after " << nsteps << " steps (residual = "
               << solver.residual_norm() << ")\n"
               << "  max error vs exact profile = " << max_error << "\n"
               << "  L2 error vs exact profile  = " << l2_relative_error << "\n";

    bool passed = l2_relative_error < 0.1;
    std::cout << (passed ? "PASS" : "FAIL")
               << ": L2 error against the exact 1D advection-diffusion profile must stay under 10%.\n";
    return passed;
}

// Runs Phase 3's own verification gate for NavierStokesFVMSolver, ahead of
// (and distinct from) Phase 4's Couette flow validation: a uniform
// freestream state, prescribed identically as the Farfield BC on every side
// of the same wavy (non-orthogonal) mesh used by Phases 0-2, must stay
// EXACTLY at that uniform state for any number of steps and any mu > 0.
//
// Why this is the right thing to check now: a spatially uniform state has
// zero gradient everywhere, so both the inviscid flux (already established
// by EulerFVMSolver) and -- newly, for this phase -- the viscous stress
// tensor and heat flux must come out to exactly zero at every face,
// regardless of the mesh's non-orthogonality. This is the direct NS
// analogue of Phase 0/1's "reconstruct a known field exactly" tests, now
// applied to the new viscous flux assembly (corrected_face_gradient_vector()
// in NavierStokesFVMSolver.cpp) rather than to a standalone scalar gradient.
// It does NOT test whether the viscous terms produce physically correct
// shear/heat transport for a NON-uniform flow -- that's what Phase 4's
// Couette test is for.
//
// Input:  none
// Output: prints the max deviation from the initial uniform state after
//         running with mu > 0
// Returns: true if that deviation stays below a tight tolerance (1e-9);
//          false otherwise, which would indicate the viscous flux assembly
//          spuriously injects flux where there is no gradient at all
bool run_verify_navier_stokes_uniform() {
    const int n = 8;
    const double wobble = 0.15;
    const double gamma = 1.4, gas_constant = 1.0, mu = 0.5, prandtl = 0.72, cfl = 0.3;
    const double rho0 = 1.0, u0 = 1.0, v0 = 0.5, p0 = 1.0;
    const int nsteps = 20;

    UnstructuredMesh mesh;
    build_skewed_verification_mesh(mesh, n, wobble);
    mesh.compute_geometry();

    EulerState uniform_state = from_primitive(rho0, u0, v0, p0, gamma);
    std::vector<NSBoundaryCondition> bcs(mesh.patches.size());
    for (auto& bc : bcs) {
        bc.type = NSBoundaryType::Farfield;
        bc.farfield_state = uniform_state;
    }

    EulerInitialCondition ic;
    ic.mode = EulerICMode::Freestream;
    ic.rho = rho0; ic.u = u0; ic.v = v0; ic.p = p0;

    NavierStokesFVMSolver solver(mesh, bcs, gamma, gas_constant, mu, prandtl, cfl, ic, NumericalFluxScheme::Rusanov,
                                  GradientScheme::LeastSquares, 1e-6, 20);
    solver.run(nsteps);

    double max_error = 0.0;
    for (const EulerState& s : solver.field()) {
        max_error = std::max(max_error, std::abs(s.rho - uniform_state.rho));
        max_error = std::max(max_error, std::abs(s.rho_u - uniform_state.rho_u));
        max_error = std::max(max_error, std::abs(s.rho_v - uniform_state.rho_v));
        max_error = std::max(max_error, std::abs(s.E - uniform_state.E));
    }

    std::cout << "Navier-Stokes uniform-flow verification: " << n << "x" << n << " wavy mesh (" << mesh.cells.size()
               << " cells), mu = " << mu << ", after " << nsteps << " steps\n"
               << "  max deviation from initial uniform state = " << max_error << "\n";

    bool passed = max_error < 1e-9;
    std::cout << (passed ? "PASS" : "FAIL")
               << ": a uniform state (zero gradient everywhere) must stay exactly uniform.\n";
    return passed;
}

// Runs Phase 4's Couette flow validation: the real viscous-physics test that
// Phase 3's uniform-flow check deliberately deferred. Classic planar Couette
// flow -- a stationary bottom wall, a top wall moving at U, no imposed
// pressure gradient -- has the exact steady analytic solution u(y) = U*y/H,
// v = 0, independent of x. "left"/"right" are set to zero-gradient outflow
// (appropriate since the exact solution has no x-dependence at all -- there
// is nothing for those boundaries to impose), so this directly exercises
// whether the corrected viscous shear stress (Phase 1's face gradient,
// generalized to a full tensor in Phase 3) reproduces real, physically
// meaningful shear transport across a genuinely non-orthogonal interior
// mesh, not just "vanishes when there's nothing to diffuse" (Phase 3's gate).
//
// Uses build_sheared_verification_mesh(), NOT the "wavy" mesh from Phases
// 0-2 -- an earlier attempt with the wavy mesh gave a ~50% L2 error that
// didn't shrink with 5x more time steps, i.e. a genuine (wrong) steady
// state, not merely a slow one. Diagnosis: the wavy mesh's own distortion
// varies with x (peaks at x=0.5), so the DISCRETIZATION isn't
// x-translation-invariant even though the continuous Couette solution is,
// and that mismatch injects a spurious x-varying truncation-error "forcing"
// -- visible as a parasitic Poiseuille-shaped (parabolic, zero at both
// walls) deviation from the exact linear profile. A pure shear doesn't have
// that problem (every cell at a given row is an identical parallelogram
// just translated in x) while still being genuinely non-orthogonal in
// Jasak's sense (see build_sheared_verification_mesh()'s comment).
//
// Low Mach number (U/c ~ 0.08 here) keeps this compressible solver's result
// close to the incompressible analytic solution; both walls are adiabatic
// (isothermal wall choice doesn't matter for the velocity profile at all --
// this solver's momentum equation doesn't depend on temperature since mu is
// a fixed constant, not temperature-dependent -- so there is deliberately no
// viscous-heating/temperature check here, matching the originally agreed scope).
//
// Input:  none
// Output: prints the steady-state max/L2 error against the exact profile
//         (both absolute and relative to U), plus max |v| as a secondary
//         sanity check (must stay ~0) and the final rho_u residual
// Returns: true if the L2 relative error is below 5%; false otherwise
bool run_verify_couette() {
    const int n_x = 1, n_y = 16;
    const double shear = 0.3;
    const double gamma = 1.4, gas_constant = 1.0, mu = 0.02, prandtl = 0.72, cfl = 0.3;
    const double rho0 = 1.0, p0 = 1.0, U = 0.1;
    const double H = 1.0;
    const int nsteps = 20000;

    UnstructuredMesh mesh;
    build_sheared_verification_mesh(mesh, n_x, n_y, shear);
    mesh.compute_geometry();

    std::vector<NSBoundaryCondition> bcs(mesh.patches.size());
    bcs[0].type = NSBoundaryType::NoSlipWall;             // bottom: stationary
    bcs[1].type = NSBoundaryType::Outflow;                // left
    bcs[2].type = NSBoundaryType::NoSlipWall;             // top: moving
    bcs[2].wall_u = U;
    bcs[3].type = NSBoundaryType::Outflow;                // right

    EulerInitialCondition ic;
    ic.mode = EulerICMode::Freestream;
    ic.rho = rho0; ic.u = 0.0; ic.v = 0.0; ic.p = p0;

    NavierStokesFVMSolver solver(mesh, bcs, gamma, gas_constant, mu, prandtl, cfl, ic, NumericalFluxScheme::Rusanov,
                                  GradientScheme::LeastSquares, 1e-6, 20);
    int last_step = 0;
    for (; last_step < nsteps; ++last_step) {
        solver.step();
        const EulerResidualNorms& r = solver.residual();
        if (std::isnan(r.rho_u) || std::isinf(r.rho_u)) {
            std::cout << "Diverged (NaN/Inf residual) at step " << last_step + 1 << "\n";
            return false;
        }
    }
    std::cout << "Simulation completed across " << last_step << " steps.\n";

    const std::vector<EulerState>& field = solver.field();
    double max_error = 0.0, sum_sq_error = 0.0, max_v = 0.0;
    for (size_t i = 0; i < mesh.cells.size(); ++i) {
        double u = field[i].rho_u / field[i].rho;
        double v = field[i].rho_v / field[i].rho;
        double exact_u = U * mesh.cells[i].y_centroid / H;
        double err = u - exact_u;
        max_error = std::max(max_error, std::abs(err));
        sum_sq_error += err * err;
        max_v = std::max(max_v, std::abs(v));
    }
    double l2_error = std::sqrt(sum_sq_error / mesh.cells.size());

    const EulerResidualNorms& r = solver.residual();
    std::cout << "Couette flow verification: " << n_x << "x" << n_y << " sheared mesh (" << mesh.cells.size()
               << " cells), U = " << U << ", mu = " << mu << ", after " << nsteps
               << " steps (rho_u residual = " << r.rho_u << ")\n"
               << "  max |u - exact| = " << max_error << " (relative to U: " << max_error / U << ")\n"
               << "  L2  |u - exact| = " << l2_error << " (relative to U: " << l2_error / U << ")\n"
               << "  max |v| (should be ~0) = " << max_v << "\n";

    bool passed = (l2_error / U) < 0.05;
    std::cout << (passed ? "PASS" : "FAIL")
               << ": velocity profile must match the analytic linear Couette profile u(y) = U*y/H within 5% (L2).\n";
    return passed;
}

// Runs Phase 1's own verification gate for WallTraction.h's force/Cf/Cp/y+/Cm
// plumbing (see docs/wall-diagnostics-plan.md), reusing run_verify_couette()'s
// exact 1x16 sheared-mesh planar Couette setup for its closed-form checks:
//   - tau_wall = mu*U/H at both walls (exact, no imposed pressure gradient)
//   - p uniform everywhere (Couette flow has uniform pressure by
//     construction -- checked directly as a raw pressure deviation from p0,
//     not through Cp: this setup's small U (needed to keep the compressible
//     solver close to the incompressible analytic solution -- see
//     run_verify_couette()'s own comment) makes the dynamic pressure tiny,
//     so Cp's normalization amplifies the SAME small residual
//     ~1-2% pressure deviation that run_verify_couette() already tolerates
//     via its own 5% velocity-profile bound (a documented, pre-existing
//     secular-drift artifact of this exact setup -- adiabatic walls plus
//     Outflow's imperfect periodicity, see NavierStokesFVMSolver.h's own
//     Outflow caveat -- not a WallTraction bug) into a Cp on the order of a
//     few, not near 0; Cp is still printed for visibility but is not what
//     this gate checks
//   - y+ computable by hand: tau_wall, rho, mu, and the first-cell
//     wall-normal distance (H/(2*n_y), exact on this mesh -- same reasoning
//     as run_verify_wall_distance()'s check) are all known
//   - Cm about the domain's vertical center (0, H/2): both walls carry equal
//     shear (tau_wall comes out numerically identical at both walls, since
//     the outward normal AND tangent direction both flip sign between them)
//     and equal-and-opposite uniform-pressure forces, offset by +-H/2 -- a
//     real, sign-sensitive nonzero analytic moment (assuming p == p0
//     exactly) that a sign error would still fail, not a degenerate zero one
//     would still pass; checked with the same relative tolerance as the raw
//     pressure check above, since it inherits that same small drift
//
// Input:  none
// Output: prints tau_wall/p/Cp/y+ at both walls and Cm's domain total
//         against their analytic values
// Returns: true if every quantity matches its analytic value within
//          tolerance; false otherwise
bool run_verify_wall_forces() {
    const int n_x = 1, n_y = 16;
    const double shear = 0.3;
    const double gamma = 1.4, gas_constant = 1.0, mu = 0.02, prandtl = 0.72, cfl = 0.3;
    const double rho0 = 1.0, p0 = 1.0, U = 0.1;
    const double H = 1.0;
    const int nsteps = 20000;

    UnstructuredMesh mesh;
    build_sheared_verification_mesh(mesh, n_x, n_y, shear);
    mesh.compute_geometry();

    std::vector<NSBoundaryCondition> bcs(mesh.patches.size());
    bcs[0].type = NSBoundaryType::NoSlipWall;             // bottom: stationary
    bcs[1].type = NSBoundaryType::Outflow;                // left
    bcs[2].type = NSBoundaryType::NoSlipWall;             // top: moving
    bcs[2].wall_u = U;
    bcs[3].type = NSBoundaryType::Outflow;                // right

    EulerInitialCondition ic;
    ic.mode = EulerICMode::Freestream;
    ic.rho = rho0; ic.u = 0.0; ic.v = 0.0; ic.p = p0;

    NavierStokesFVMSolver solver(mesh, bcs, gamma, gas_constant, mu, prandtl, cfl, ic, NumericalFluxScheme::Rusanov,
                                  GradientScheme::LeastSquares, 1e-6, 20);
    for (int step = 0; step < nsteps; ++step) {
        solver.step();
    }

    std::vector<int> wall_faces = collect_wall_faces(mesh, bcs);

    WallReferenceQuantities ref;
    ref.rho_ref = rho0;
    ref.velocity_x_ref = U;
    ref.velocity_y_ref = 0.0;
    ref.p_ref = p0;
    ref.length_ref = H;
    ref.moment_reference_x = 0.0;
    ref.moment_reference_y = H / 2.0; // domain's vertical center

    std::vector<WallFaceSample> samples = solver.compute_wall_traction_samples(wall_faces);

    double tau_exact = mu * U / H;
    double y_wall_normal_exact = H / (2.0 * n_y);
    double u_tau_exact = std::sqrt(tau_exact / rho0);
    double yplus_exact = rho0 * u_tau_exact * y_wall_normal_exact / mu;

    double max_tau_rel_error = 0.0, max_p_rel_error = 0.0, max_cp = 0.0, max_yplus_rel_error = 0.0;
    for (const WallFaceSample& s : samples) {
        double cp = pressure_coefficient(s, ref);
        double yplus = wall_y_plus(s);
        max_tau_rel_error = std::max(max_tau_rel_error, std::abs(std::fabs(s.tau_wall) - tau_exact) / tau_exact);
        max_p_rel_error = std::max(max_p_rel_error, std::abs(s.p - p0) / p0);
        max_cp = std::max(max_cp, std::abs(cp));
        max_yplus_rel_error = std::max(max_yplus_rel_error, std::abs(yplus - yplus_exact) / yplus_exact);
    }

    std::vector<WallForceReport> reports = compute_wall_forces(mesh, samples, ref);
    const WallForceReport* total = nullptr;
    for (const auto& r : reports) {
        if (r.patch_id == -1) total = &r;
    }

    // Exact moment: both walls' friction force is (area * tau_exact,
    // signed per the traction convention above) plus each wall's own
    // uniform-pressure force, taken about (0, H/2). Both walls span x in
    // [0,1] displaced by shear*Y, so their midpoints/normals/tangents can be
    // worked out directly from build_sheared_verification_mesh()'s mapping
    // x = X + shear*Y, y = Y.
    double bottom_mid_x = 0.5, bottom_mid_y = 0.0;
    double top_mid_x = 0.5 + shear, top_mid_y = H;
    double signed_tau = -tau_exact; // sign convention derived from this traction/tangent setup (see class comment)
    auto moment_of = [&](double mx, double my, double nx, double ny, double tx, double ty) {
        double rx = mx - ref.moment_reference_x, ry = my - ref.moment_reference_y;
        double fx = signed_tau * tx + (-p0) * nx, fy = signed_tau * ty + (-p0) * ny;
        return rx * fy - ry * fx;
    };
    double moment_exact = moment_of(bottom_mid_x, bottom_mid_y, 0.0, -1.0, 1.0, 0.0) +
                           moment_of(top_mid_x, top_mid_y, 0.0, 1.0, -1.0, 0.0);
    double moment_rel_error = total ? std::abs(total->moment - moment_exact) / std::abs(moment_exact)
                                     : std::numeric_limits<double>::infinity();

    std::cout << "Wall-forces verification: " << n_x << "x" << n_y << " sheared mesh (" << mesh.cells.size()
               << " cells), U = " << U << ", mu = " << mu << ", after " << nsteps << " steps\n"
               << "  max |tau_wall - exact| / exact = " << max_tau_rel_error << " (exact = " << tau_exact << ")\n"
               << "  max |p - p0| / p0 = " << max_p_rel_error << " (max |Cp| = " << max_cp << ", printed for "
                  "visibility only -- see class comment on why it isn't near 0 for this reference velocity)\n"
               << "  max |y+ - exact| / exact = " << max_yplus_rel_error << " (exact = " << yplus_exact << ")\n"
               << "  |moment - exact| / |exact| = " << moment_rel_error << " (computed = "
               << (total ? total->moment : 0.0) << ", exact = " << moment_exact << ")\n";

    bool passed = max_tau_rel_error < 1e-6 && max_p_rel_error < 0.05 && max_yplus_rel_error < 1e-6 &&
                  moment_rel_error < 0.05;
    std::cout << (passed ? "PASS" : "FAIL")
               << ": tau_wall/p/y+/moment must match their closed-form Couette-flow values (5% relative "
                  "tolerance on p/moment, matching run_verify_couette()'s own bound; near machine precision "
                  "on tau_wall/y+, which don't depend on the pressure field at all).\n";
    return passed;
}

// Regression check for the compute_dt() viscous-length-scale fix (see
// docs/ns-cfl-margin-and-farfield-bc-findings.md and
// docs/navier-stokes-tracker.md Phase 6): a laminar NavierStokesFVMSolver run
// on a wall-normal-stretched, boundary-layer-style mesh -- reusing
// build_flat_plate_mesh() (RANS Phase 4's already-validated anisotropic-mesh
// generator) rather than this file's other meshes, none of which have a
// genuinely thin first cell -- at a cfl/first_cell_height/mu combination in
// the same regime as the real airfoil case that diverged before the fix
// (mu = 0.02, cfl = 0.5). This isn't a quantitative accuracy check like
// Couette flow; it only checks that a run this anisotropic doesn't diverge.
//
// Input:  none
// Output: prints the mesh/case parameters and whether the run diverged
// Returns: true if the run completes all steps without a NaN/Inf residual;
//          false (with the step at which it diverged) otherwise
bool run_verify_ns_stretched_cfl() {
    const double L = 1.0, H = 0.05;
    const int n_x = 10;
    const double first_cell_height = 1.0e-4;
    const double growth_ratio = 1.05;

    const double gamma = 1.4, gas_constant = 1.0, prandtl = 0.72, cfl = 2.0;
    const double rho_inf = 1.0, p_inf = 1.0, u_inf = 0.5, mu = 0.02;
    const int nsteps = 2000;

    UnstructuredMesh mesh;
    int n_y = build_flat_plate_mesh(mesh, n_x, L, H, first_cell_height, growth_ratio);
    mesh.compute_geometry();

    EulerState freestream = from_primitive(rho_inf, u_inf, 0.0, p_inf, gamma);

    std::vector<NSBoundaryCondition> bcs(mesh.patches.size());
    bcs[0].type = NSBoundaryType::NoSlipWall; // bottom: the plate
    bcs[1].type = NSBoundaryType::Farfield;   // left: inlet
    bcs[1].farfield_state = freestream;
    bcs[2].type = NSBoundaryType::Farfield;   // top: undisturbed freestream
    bcs[2].farfield_state = freestream;
    bcs[3].type = NSBoundaryType::Outflow;    // right: outlet

    EulerInitialCondition ic;
    ic.mode = EulerICMode::Freestream;
    ic.rho = rho_inf; ic.u = u_inf; ic.v = 0.0; ic.p = p_inf;

    NavierStokesFVMSolver solver(mesh, bcs, gamma, gas_constant, mu, prandtl, cfl, ic, NumericalFluxScheme::Rusanov,
                                  GradientScheme::LeastSquares, 1e-6, 20);

    std::cout << "NS stretched-mesh CFL verification: " << n_x << "x" << n_y << " stretched mesh ("
               << mesh.cells.size() << " cells), first cell height = " << (mesh.cells[0].y_centroid * 2.0)
               << ", mu = " << mu << ", cfl = " << cfl << "\n";

    int last_step = 0;
    for (; last_step < nsteps; ++last_step) {
        solver.step();
        const EulerResidualNorms& r = solver.residual();
        if (std::isnan(r.rho_u) || std::isinf(r.rho_u)) {
            std::cout << "Diverged (NaN/Inf residual) at step " << last_step + 1 << "\n";
            std::cout << "FAIL: run must complete without diverging at this cfl on an anisotropic mesh.\n";
            return false;
        }
    }
    std::cout << "Simulation completed across " << last_step << " steps.\n";
    std::cout << "PASS: run must complete without diverging at this cfl on an anisotropic mesh.\n";
    return true;
}

// Runs Phase 1 of the RANS (Spalart-Allmaras) plan's own verification gate:
// checks compute_wall_distance() (WallDistance.h) against an exact analytic
// case, per docs/archive/rans-spalart-allmaras-tracker.md Phase 1's "planned
// verification". Reuses run_verify_couette()'s sheared mesh purely for its
// geometry -- no solver is involved here, since the wall-distance module has
// no dependency on NavierStokesFVMSolver/RANSFVMSolver, only on mesh
// geometry and a list of wall face indices.
//
// The bottom wall ("bottom" patch) is a straight horizontal line at y=0 for
// every x (the shear x = X + shear*Y, y = Y leaves Y=0's row exactly at y=0
// regardless of shear), so the exact distance from any cell centroid to it
// is simply that centroid's y_centroid -- no approximation, an exact check
// like the NS tracker's Phase 0/1 linear-field gradient tests.
//
// Input:  none
// Output: prints the max error against the exact distance
// Returns: true if the max error is below a tight tolerance (1e-9); false otherwise
bool run_verify_wall_distance() {
    const int n_x = 1, n_y = 16;
    const double shear = 0.3;

    UnstructuredMesh mesh;
    build_sheared_verification_mesh(mesh, n_x, n_y, shear);
    mesh.compute_geometry();

    const int PATCH_BOTTOM = 0;
    std::vector<int> wall_faces;
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        if (mesh.faces[i].patch_id == PATCH_BOTTOM) wall_faces.push_back((int)i);
    }

    std::vector<double> distance = compute_wall_distance(mesh, wall_faces);

    double max_error = 0.0;
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        double exact = mesh.cells[c].y_centroid;
        max_error = std::max(max_error, std::abs(distance[c] - exact));
    }

    std::cout << "Wall-distance verification: " << n_x << "x" << n_y << " sheared mesh (" << mesh.cells.size()
               << " cells), " << wall_faces.size() << " bottom-wall face(s)\n"
               << "  max |distance - y_centroid| = " << max_error << "\n";

    bool passed = max_error < 1e-9;
    std::cout << (passed ? "PASS" : "FAIL")
               << ": distance to a flat horizontal wall at y=0 must equal y_centroid exactly.\n";
    return passed;
}

// Runs Phase 2 of the RANS (Spalart-Allmaras) plan's own verification gate:
// checks compute_sa_source_terms() (SpalartAllmaras.h) in isolation, per
// docs/archive/rans-spalart-allmaras-tracker.md Phase 2's "planned verification" --
// independent of whether these terms are yet wired into any transport
// equation or coupled to a real flow field (that's Phase 3's job).
//
// Test 1: the zero state is a genuine fixed point. On the wall-distance
// verification's 1x16 sheared mesh (reusing compute_wall_distance() from
// Phase 1), with a uniform (synthetic) zero vorticity field and nut
// initialized to 0 everywhere, every source term must come out to exactly
// 0 -- and since nothing changes, an explicit-Euler time loop of nut using
// those source terms must leave nut at exactly 0 after any number of steps.
// This is the SA-tracker analogue of docs/navier-stokes-tracker.md Phase 3's
// "a uniform state stays uniform" check.
//
// Test 2: a manufactured (not physically self-consistent, per the tracker's
// own wording) smooth vorticity field and a smooth, always-positive nut
// field on the same mesh, with grad(nut) reconstructed via the existing
// GradientCalculator (Phase 0 of the NS tracker) and wall distance via
// compute_wall_distance() (Phase 1 of this tracker). Chosen so S~ stays
// comfortably positive everywhere (see SpalartAllmaras.h's note on the
// unimplemented negative-S~ fix) -- checks every cell's production/
// destruction/cross_diffusion is finite and non-negative, i.e. the model's
// nonlinear closure (fv1/fv2/r/g/fw) doesn't blow up or invert sign on a
// smooth field before ever being run through a real coupled solve.
//
// Input:  none
// Output: prints both tests' results
// Returns: true if both tests pass; false otherwise
bool run_verify_sa_source_terms() {
    const int n_x = 1, n_y = 16;
    const double shear = 0.3;
    const double nu = 0.02; // molecular kinematic viscosity, same order as run_verify_couette()'s mu

    UnstructuredMesh mesh;
    build_sheared_verification_mesh(mesh, n_x, n_y, shear);
    mesh.compute_geometry();

    const int PATCH_BOTTOM = 0;
    std::vector<int> wall_faces;
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        if (mesh.faces[i].patch_id == PATCH_BOTTOM) wall_faces.push_back((int)i);
    }
    std::vector<double> wall_distance = compute_wall_distance(mesh, wall_faces);

    // Test 1: zero state stays exactly zero under time marching.
    std::vector<double> nut(mesh.cells.size(), 0.0);
    const int nsteps = 1000;
    const double dt = 0.01;
    bool test1_zero_source = true;
    for (int step = 0; step < nsteps; ++step) {
        for (size_t c = 0; c < mesh.cells.size(); ++c) {
            SASourceTerms s = compute_sa_source_terms(nut[c], nu, /*omega=*/0.0, wall_distance[c], Gradient2{});
            if (s.production != 0.0 || s.destruction != 0.0 || s.cross_diffusion != 0.0) {
                test1_zero_source = false;
            }
            nut[c] += dt * (s.production - s.destruction + s.cross_diffusion);
        }
    }
    bool test1_stayed_zero = std::all_of(nut.begin(), nut.end(), [](double v) { return v == 0.0; });

    std::cout << "SA source-term verification: " << n_x << "x" << n_y << " sheared mesh (" << mesh.cells.size()
               << " cells)\n"
               << "  Test 1 (zero vorticity, nut=0): all source terms exactly 0 every step = "
               << (test1_zero_source ? "yes" : "no") << "; nut after " << nsteps << " steps = "
               << (test1_stayed_zero ? "exactly 0" : "nonzero") << "\n";
    bool test1_passed = test1_zero_source && test1_stayed_zero;
    std::cout << (test1_passed ? "PASS" : "FAIL")
               << ": nut=0 with zero vorticity must be an exact fixed point of production/destruction/diffusion.\n";

    // Test 2: manufactured smooth vorticity + nut fields; check finite, sane, non-negative source terms.
    std::vector<double> omega(mesh.cells.size());
    std::vector<double> nut2(mesh.cells.size());
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        double y = mesh.cells[c].y_centroid;
        omega[c] = 1.0 + y;         // smooth, always-positive synthetic vorticity, O(1)
        nut2[c] = 0.01 * (1.0 + y); // smooth, always-positive synthetic nut, same order as nu
    }
    std::vector<double> boundary_nut(mesh.faces.size(), 0.0);
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
        if (mesh.faces[i].cell_right == -1) boundary_nut[i] = 0.01 * (1.0 + mesh.faces[i].y_mid);
    }
    GradientCalculator grad_calc(mesh, GradientScheme::LeastSquares);
    std::vector<Gradient2> grad_nut = grad_calc.compute(mesh, nut2, boundary_nut);

    bool test2_passed = true;
    double max_production = 0.0, max_destruction = 0.0, max_cross_diffusion = 0.0;
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        SASourceTerms s = compute_sa_source_terms(nut2[c], nu, omega[c], wall_distance[c], grad_nut[c]);
        bool finite = std::isfinite(s.production) && std::isfinite(s.destruction) && std::isfinite(s.cross_diffusion);
        bool sane_sign = (s.production >= 0.0) && (s.destruction >= 0.0) && (s.cross_diffusion >= 0.0);
        if (!finite || !sane_sign) test2_passed = false;
        max_production = std::max(max_production, s.production);
        max_destruction = std::max(max_destruction, s.destruction);
        max_cross_diffusion = std::max(max_cross_diffusion, s.cross_diffusion);
    }
    std::cout << "  Test 2 (manufactured omega/nut fields): max production = " << max_production
               << ", max destruction = " << max_destruction << ", max cross-diffusion = " << max_cross_diffusion << "\n";
    std::cout << (test2_passed ? "PASS" : "FAIL")
               << ": every cell's source terms must be finite and non-negative on a smooth manufactured field.\n";

    return test1_passed && test2_passed;
}

// Runs Phase 3 of the RANS (Spalart-Allmaras) plan's own verification gate:
// RANSFVMSolver's "does the plumbing hold together" stability check, per
// docs/archive/rans-spalart-allmaras-tracker.md Phase 3's "planned verification" --
// deliberately NOT a log-law accuracy check (that's Phase 4's job). Reuses
// run_verify_couette()'s exact 1x16 sheared-mesh planar Couette setup
// (stationary bottom wall, top wall moving at U, zero-gradient outflow
// left/right) with SA now genuinely turned on (initial/freestream
// nut = 3*nu, a typical literature-recommended freestream value), checking
// only that the run stays finite and settles rather than diverging or
// oscillating persistently -- exactly the caution the tracker's Phase 3
// section calls for, since nu_t can be 10-1000x molecular nu and an
// under-conservative compute_dt() would very plausibly reproduce a
// Couette-flow-Phase-4-style divergence.
//
// Input:  none
// Output: prints the final mean-flow and nut residuals, plus nut's min/max
//         across the mesh (a sane run should show nut departing from its
//         uniform initial value, not exploding or going negative)
// Returns: true if the run completes all steps without NaN/Inf and the
//          final rho_u/nut residuals are no larger than their values partway
//          through the run (i.e. settling, not persistently growing); false otherwise
bool run_verify_rans_stability() {
    const int n_x = 1, n_y = 16;
    const double shear = 0.3;
    const double gamma = 1.4, gas_constant = 1.0, mu = 0.02, prandtl = 0.72, prandtl_t = 0.9, cfl = 0.3;
    const double rho0 = 1.0, p0 = 1.0, U = 0.1;
    const double nu = mu / rho0;
    const double initial_nut = 3.0 * nu; // literature-typical freestream nut/nu
    const int nsteps = 20000;

    UnstructuredMesh mesh;
    build_sheared_verification_mesh(mesh, n_x, n_y, shear);
    mesh.compute_geometry();

    std::vector<RANSBoundaryCondition> bcs(mesh.patches.size());
    bcs[0].ns.type = NSBoundaryType::NoSlipWall; // bottom: stationary
    bcs[1].ns.type = NSBoundaryType::Outflow;    // left
    bcs[2].ns.type = NSBoundaryType::NoSlipWall; // top: moving
    bcs[2].ns.wall_u = U;
    bcs[3].ns.type = NSBoundaryType::Outflow;    // right

    EulerInitialCondition ic;
    ic.mode = EulerICMode::Freestream;
    ic.rho = rho0; ic.u = 0.0; ic.v = 0.0; ic.p = p0;

    RANSFVMSolver solver(mesh, bcs, gamma, gas_constant, mu, prandtl, prandtl_t, cfl, ic, initial_nut,
                          NumericalFluxScheme::Rusanov, GradientScheme::LeastSquares, 1e-6, 20);

    double residual_at_half = 0.0, nut_residual_at_half = 0.0;
    int last_step = 0;
    for (; last_step < nsteps; ++last_step) {
        solver.step();
        const EulerResidualNorms& r = solver.residual();
        double nr = solver.nut_residual();
        if (std::isnan(r.rho_u) || std::isinf(r.rho_u) || std::isnan(nr) || std::isinf(nr)) {
            std::cout << "Diverged (NaN/Inf residual) at step " << last_step + 1 << "\n";
            return false;
        }
        if (last_step + 1 == nsteps / 2) {
            residual_at_half = r.rho_u;
            nut_residual_at_half = nr;
        }
    }
    std::cout << "Simulation completed across " << last_step << " steps.\n";

    const std::vector<double>& nut_field = solver.nut_field();
    double min_nut = *std::min_element(nut_field.begin(), nut_field.end());
    double max_nut = *std::max_element(nut_field.begin(), nut_field.end());

    const EulerResidualNorms& r = solver.residual();
    double nr = solver.nut_residual();
    std::cout << "RANS stability verification: " << n_x << "x" << n_y << " sheared mesh (" << mesh.cells.size()
               << " cells), U = " << U << ", mu = " << mu << ", initial nut = " << initial_nut << ", after "
               << nsteps << " steps\n"
               << "  rho_u residual: " << residual_at_half << " (step " << nsteps / 2 << ") -> " << r.rho_u
               << " (final)\n"
               << "  nut residual:   " << nut_residual_at_half << " (step " << nsteps / 2 << ") -> " << nr
               << " (final)\n"
               << "  nut range across mesh: [" << min_nut << ", " << max_nut << "]\n";

    // A tight relative comparison between residual_at_half and the final
    // residual is not meaningful once both have decayed into floating-point
    // noise (~1e-16-1e-17 here) -- an absolute "did it actually reach a
    // settled state" threshold is what "no persistent oscillation" means in
    // practice, matching run_verify_couette()'s use of a fixed pass threshold.
    bool settling = (r.rho_u < 1e-6) && (nr < 1e-6);
    std::cout << (settling ? "PASS" : "FAIL")
               << ": both residuals must decay to a settled state (not diverge or persistently oscillate).\n";
    return settling;
}

// Runs the scalar diffusion solver, driving its step loop directly so a
// residual history, periodic VTK snapshots, stopping criteria, and
// checkpointing can all be handled along the way.
//
// Methodology: resolve the case file's per-patch boundary conditions onto
// the mesh (by patch name), construct UnstructuredFVMSolver, auto-resume
// from case_input.checkpoint_file if it already exists (via
// UnstructuredFVMSolver::set_field()), then advance the solver one step at a
// time up to the absolute target case_input.nsteps. Every step:
//   - if the residual has gone NaN/Inf, the run has diverged: stop
//     immediately (skipping all of the below), write a final VTK snapshot
//     for post-mortem inspection, and report failure -- no checkpoint.
//   - every residual_interval-th step (if case_input.residual_file is set),
//     append the step's residual norm as a CSV row.
//   - every write_interval-th step (if set), write a numbered VTK snapshot.
//   - if case_input.residual_tolerance is set and satisfied, or the user
//     pressed Ctrl+C, stop after this step (a "natural" stop).
// After the loop (by any natural means -- nsteps reached, converged, or
// interrupted), a checkpoint is written (if requested and at least one step
// ran) before the final output_file write, matching the pre-existing
// (untracked) behavior exactly when none of these keys are set.
//
// Input:
//   case_input - run parameters, diffusion fields (alpha/dt/initial_value/
//                initial_radius/boundary_conditions), monitoring fields
//                (residual_file/residual_interval/write_interval), stopping
//                criteria (residual_tolerance), checkpoint_file, and
//                output_file path
//   mesh       - mesh to solve on; its patches are mutated in place with the
//                resolved boundary condition type/value (see UnstructuredMesh.h)
// Output: writes case_input.output_file, case_input.checkpoint_file, and
//         case_input.residual_file/numbered snapshots as requested; prints a
//         warning to stderr per unmatched patch and a message identifying
//         why the run stopped
// Returns: true if the run completed (naturally, converged, or interrupted)
//          and its output was written successfully; false if the run
//          diverged or the output/checkpoint could not be written
bool run_diffusion(const CaseInput& case_input, UnstructuredMesh& mesh) {
    ensure_parent_directory(case_input.output_file);
    ensure_parent_directory(case_input.checkpoint_file);
    ensure_parent_directory(case_input.residual_file);

    // Apply boundary conditions from the case file to the mesh patches, matched by
    // name; a patch with no matching "boundary" entry keeps its default
    // (Dirichlet, 0.0) and a warning is printed.
    for (auto& patch : mesh.patches) {
        bool matched = false;
        for (const auto& bc : case_input.boundary_conditions) {
            if (bc.patch_name == patch.name) {
                patch.type = bc.type;
                patch.value = bc.value;
                matched = true;
                break;
            }
        }
        if (!matched) {
            std::cerr << "Warning: no boundary condition specified for patch '" << patch.name
                       << "', defaulting to dirichlet 0.0\n";
        }
    }

    UnstructuredFVMSolver solver(mesh, case_input.alpha, case_input.dt,
                                  case_input.initial_value, case_input.initial_radius);

    // Auto-resume: if a checkpoint from a previous invocation of this same
    // case file already exists, pick up where it left off instead of
    // starting from the initial condition.
    bool checkpointing = !case_input.checkpoint_file.empty();
    long long resume_start = 1;
    if (checkpointing && Checkpoint::exists(case_input.checkpoint_file)) {
        long long resumed_step = 0;
        unsigned long long resumed_build = 0;
        std::vector<double> resumed_field;
        if (!Checkpoint::read(case_input.checkpoint_file, CheckpointEquation::Diffusion,
                               mesh.cells.size(), resumed_step, resumed_build, resumed_field)) {
            return false;
        }
        solver.set_field(resumed_field);
        resume_start = resumed_step + 1;
        std::cout << "Resuming from checkpoint at step " << resumed_step
                   << " (written by build " << resumed_build << ")\n";
    }

    bool tracking_residual = !case_input.residual_file.empty();
    std::ofstream residual_out;
    if (tracking_residual) {
        residual_out.open(case_input.residual_file);
        residual_out << std::setprecision(case_input.output_precision);
        residual_out << "step,residual\n";
    }

    bool checking_convergence = case_input.residual_tolerance >= 0.0;
    bool diverged = false, converged = false, interrupted = false;
    long long last_completed_step = resume_start - 1;

    for (long long step_index = resume_start; step_index <= case_input.nsteps; ++step_index) {
        solver.step();
        last_completed_step = step_index;

        double residual = solver.residual_norm();
        if (std::isnan(residual) || std::isinf(residual)) {
            diverged = true;
            break;
        }

        if (tracking_residual && step_index % case_input.residual_interval == 0) {
            residual_out << step_index << "," << residual << "\n";
        }
        if (case_input.write_interval > 0 && step_index % case_input.write_interval == 0) {
            VtkWriter::write(numbered_filename(case_input.output_file, step_index), mesh, solver.field(),
                              case_input.output_precision);
        }

        if (checking_convergence && residual < case_input.residual_tolerance) {
            converged = true;
        }
        if (g_interrupt_requested) {
            interrupted = true;
        }
        if (converged || interrupted) {
            break;
        }
    }

    if (diverged) {
        std::cerr << "Error: solver diverged (NaN/Inf residual) at step " << last_completed_step << "\n";
        VtkWriter::write(case_input.output_file, mesh, solver.field(), case_input.output_precision);
        return false;
    }

    if (converged) {
        std::cerr << "Residual converged at step " << last_completed_step << "\n";
    } else if (interrupted) {
        std::cerr << "Interrupted by user at step " << last_completed_step << "\n";
    } else {
        std::cout << "Simulation completed across " << last_completed_step << " steps.\n";
    }

    bool ran_any_step = last_completed_step >= resume_start;
    if (checkpointing && ran_any_step) {
        if (!Checkpoint::write(case_input.checkpoint_file, CheckpointEquation::Diffusion,
                                last_completed_step, FV_BUILD_NUMBER, solver.field())) {
            std::cerr << "Warning: failed to write checkpoint file '" << case_input.checkpoint_file << "'\n";
        }
    }

    if (!VtkWriter::write(case_input.output_file, mesh, solver.field(), case_input.output_precision)) {
        std::cerr << "Failed to write output file: " << case_input.output_file << "\n";
        return false;
    }
    return true;
}

// Runs the advection-diffusion solver for a passive scalar. Identical
// structure to run_diffusion() above (same boundary-condition matching,
// checkpointing, residual tracking/convergence, and stopping-criteria
// handling -- see that function's comment for the full contract), just
// constructing AdvectionDiffusionFVMSolver with its extra u_adv/v_adv/
// gradient_scheme parameters and tagging checkpoints with
// CheckpointEquation::AdvectionDiffusion instead of ::Diffusion.
//
// Input/Output/Returns: same contract as run_diffusion().
bool run_advection_diffusion(const CaseInput& case_input, UnstructuredMesh& mesh) {
    ensure_parent_directory(case_input.output_file);
    ensure_parent_directory(case_input.checkpoint_file);
    ensure_parent_directory(case_input.residual_file);

    for (auto& patch : mesh.patches) {
        bool matched = false;
        for (const auto& bc : case_input.boundary_conditions) {
            if (bc.patch_name == patch.name) {
                patch.type = bc.type;
                patch.value = bc.value;
                matched = true;
                break;
            }
        }
        if (!matched) {
            std::cerr << "Warning: no boundary condition specified for patch '" << patch.name
                       << "', defaulting to dirichlet 0.0\n";
        }
    }

    AdvectionDiffusionFVMSolver solver(mesh, case_input.alpha, case_input.u_adv, case_input.v_adv, case_input.dt,
                                        case_input.gradient_scheme, case_input.initial_value,
                                        case_input.initial_radius);

    bool checkpointing = !case_input.checkpoint_file.empty();
    long long resume_start = 1;
    if (checkpointing && Checkpoint::exists(case_input.checkpoint_file)) {
        long long resumed_step = 0;
        unsigned long long resumed_build = 0;
        std::vector<double> resumed_field;
        if (!Checkpoint::read(case_input.checkpoint_file, CheckpointEquation::AdvectionDiffusion,
                               mesh.cells.size(), resumed_step, resumed_build, resumed_field)) {
            return false;
        }
        solver.set_field(resumed_field);
        resume_start = resumed_step + 1;
        std::cout << "Resuming from checkpoint at step " << resumed_step
                   << " (written by build " << resumed_build << ")\n";
    }

    bool tracking_residual = !case_input.residual_file.empty();
    std::ofstream residual_out;
    if (tracking_residual) {
        residual_out.open(case_input.residual_file);
        residual_out << std::setprecision(case_input.output_precision);
        residual_out << "step,residual\n";
    }

    bool checking_convergence = case_input.residual_tolerance >= 0.0;
    bool diverged = false, converged = false, interrupted = false;
    long long last_completed_step = resume_start - 1;

    for (long long step_index = resume_start; step_index <= case_input.nsteps; ++step_index) {
        solver.step();
        last_completed_step = step_index;

        double residual = solver.residual_norm();
        if (std::isnan(residual) || std::isinf(residual)) {
            diverged = true;
            break;
        }

        if (tracking_residual && step_index % case_input.residual_interval == 0) {
            residual_out << step_index << "," << residual << "\n";
        }
        if (case_input.write_interval > 0 && step_index % case_input.write_interval == 0) {
            VtkWriter::write(numbered_filename(case_input.output_file, step_index), mesh, solver.field(),
                              case_input.output_precision);
        }

        if (checking_convergence && residual < case_input.residual_tolerance) {
            converged = true;
        }
        if (g_interrupt_requested) {
            interrupted = true;
        }
        if (converged || interrupted) {
            break;
        }
    }

    if (diverged) {
        std::cerr << "Error: solver diverged (NaN/Inf residual) at step " << last_completed_step << "\n";
        VtkWriter::write(case_input.output_file, mesh, solver.field(), case_input.output_precision);
        return false;
    }

    if (converged) {
        std::cerr << "Residual converged at step " << last_completed_step << "\n";
    } else if (interrupted) {
        std::cerr << "Interrupted by user at step " << last_completed_step << "\n";
    } else {
        std::cout << "Simulation completed across " << last_completed_step << " steps.\n";
    }

    bool ran_any_step = last_completed_step >= resume_start;
    if (checkpointing && ran_any_step) {
        if (!Checkpoint::write(case_input.checkpoint_file, CheckpointEquation::AdvectionDiffusion,
                                last_completed_step, FV_BUILD_NUMBER, solver.field())) {
            std::cerr << "Warning: failed to write checkpoint file '" << case_input.checkpoint_file << "'\n";
        }
    }

    if (!VtkWriter::write(case_input.output_file, mesh, solver.field(), case_input.output_precision)) {
        std::cerr << "Failed to write output file: " << case_input.output_file << "\n";
        return false;
    }
    return true;
}

// Derives primitive fields (rho, u, v, p) and Mach number from an Euler
// solver's conserved state field and writes them to a VTK file. Shared by
// run_euler()'s periodic snapshots and its final output_file write, so the
// per-cell field-derivation loop isn't duplicated between the two.
//
// Input:
//   filename  - path to write the .vtk file to
//   mesh      - mesh to write geometry from
//   U         - conserved state field, one entry per mesh cell (see EulerState.h)
//   gamma     - ratio of specific heats, dimensionless
//   precision - significant digits for every written double; valid range 1-17
// Output: none (writes to disk)
// Returns: true if the file was written successfully; false otherwise
bool write_euler_fields(const std::string& filename, const UnstructuredMesh& mesh,
                         const std::vector<EulerState>& U, double gamma, int precision) {
    // Recover primitive variables (rho: density; u,v: velocity components;
    // p: pressure -- all in mesh-consistent units, see EulerState.h) plus the
    // dimensionless Mach number, one value per cell.
    NamedField rho{"rho", {}}, u{"u", {}}, v{"v", {}}, p{"p", {}}, mach{"mach", {}};
    for (const auto& state : U) {
        double vel_u = state.rho_u / state.rho;
        double vel_v = state.rho_v / state.rho;
        double pressure_value = pressure(state, gamma);
        rho.values.push_back(state.rho);
        u.values.push_back(vel_u);
        v.values.push_back(vel_v);
        p.values.push_back(pressure_value);
        mach.values.push_back(std::sqrt(vel_u * vel_u + vel_v * vel_v) / sound_speed(state, gamma));
    }

    return VtkWriter::write(filename, mesh, std::vector<NamedField>{rho, u, v, p, mach}, precision);
}

// Runs the compressible Euler solver, driving its step loop directly so a
// residual history, periodic VTK snapshots, stopping criteria, and
// checkpointing can all be handled along the way.
//
// Methodology: resolve the case file's per-patch Euler boundary conditions
// (by patch name) into a vector indexed like mesh.patches, converting any
// Farfield spec's primitive state to conserved form once via from_primitive()
// (see EulerState.h); construct EulerFVMSolver, auto-resume from
// case_input.checkpoint_file if it already exists (via
// EulerFVMSolver::set_field()), then advance the solver one adaptive-dt step
// at a time up to the absolute target case_input.nsteps. Every step:
//   - if any conserved variable's residual has gone NaN/Inf, the run has
//     diverged: stop immediately (skipping all of the below), write a final
//     VTK snapshot for post-mortem inspection, and report failure -- no
//     checkpoint.
//   - every residual_interval-th step (if case_input.residual_file is set),
//     append the step's per-conserved-variable residual norms as a CSV row.
//   - every write_interval-th step (if set), write a numbered VTK snapshot.
//   - if converged (every residual_tolerance_* the user set is satisfied
//     simultaneously) or the user pressed Ctrl+C, stop after this step.
// After the loop (by any natural means -- nsteps reached, converged, or
// interrupted), a checkpoint is written (if requested and at least one step
// ran) before the final output_file write, matching the pre-existing
// (untracked) behavior exactly when none of these keys are set.
//
// Input:
//   case_input - run parameters, Euler fields (gamma/cfl/flux_scheme/
//                exact_riemann_tol/exact_riemann_max_iter/euler_ic/
//                euler_boundary_conditions), monitoring fields
//                (residual_file/residual_interval/write_interval), stopping
//                criteria (residual_tolerance_rho/_rho_u/_rho_v/_E),
//                checkpoint_file, and output_file path
//   mesh       - mesh to solve on (read-only here; unlike run_diffusion,
//                boundary conditions are kept in a separate vector rather
//                than written onto mesh.patches, since EulerBoundaryType has
//                no home in the shared, diffusion-oriented BoundaryPatch struct)
// Output: writes case_input.output_file, case_input.checkpoint_file, and
//         case_input.residual_file/numbered snapshots as requested; prints a
//         warning to stderr per unmatched patch and a message identifying
//         why the run stopped
// Returns: true if the run completed (naturally, converged, or interrupted)
//          and its output was written successfully; false if the run
//          diverged or the output/checkpoint could not be written
bool run_euler(const CaseInput& case_input, const UnstructuredMesh& mesh) {
    ensure_parent_directory(case_input.output_file);
    ensure_parent_directory(case_input.checkpoint_file);
    ensure_parent_directory(case_input.residual_file);

    // Resolve boundary conditions from the case file onto the mesh patches,
    // matched by name; a patch with no matching "boundary" entry defaults to
    // Wall, the conservative choice (it cannot leak or inject mass/energy
    // through an unconfigured boundary), and a warning is printed.
    std::vector<EulerBoundaryCondition> bcs(mesh.patches.size());
    for (size_t i = 0; i < mesh.patches.size(); ++i) {
        bool matched = false;
        for (const auto& bc : case_input.euler_boundary_conditions) {
            if (bc.patch_name == mesh.patches[i].name) {
                bcs[i].type = bc.type;
                if (bc.type == EulerBoundaryType::Farfield) {
                    bcs[i].farfield_state = from_primitive(bc.rho, bc.u, bc.v, bc.p, case_input.gamma);
                }
                matched = true;
                break;
            }
        }
        if (!matched) {
            std::cerr << "Warning: no boundary condition specified for patch '" << mesh.patches[i].name
                       << "', defaulting to wall\n";
        }
    }

    EulerFVMSolver solver(mesh, bcs, case_input.gamma, case_input.cfl, case_input.euler_ic, case_input.flux_scheme,
                           case_input.exact_riemann_tol, case_input.exact_riemann_max_iter);

    // Auto-resume: if a checkpoint from a previous invocation of this same
    // case file already exists, pick up where it left off instead of
    // starting from the initial condition.
    bool checkpointing = !case_input.checkpoint_file.empty();
    long long resume_start = 1;
    if (checkpointing && Checkpoint::exists(case_input.checkpoint_file)) {
        long long resumed_step = 0;
        unsigned long long resumed_build = 0;
        std::vector<EulerState> resumed_field;
        if (!Checkpoint::read(case_input.checkpoint_file, CheckpointEquation::Euler,
                               mesh.cells.size(), resumed_step, resumed_build, resumed_field)) {
            return false;
        }
        solver.set_field(resumed_field);
        resume_start = resumed_step + 1;
        std::cout << "Resuming from checkpoint at step " << resumed_step
                   << " (written by build " << resumed_build << ")\n";
    }

    bool tracking_residual = !case_input.residual_file.empty();
    std::ofstream residual_out;
    if (tracking_residual) {
        residual_out.open(case_input.residual_file);
        residual_out << std::setprecision(case_input.output_precision);
        residual_out << "step,residual_rho,residual_rho_u,residual_rho_v,residual_E\n";
    }

    bool check_rho = case_input.residual_tolerance_rho >= 0.0;
    bool check_rho_u = case_input.residual_tolerance_rho_u >= 0.0;
    bool check_rho_v = case_input.residual_tolerance_rho_v >= 0.0;
    bool check_E = case_input.residual_tolerance_E >= 0.0;
    bool checking_convergence = check_rho || check_rho_u || check_rho_v || check_E;

    bool diverged = false, converged = false, interrupted = false;
    long long last_completed_step = resume_start - 1;

    for (long long step_index = resume_start; step_index <= case_input.nsteps; ++step_index) {
        solver.step();
        last_completed_step = step_index;

        const EulerResidualNorms& r = solver.residual();
        if (std::isnan(r.rho) || std::isinf(r.rho) || std::isnan(r.rho_u) || std::isinf(r.rho_u) ||
            std::isnan(r.rho_v) || std::isinf(r.rho_v) || std::isnan(r.E) || std::isinf(r.E)) {
            diverged = true;
            break;
        }

        if (tracking_residual && step_index % case_input.residual_interval == 0) {
            residual_out << step_index << "," << r.rho << "," << r.rho_u << ","
                         << r.rho_v << "," << r.E << "\n";
        }
        if (case_input.write_interval > 0 && step_index % case_input.write_interval == 0) {
            write_euler_fields(numbered_filename(case_input.output_file, step_index), mesh,
                                solver.field(), case_input.gamma, case_input.output_precision);
        }

        if (checking_convergence) {
            bool rho_ok = !check_rho || r.rho < case_input.residual_tolerance_rho;
            bool rho_u_ok = !check_rho_u || r.rho_u < case_input.residual_tolerance_rho_u;
            bool rho_v_ok = !check_rho_v || r.rho_v < case_input.residual_tolerance_rho_v;
            bool E_ok = !check_E || r.E < case_input.residual_tolerance_E;
            converged = rho_ok && rho_u_ok && rho_v_ok && E_ok;
        }
        if (g_interrupt_requested) {
            interrupted = true;
        }
        if (converged || interrupted) {
            break;
        }
    }

    if (diverged) {
        std::cerr << "Error: solver diverged (NaN/Inf residual) at step " << last_completed_step << "\n";
        write_euler_fields(case_input.output_file, mesh, solver.field(), case_input.gamma, case_input.output_precision);
        return false;
    }

    if (converged) {
        std::cerr << "Residual converged at step " << last_completed_step << "\n";
    } else if (interrupted) {
        std::cerr << "Interrupted by user at step " << last_completed_step << "\n";
    } else {
        std::cout << "Simulation completed across " << last_completed_step << " steps.\n";
    }

    bool ran_any_step = last_completed_step >= resume_start;
    if (checkpointing && ran_any_step) {
        if (!Checkpoint::write(case_input.checkpoint_file, CheckpointEquation::Euler,
                                last_completed_step, FV_BUILD_NUMBER, solver.field())) {
            std::cerr << "Warning: failed to write checkpoint file '" << case_input.checkpoint_file << "'\n";
        }
    }

    if (!write_euler_fields(case_input.output_file, mesh, solver.field(), case_input.gamma, case_input.output_precision)) {
        std::cerr << "Failed to write output file: " << case_input.output_file << "\n";
        return false;
    }
    return true;
}

// Same derivation as write_euler_fields(), plus a "T" temperature field
// (see EulerState.h's temperature()) since that's directly relevant to a
// viscous/heat-conducting solver's results in a way it isn't for pure Euler.
//
// Input:  filename, mesh, U, gamma, precision - same as write_euler_fields()
//         gas_constant - specific gas constant R in p = rho*R*T
// Output/Returns: same contract as write_euler_fields()
bool write_navier_stokes_fields(const std::string& filename, const UnstructuredMesh& mesh,
                                  const std::vector<EulerState>& U, double gamma, double gas_constant,
                                  int precision) {
    NamedField rho{"rho", {}}, u{"u", {}}, v{"v", {}}, p{"p", {}}, mach{"mach", {}}, T{"T", {}};
    for (const auto& state : U) {
        double vel_u = state.rho_u / state.rho;
        double vel_v = state.rho_v / state.rho;
        double pressure_value = pressure(state, gamma);
        rho.values.push_back(state.rho);
        u.values.push_back(vel_u);
        v.values.push_back(vel_v);
        p.values.push_back(pressure_value);
        mach.values.push_back(std::sqrt(vel_u * vel_u + vel_v * vel_v) / sound_speed(state, gamma));
        T.values.push_back(temperature(state, gamma, gas_constant));
    }

    return VtkWriter::write(filename, mesh, std::vector<NamedField>{rho, u, v, p, mach, T}, precision);
}

// Runs the compressible Navier-Stokes solver. Identical structure to
// run_euler() above (same checkpointing, residual tracking/convergence, and
// stopping-criteria handling, reusing the same residual_tolerance_rho/_rho_u/
// _rho_v/_E fields -- see that function's comment for the full contract),
// just constructing NavierStokesFVMSolver with its extra gas_constant/mu/
// prandtl/gradient_scheme parameters, matching per-patch boundary conditions
// from case_input.ns_boundary_conditions (defaulting an unmatched patch to
// NoSlipWall, adiabatic -- the conservative choice for a viscous solver, same
// spirit as Euler defaulting to Wall), and tagging checkpoints with
// CheckpointEquation::NavierStokes instead of ::Euler.
//
// Input/Output/Returns: same contract as run_euler().
bool run_navier_stokes(const CaseInput& case_input, const UnstructuredMesh& mesh) {
    ensure_parent_directory(case_input.output_file);
    ensure_parent_directory(case_input.checkpoint_file);
    ensure_parent_directory(case_input.residual_file);
    ensure_parent_directory(case_input.resolution_report_file);
    ensure_parent_directory(case_input.wall_forces_file);
    ensure_parent_directory(case_input.wall_profile_file);

    std::vector<NSBoundaryCondition> bcs(mesh.patches.size());
    for (size_t i = 0; i < mesh.patches.size(); ++i) {
        bool matched = false;
        for (const auto& bc : case_input.ns_boundary_conditions) {
            if (bc.patch_name == mesh.patches[i].name) {
                bcs[i].type = bc.type;
                bcs[i].wall_u = bc.wall_u;
                bcs[i].wall_v = bc.wall_v;
                bcs[i].is_isothermal_wall = bc.is_isothermal_wall;
                bcs[i].wall_temperature = bc.wall_temperature;
                if (bc.type == NSBoundaryType::Farfield) {
                    bcs[i].farfield_state = from_primitive(bc.rho, bc.u, bc.v, bc.p, case_input.gamma);
                }
                matched = true;
                break;
            }
        }
        if (!matched) {
            std::cerr << "Warning: no boundary condition specified for patch '" << mesh.patches[i].name
                       << "', defaulting to ns_wall (adiabatic no-slip)\n";
        }
    }

    NavierStokesFVMSolver solver(mesh, bcs, case_input.gamma, case_input.gas_constant, case_input.mu,
                                  case_input.prandtl, case_input.cfl, case_input.ns_ic, case_input.flux_scheme,
                                  case_input.gradient_scheme, case_input.exact_riemann_tol,
                                  case_input.exact_riemann_max_iter);

    bool checkpointing = !case_input.checkpoint_file.empty();
    long long resume_start = 1;
    if (checkpointing && Checkpoint::exists(case_input.checkpoint_file)) {
        long long resumed_step = 0;
        unsigned long long resumed_build = 0;
        std::vector<EulerState> resumed_field;
        if (!Checkpoint::read(case_input.checkpoint_file, CheckpointEquation::NavierStokes,
                               mesh.cells.size(), resumed_step, resumed_build, resumed_field)) {
            return false;
        }
        solver.set_field(resumed_field);
        resume_start = resumed_step + 1;
        std::cout << "Resuming from checkpoint at step " << resumed_step
                   << " (written by build " << resumed_build << ")\n";
    }

    bool tracking_residual = !case_input.residual_file.empty();
    std::ofstream residual_out;
    if (tracking_residual) {
        residual_out.open(case_input.residual_file);
        residual_out << std::setprecision(case_input.output_precision);
        residual_out << "step,residual_rho,residual_rho_u,residual_rho_v,residual_E\n";
    }

    bool tracking_resolution = !case_input.resolution_report_file.empty();
    std::ofstream resolution_out;
    if (tracking_resolution) {
        resolution_out.open(case_input.resolution_report_file);
        resolution_out << std::setprecision(case_input.output_precision);
        resolution_out << "step,min_h_over_eta,max_h_over_eta,mean_h_over_eta,n_active_cells\n";
    }

    bool tracking_wall_forces = !case_input.wall_forces_file.empty();
    bool tracking_wall_profile = !case_input.wall_profile_file.empty();
    std::vector<int> wall_faces;
    WallReferenceQuantities wall_ref;
    std::ofstream wall_forces_out;
    if (tracking_wall_forces || tracking_wall_profile) {
        wall_faces = collect_wall_faces(mesh, bcs);
        wall_ref = build_wall_reference_quantities(case_input);
    }
    if (tracking_wall_forces) {
        wall_forces_out.open(case_input.wall_forces_file);
        wall_forces_out << std::setprecision(case_input.output_precision);
        wall_forces_out << "step,patch,friction_drag,pressure_drag,total_drag,cd_friction,cd_pressure,cd_total,"
                            "lift,cl,moment,cm\n";
    }

    // Boundary-layer edge velocity is taken as the reference velocity's
    // magnitude -- physically the same "freestream" state in the common case
    // of one farfield patch (see CaseInput's reference_velocity_* auto-fallback).
    // max_cells_per_march has no case-file key (not part of this plan's
    // scope) -- 200 is a generous cap for any boundary-layer-clustered mesh
    // this project targets; compute_boundary_layer_profiles() reports
    // n_cells_marched if a march ever hits it before plateauing.
    double wall_bl_u_edge = std::hypot(wall_ref.velocity_x_ref, wall_ref.velocity_y_ref);
    const int wall_bl_max_cells_per_march = 200;
    double wall_bl_max_distance = (case_input.boundary_layer_max_distance > 0.0)
                                       ? case_input.boundary_layer_max_distance
                                       : mesh_bounding_box_diagonal(mesh);
    auto compute_wall_node_snapshot = [&]() {
        std::vector<WallFaceSample> samples = solver.compute_wall_traction_samples(wall_faces);
        std::vector<BoundaryLayerProfile> profiles =
            (case_input.boundary_layer_method == BoundaryLayerMethod::PointLocation)
                ? solver.compute_boundary_layer_profile_samples_point_location(
                      wall_faces, wall_bl_u_edge, wall_bl_max_distance, case_input.boundary_layer_n_samples)
                : solver.compute_boundary_layer_profile_samples(wall_faces, wall_bl_u_edge,
                                                                  wall_bl_max_cells_per_march);
        return average_wall_samples_to_nodes(mesh, samples, wall_ref, profiles);
    };

    bool check_rho = case_input.residual_tolerance_rho >= 0.0;
    bool check_rho_u = case_input.residual_tolerance_rho_u >= 0.0;
    bool check_rho_v = case_input.residual_tolerance_rho_v >= 0.0;
    bool check_E = case_input.residual_tolerance_E >= 0.0;
    bool checking_convergence = check_rho || check_rho_u || check_rho_v || check_E;

    bool diverged = false, converged = false, interrupted = false;
    long long last_completed_step = resume_start - 1;

    for (long long step_index = resume_start; step_index <= case_input.nsteps; ++step_index) {
        solver.step();
        last_completed_step = step_index;

        const EulerResidualNorms& r = solver.residual();
        if (std::isnan(r.rho) || std::isinf(r.rho) || std::isnan(r.rho_u) || std::isinf(r.rho_u) ||
            std::isnan(r.rho_v) || std::isinf(r.rho_v) || std::isnan(r.E) || std::isinf(r.E)) {
            diverged = true;
            break;
        }

        if (tracking_residual && step_index % case_input.residual_interval == 0) {
            residual_out << step_index << "," << r.rho << "," << r.rho_u << ","
                         << r.rho_v << "," << r.E << "\n";
        }
        if (tracking_resolution && step_index % case_input.resolution_report_interval == 0) {
            ResolutionDiagnostics diag = solver.compute_resolution_diagnostics();
            resolution_out << step_index << "," << diag.min_ratio << "," << diag.max_ratio << ","
                            << diag.mean_ratio << "," << diag.n_active << "\n";
        }
        if (tracking_wall_forces && step_index % case_input.wall_forces_interval == 0) {
            std::vector<WallFaceSample> samples = solver.compute_wall_traction_samples(wall_faces);
            std::vector<WallForceReport> reports = compute_wall_forces(mesh, samples, wall_ref);
            write_wall_forces_rows(wall_forces_out, step_index, mesh, reports);
        }
        if (tracking_wall_profile && case_input.wall_profile_interval > 0 &&
            step_index % case_input.wall_profile_interval == 0) {
            write_wall_profile_snapshot(numbered_filename(case_input.wall_profile_file, (int)step_index),
                                          case_input.output_precision, mesh, compute_wall_node_snapshot());
        }
        if (case_input.write_interval > 0 && step_index % case_input.write_interval == 0) {
            write_navier_stokes_fields(numbered_filename(case_input.output_file, step_index), mesh, solver.field(),
                                        case_input.gamma, case_input.gas_constant, case_input.output_precision);
        }

        if (checking_convergence) {
            bool rho_ok = !check_rho || r.rho < case_input.residual_tolerance_rho;
            bool rho_u_ok = !check_rho_u || r.rho_u < case_input.residual_tolerance_rho_u;
            bool rho_v_ok = !check_rho_v || r.rho_v < case_input.residual_tolerance_rho_v;
            bool E_ok = !check_E || r.E < case_input.residual_tolerance_E;
            converged = rho_ok && rho_u_ok && rho_v_ok && E_ok;
        }
        if (g_interrupt_requested) {
            interrupted = true;
        }
        if (converged || interrupted) {
            break;
        }
    }

    if (diverged) {
        std::cerr << "Error: solver diverged (NaN/Inf residual) at step " << last_completed_step << "\n";
        write_navier_stokes_fields(case_input.output_file, mesh, solver.field(), case_input.gamma,
                                    case_input.gas_constant, case_input.output_precision);
        if (tracking_wall_profile) {
            write_wall_profile_snapshot(case_input.wall_profile_file, case_input.output_precision, mesh,
                                          compute_wall_node_snapshot());
        }
        return false;
    }

    if (converged) {
        std::cerr << "Residual converged at step " << last_completed_step << "\n";
    } else if (interrupted) {
        std::cerr << "Interrupted by user at step " << last_completed_step << "\n";
    } else {
        std::cout << "Simulation completed across " << last_completed_step << " steps.\n";
    }

    bool ran_any_step = last_completed_step >= resume_start;
    if (checkpointing && ran_any_step) {
        if (!Checkpoint::write(case_input.checkpoint_file, CheckpointEquation::NavierStokes,
                                last_completed_step, FV_BUILD_NUMBER, solver.field())) {
            std::cerr << "Warning: failed to write checkpoint file '" << case_input.checkpoint_file << "'\n";
        }
    }

    if (!write_navier_stokes_fields(case_input.output_file, mesh, solver.field(), case_input.gamma,
                                     case_input.gas_constant, case_input.output_precision)) {
        std::cerr << "Failed to write output file: " << case_input.output_file << "\n";
        return false;
    }
    if (tracking_wall_profile) {
        write_wall_profile_snapshot(case_input.wall_profile_file, case_input.output_precision, mesh,
                                      compute_wall_node_snapshot());
    }
    return true;
}

// Same derivation as write_navier_stokes_fields(), plus "nut" (the transported
// SA scalar) and "nu_t" (the derived turbulent eddy viscosity, see
// sa_eddy_viscosity() in SpalartAllmaras.h) -- RANS's one addition over the
// mean-flow fields Navier-Stokes already writes.
//
// Input:  filename, mesh, U, gamma, gas_constant, precision - same as write_navier_stokes_fields()
//         nut          - the transported nu-tilde field, one per mesh.cells entry
//         mu           - molecular dynamic viscosity, needed to derive nu_t from nut
//         sa_constants - SA model constants used to derive nu_t (see SpalartAllmaras.h)
// Output/Returns: same contract as write_euler_fields()
bool write_rans_fields(const std::string& filename, const UnstructuredMesh& mesh, const std::vector<EulerState>& U,
                         const std::vector<double>& nut, double gamma, double gas_constant, double mu,
                         const SAModelConstants& sa_constants, int precision) {
    NamedField rho{"rho", {}}, u{"u", {}}, v{"v", {}}, p{"p", {}}, mach{"mach", {}}, T{"T", {}};
    NamedField nut_field{"nut", {}}, nu_t_field{"nu_t", {}};
    for (size_t c = 0; c < U.size(); ++c) {
        const EulerState& state = U[c];
        double vel_u = state.rho_u / state.rho;
        double vel_v = state.rho_v / state.rho;
        double pressure_value = pressure(state, gamma);
        rho.values.push_back(state.rho);
        u.values.push_back(vel_u);
        v.values.push_back(vel_v);
        p.values.push_back(pressure_value);
        mach.values.push_back(std::sqrt(vel_u * vel_u + vel_v * vel_v) / sound_speed(state, gamma));
        T.values.push_back(temperature(state, gamma, gas_constant));
        nut_field.values.push_back(nut[c]);
        nu_t_field.values.push_back(sa_eddy_viscosity(nut[c], mu / state.rho, sa_constants));
    }

    return VtkWriter::write(filename, mesh, std::vector<NamedField>{rho, u, v, p, mach, T, nut_field, nu_t_field},
                              precision);
}

// Runs the RANS (Spalart-Allmaras) solver. Identical structure to
// run_navier_stokes() above (same checkpointing, residual tracking/
// convergence, wall-diagnostics wiring, and stopping-criteria handling),
// plus: matching per-patch boundary conditions from
// case_input.rans_boundary_conditions (defaulting an unmatched patch to
// rans_wall, adiabatic -- same conservative default as Navier-Stokes),
// constructing RANSFVMSolver with its extra prandtl_t/initial_nut/sa_constants
// parameters, tracking the nut residual alongside rho/rho_u/rho_v/E, and
// tagging checkpoints with CheckpointEquation::RANS (which round-trips both
// U and nut, unlike ::NavierStokes's U-only payload).
//
// Input/Output/Returns: same contract as run_navier_stokes().
bool run_rans(const CaseInput& case_input, const UnstructuredMesh& mesh) {
    ensure_parent_directory(case_input.output_file);
    ensure_parent_directory(case_input.checkpoint_file);
    ensure_parent_directory(case_input.residual_file);
    ensure_parent_directory(case_input.wall_forces_file);
    ensure_parent_directory(case_input.wall_profile_file);

    std::vector<RANSBoundaryCondition> bcs(mesh.patches.size());
    for (size_t i = 0; i < mesh.patches.size(); ++i) {
        bool matched = false;
        for (const auto& bc : case_input.rans_boundary_conditions) {
            if (bc.patch_name == mesh.patches[i].name) {
                bcs[i].ns.type = bc.type;
                bcs[i].ns.wall_u = bc.wall_u;
                bcs[i].ns.wall_v = bc.wall_v;
                bcs[i].ns.is_isothermal_wall = bc.is_isothermal_wall;
                bcs[i].ns.wall_temperature = bc.wall_temperature;
                if (bc.type == NSBoundaryType::Farfield) {
                    bcs[i].ns.farfield_state = from_primitive(bc.rho, bc.u, bc.v, bc.p, case_input.gamma);
                    bcs[i].farfield_nut = bc.farfield_nut;
                }
                matched = true;
                break;
            }
        }
        if (!matched) {
            std::cerr << "Warning: no boundary condition specified for patch '" << mesh.patches[i].name
                       << "', defaulting to rans_wall (adiabatic no-slip)\n";
        }
    }

    RANSFVMSolver solver(mesh, bcs, case_input.gamma, case_input.gas_constant, case_input.mu, case_input.prandtl,
                          case_input.prandtl_t, case_input.cfl, case_input.rans_ic, case_input.initial_nut,
                          case_input.flux_scheme, case_input.gradient_scheme, case_input.exact_riemann_tol,
                          case_input.exact_riemann_max_iter, case_input.sa_constants);

    bool checkpointing = !case_input.checkpoint_file.empty();
    long long resume_start = 1;
    if (checkpointing && Checkpoint::exists(case_input.checkpoint_file)) {
        long long resumed_step = 0;
        unsigned long long resumed_build = 0;
        std::vector<EulerState> resumed_field;
        std::vector<double> resumed_nut;
        if (!Checkpoint::read(case_input.checkpoint_file, CheckpointEquation::RANS, mesh.cells.size(), resumed_step,
                               resumed_build, resumed_field, resumed_nut)) {
            return false;
        }
        solver.set_field(resumed_field, resumed_nut);
        resume_start = resumed_step + 1;
        std::cout << "Resuming from checkpoint at step " << resumed_step
                   << " (written by build " << resumed_build << ")\n";
    }

    bool tracking_residual = !case_input.residual_file.empty();
    std::ofstream residual_out;
    if (tracking_residual) {
        residual_out.open(case_input.residual_file);
        residual_out << std::setprecision(case_input.output_precision);
        residual_out << "step,residual_rho,residual_rho_u,residual_rho_v,residual_E,residual_nut\n";
    }

    bool tracking_wall_forces = !case_input.wall_forces_file.empty();
    bool tracking_wall_profile = !case_input.wall_profile_file.empty();
    std::vector<int> wall_faces;
    WallReferenceQuantities wall_ref;
    std::ofstream wall_forces_out;
    if (tracking_wall_forces || tracking_wall_profile) {
        wall_faces = collect_wall_faces(mesh, bcs);
        wall_ref = build_wall_reference_quantities(case_input);
    }
    if (tracking_wall_forces) {
        wall_forces_out.open(case_input.wall_forces_file);
        wall_forces_out << std::setprecision(case_input.output_precision);
        wall_forces_out << "step,patch,friction_drag,pressure_drag,total_drag,cd_friction,cd_pressure,cd_total,"
                            "lift,cl,moment,cm\n";
    }

    double wall_bl_u_edge = std::hypot(wall_ref.velocity_x_ref, wall_ref.velocity_y_ref);
    const int wall_bl_max_cells_per_march = 200;
    double wall_bl_max_distance = (case_input.boundary_layer_max_distance > 0.0)
                                       ? case_input.boundary_layer_max_distance
                                       : mesh_bounding_box_diagonal(mesh);
    auto compute_wall_node_snapshot = [&]() {
        std::vector<WallFaceSample> samples = solver.compute_wall_traction_samples(wall_faces);
        std::vector<BoundaryLayerProfile> profiles =
            (case_input.boundary_layer_method == BoundaryLayerMethod::PointLocation)
                ? solver.compute_boundary_layer_profile_samples_point_location(
                      wall_faces, wall_bl_u_edge, wall_bl_max_distance, case_input.boundary_layer_n_samples)
                : solver.compute_boundary_layer_profile_samples(wall_faces, wall_bl_u_edge,
                                                                  wall_bl_max_cells_per_march);
        return average_wall_samples_to_nodes(mesh, samples, wall_ref, profiles);
    };

    bool check_rho = case_input.residual_tolerance_rho >= 0.0;
    bool check_rho_u = case_input.residual_tolerance_rho_u >= 0.0;
    bool check_rho_v = case_input.residual_tolerance_rho_v >= 0.0;
    bool check_E = case_input.residual_tolerance_E >= 0.0;
    bool check_nut = case_input.residual_tolerance_nut >= 0.0;
    bool checking_convergence = check_rho || check_rho_u || check_rho_v || check_E || check_nut;

    bool diverged = false, converged = false, interrupted = false;
    long long last_completed_step = resume_start - 1;

    for (long long step_index = resume_start; step_index <= case_input.nsteps; ++step_index) {
        solver.step();
        last_completed_step = step_index;

        const EulerResidualNorms& r = solver.residual();
        double nr = solver.nut_residual();
        if (std::isnan(r.rho) || std::isinf(r.rho) || std::isnan(r.rho_u) || std::isinf(r.rho_u) ||
            std::isnan(r.rho_v) || std::isinf(r.rho_v) || std::isnan(r.E) || std::isinf(r.E) ||
            std::isnan(nr) || std::isinf(nr)) {
            diverged = true;
            break;
        }

        if (tracking_residual && step_index % case_input.residual_interval == 0) {
            residual_out << step_index << "," << r.rho << "," << r.rho_u << ","
                         << r.rho_v << "," << r.E << "," << nr << "\n";
        }
        if (tracking_wall_forces && step_index % case_input.wall_forces_interval == 0) {
            std::vector<WallFaceSample> samples = solver.compute_wall_traction_samples(wall_faces);
            std::vector<WallForceReport> reports = compute_wall_forces(mesh, samples, wall_ref);
            write_wall_forces_rows(wall_forces_out, step_index, mesh, reports);
        }
        if (tracking_wall_profile && case_input.wall_profile_interval > 0 &&
            step_index % case_input.wall_profile_interval == 0) {
            write_wall_profile_snapshot(numbered_filename(case_input.wall_profile_file, (int)step_index),
                                          case_input.output_precision, mesh, compute_wall_node_snapshot());
        }
        if (case_input.write_interval > 0 && step_index % case_input.write_interval == 0) {
            write_rans_fields(numbered_filename(case_input.output_file, step_index), mesh, solver.field(),
                                solver.nut_field(), case_input.gamma, case_input.gas_constant, case_input.mu,
                                case_input.sa_constants, case_input.output_precision);
        }

        if (checking_convergence) {
            bool rho_ok = !check_rho || r.rho < case_input.residual_tolerance_rho;
            bool rho_u_ok = !check_rho_u || r.rho_u < case_input.residual_tolerance_rho_u;
            bool rho_v_ok = !check_rho_v || r.rho_v < case_input.residual_tolerance_rho_v;
            bool E_ok = !check_E || r.E < case_input.residual_tolerance_E;
            bool nut_ok = !check_nut || nr < case_input.residual_tolerance_nut;
            converged = rho_ok && rho_u_ok && rho_v_ok && E_ok && nut_ok;
        }
        if (g_interrupt_requested) {
            interrupted = true;
        }
        if (converged || interrupted) {
            break;
        }
    }

    if (diverged) {
        std::cerr << "Error: solver diverged (NaN/Inf residual) at step " << last_completed_step << "\n";
        write_rans_fields(case_input.output_file, mesh, solver.field(), solver.nut_field(), case_input.gamma,
                            case_input.gas_constant, case_input.mu, case_input.sa_constants,
                            case_input.output_precision);
        if (tracking_wall_profile) {
            write_wall_profile_snapshot(case_input.wall_profile_file, case_input.output_precision, mesh,
                                          compute_wall_node_snapshot());
        }
        return false;
    }

    if (converged) {
        std::cerr << "Residual converged at step " << last_completed_step << "\n";
    } else if (interrupted) {
        std::cerr << "Interrupted by user at step " << last_completed_step << "\n";
    } else {
        std::cout << "Simulation completed across " << last_completed_step << " steps.\n";
    }

    bool ran_any_step = last_completed_step >= resume_start;
    if (checkpointing && ran_any_step) {
        if (!Checkpoint::write(case_input.checkpoint_file, CheckpointEquation::RANS, last_completed_step,
                                FV_BUILD_NUMBER, solver.field(), solver.nut_field())) {
            std::cerr << "Warning: failed to write checkpoint file '" << case_input.checkpoint_file << "'\n";
        }
    }

    if (!write_rans_fields(case_input.output_file, mesh, solver.field(), solver.nut_field(), case_input.gamma,
                             case_input.gas_constant, case_input.mu, case_input.sa_constants,
                             case_input.output_precision)) {
        std::cerr << "Failed to write output file: " << case_input.output_file << "\n";
        return false;
    }
    if (tracking_wall_profile) {
        write_wall_profile_snapshot(case_input.wall_profile_file, case_input.output_precision, mesh,
                                      compute_wall_node_snapshot());
    }
    return true;
}

} // namespace

// Entry point: standalone 2D unstructured finite-volume solver, supporting
// scalar diffusion, advection-diffusion, the compressible Euler equations,
// the compressible Navier-Stokes equations, or RANS (Spalart-Allmaras)
// (selected by the case file's "equation" key).
//
// Usage: FiniteVolume.exe <case_file>
//        FiniteVolume.exe --version
//
// Pipeline: parse the case file (run parameters + boundary conditions) ->
// parse the Gmsh mesh it references -> assign boundary conditions to the
// mesh's named patches -> run the explicit solver, stopping at nsteps, on
// residual convergence, on a Ctrl+C interrupt, or (immediately) on
// divergence -> write the resulting field(s) to a VTK file. If
// case_input.checkpoint_file already exists, the run auto-resumes from it
// instead of starting over; see run_diffusion()/run_euler() for the full
// stopping/checkpoint contract.
//
// Input:   argv[1] - path to a case file (see CaseInput.h for its format);
//          or "--version" to print identifying build info and exit;
//          or "--validate-mesh <mesh_file>" to parse a .msh/.fvmesh file
//          (via MeshReader::read) and print its node/cell/face/patch counts
//          and min/max cell volume, without running any solver;
//          or "--verify-gradient" to check GradientCalculator's Green-Gauss
//          and Least-Squares schemes against a known linear analytic field
//          on a self-generated skewed mesh (see run_verify_gradient());
//          or "--verify-advdiff" to check AdvectionDiffusionFVMSolver
//          against the known steady 1D advection-diffusion profile on the
//          same kind of mesh (see run_verify_advection_diffusion());
//          or "--verify-ns-uniform" to check that NavierStokesFVMSolver
//          leaves a uniform freestream state exactly unperturbed with
//          viscosity turned on (see run_verify_navier_stokes_uniform());
//          or "--verify-couette" to check NavierStokesFVMSolver's steady
//          shear profile against the analytic planar Couette flow solution
//          (see run_verify_couette());
//          or "--verify-wall-distance" to check compute_wall_distance()
//          against an exact analytic flat-wall distance (see
//          run_verify_wall_distance());
//          or "--verify-wall-forces" to check WallTraction.h's
//          tau_wall/Cf/Cp/y+/Cm against closed-form planar Couette flow
//          values (see run_verify_wall_forces());
//          or "--verify-bl-marching-unstructured" to stress-test
//          compute_boundary_layer_profiles()'s marching heuristic against a
//          genuinely unstructured (checkerboard-triangulated) near-wall
//          mesh (see run_verify_bl_marching_unstructured());
//          or "--verify-bl-point-location" to check that
//          compute_boundary_layer_profiles_point_location() does not share
//          that failure mode on the same mesh (see
//          run_verify_bl_point_location());
//          or "--verify-sa-source" to check compute_sa_source_terms() in
//          isolation against a zero-vorticity fixed point and a manufactured
//          smooth field (see run_verify_sa_source_terms());
//          or "--verify-rans-stability" to check RANSFVMSolver runs stably
//          (no divergence, residuals settle) on a simple sheared/turbulent
//          setup (see run_verify_rans_stability());
//          or "--verify-flat-plate" to check RANSFVMSolver's flat-plate
//          case against a laminar reference profile at a tractable
//          (sub-transition) Reynolds number -- see
//          run_verify_flat_plate_boundary_layer() for why this isn't the
//          turbulent log-law Phase 4 originally set out to reach
// Output:  writes the result to the .vtk path given by the case file's
//          output_file key (and its checkpoint_file/residual_file, if set);
//          diagnostics/errors go to stderr
// Returns: 0 if the run completed (naturally, converged, or interrupted) and
//          its output was written, --version was requested, --validate-mesh
//          succeeded, or one of the --verify-* checks passed; 1 if the case
//          file, mesh file, or output file could not be loaded/written, the
//          run diverged, or one of the --verify-* checks failed
int main(int argc, char** argv) {
#ifdef _WIN32
    // Force this process's console output to UTF-8, regardless of the
    // codepage the console session started in, so accented characters
    // (e.g. FV_PUBLISHER/FV_COPYRIGHT in --version's output) always display
    // correctly rather than depending on the user having run "chcp 65001".
    SetConsoleOutputCP(CP_UTF8);
#endif
    std::signal(SIGINT, handle_sigint);

    if (argc >= 2 && std::string(argv[1]) == "--version") {
        std::cout << "FiniteVolume " << FV_VERSION_FULL << "\n"
                   << FV_DESCRIPTION << "\n"
                   << FV_COPYRIGHT << "\n";
        return 0;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-gradient") {
        return run_verify_gradient() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-advdiff") {
        return run_verify_advection_diffusion() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-ns-uniform") {
        return run_verify_navier_stokes_uniform() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-couette") {
        return run_verify_couette() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-ns-stretched-cfl") {
        return run_verify_ns_stretched_cfl() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-wall-forces") {
        return run_verify_wall_forces() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-bl-marching-unstructured") {
        return run_verify_bl_marching_unstructured() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-bl-point-location") {
        return run_verify_bl_point_location() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-wall-distance") {
        return run_verify_wall_distance() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-sa-source") {
        return run_verify_sa_source_terms() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-rans-stability") {
        return run_verify_rans_stability() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--verify-flat-plate") {
        return run_verify_flat_plate_boundary_layer() ? 0 : 1;
    }

    if (argc >= 2 && std::string(argv[1]) == "--validate-mesh") {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " --validate-mesh <mesh_file>\n";
            return 1;
        }

        UnstructuredMesh mesh;
        if (!MeshReader::read(argv[2], mesh)) {
            std::cerr << "Failed to load mesh file: " << argv[2] << "\n";
            return 1;
        }

        double min_volume = mesh.cells.empty() ? 0.0 : mesh.cells[0].volume;
        double max_volume = min_volume;
        for (const auto& cell : mesh.cells) {
            if (cell.volume < min_volume) min_volume = cell.volume;
            if (cell.volume > max_volume) max_volume = cell.volume;
        }

        std::cout << "Mesh OK: " << argv[2] << "\n"
                   << "  nodes:   " << mesh.nodes.size() << "\n"
                   << "  cells:   " << mesh.cells.size() << "\n"
                   << "  faces:   " << mesh.faces.size() << "\n"
                   << "  patches: " << mesh.patches.size() << "\n";
        for (const auto& patch : mesh.patches) {
            std::cout << "    - " << patch.name << "\n";
        }
        std::cout << "  cell volume: min=" << min_volume << " max=" << max_volume << "\n";

        return 0;
    }

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <case_file>\n";
        return 1;
    }

    CaseInput case_input;
    if (!case_input.load(argv[1])) {
        std::cerr << "Failed to load case file: " << argv[1] << "\n";
        return 1;
    }

    // Rebase relative output paths under scratch_dir, if set; mesh_file is
    // never affected (see resolve_output_path).
    case_input.output_file = resolve_output_path(case_input.output_file, case_input.scratch_dir);
    case_input.checkpoint_file = resolve_output_path(case_input.checkpoint_file, case_input.scratch_dir);
    case_input.residual_file = resolve_output_path(case_input.residual_file, case_input.scratch_dir);
    case_input.resolution_report_file = resolve_output_path(case_input.resolution_report_file, case_input.scratch_dir);
    case_input.wall_forces_file = resolve_output_path(case_input.wall_forces_file, case_input.scratch_dir);
    case_input.wall_profile_file = resolve_output_path(case_input.wall_profile_file, case_input.scratch_dir);

#ifdef _OPENMP
    if (case_input.num_threads > 0) {
        omp_set_num_threads(case_input.num_threads);
    }
#endif

    UnstructuredMesh mesh;
    if (!MeshReader::read(case_input.mesh_file, mesh)) {
        std::cerr << "Failed to load mesh file: " << case_input.mesh_file << "\n";
        return 1;
    }
    // Every solver constructor computes face geometry (x_mid/y_mid/area/nx/ny)
    // on its OWN private mesh copy (mesh_with_geometry()-style pattern), never
    // on this outer 'mesh' -- harmless for every existing consumer (VtkWriter
    // and cell.volume/centroid-based code only need geometry MeshReader
    // already fills in during parsing), but run_navier_stokes()'s
    // compute_wall_forces() call needs valid face.area/nx/ny on THIS mesh
    // directly. Idempotent, so calling it here doesn't change anything for
    // equation sets that never needed it.
    mesh.compute_geometry();

    bool ok;
    if (case_input.equation == EquationSet::Diffusion) {
        ok = run_diffusion(case_input, mesh);
    } else if (case_input.equation == EquationSet::AdvectionDiffusion) {
        ok = run_advection_diffusion(case_input, mesh);
    } else if (case_input.equation == EquationSet::NavierStokes) {
        ok = run_navier_stokes(case_input, mesh);
    } else if (case_input.equation == EquationSet::RANS) {
        ok = run_rans(case_input, mesh);
    } else {
        ok = run_euler(case_input, mesh);
    }

    if (!ok) {
        return 1;
    }

    std::cout << "Results written to " << case_input.output_file << "\n";
    return 0;
}
