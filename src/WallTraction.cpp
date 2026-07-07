// SPDX-License-Identifier: GPL-3.0-only
#include "WallTraction.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace {

// Rotates (nx, ny) by +90 degrees: (x, y) -> (-y, x).
std::pair<double, double> rotate90(double nx, double ny) {
    return {-ny, nx};
}

// Derives a BoundaryLayerProfile from a raw (distance, u_tangential) sample
// list, shared by both compute_boundary_layer_profiles() (cell-marching) and
// compute_boundary_layer_profiles_point_location() -- everything downstream
// of "how the samples were obtained" (the delta_99 crossing interpolation
// and the displacement/momentum trapezoidal integration) is identical
// between the two methods. 'distance'/'u_tangential' must be the same
// length, sorted by increasing distance, with distance[0] == 0.0 (the
// wall's own prescribed tangential speed).
BoundaryLayerProfile derive_boundary_layer_profile(int wall_face_idx, double x_mid, double y_mid,
                                                     const std::vector<double>& distance,
                                                     const std::vector<double>& u_tangential, double u_edge,
                                                     int n_samples) {
    BoundaryLayerProfile profile;
    profile.wall_face_index = wall_face_idx;
    profile.x_mid = x_mid;
    profile.y_mid = y_mid;
    profile.n_cells_marched = n_samples;

    // delta_99: linear interpolation where u/u_edge first crosses 0.99.
    double delta99 = distance.back();
    for (size_t i = 0; i + 1 < distance.size(); ++i) {
        double ratio_i = u_tangential[i] / u_edge, ratio_ip1 = u_tangential[i + 1] / u_edge;
        if (ratio_i < 0.99 && ratio_ip1 >= 0.99) {
            double t = (0.99 - ratio_i) / (ratio_ip1 - ratio_i);
            delta99 = distance[i] + t * (distance[i + 1] - distance[i]);
            break;
        }
    }
    profile.delta_99 = delta99;

    // Trapezoidal integration of (1 - u/u_edge) and (u/u_edge)*(1 - u/u_edge).
    double displacement = 0.0, momentum = 0.0;
    for (size_t i = 0; i + 1 < distance.size(); ++i) {
        double h = distance[i + 1] - distance[i];
        double r0 = u_tangential[i] / u_edge, r1 = u_tangential[i + 1] / u_edge;
        double f0_disp = 1.0 - r0, f1_disp = 1.0 - r1;
        double f0_mom = r0 * (1.0 - r0), f1_mom = r1 * (1.0 - r1);
        displacement += 0.5 * h * (f0_disp + f1_disp);
        momentum += 0.5 * h * (f0_mom + f1_mom);
    }
    profile.displacement_thickness = displacement;
    profile.momentum_thickness = momentum;
    profile.shape_factor = (momentum != 0.0) ? displacement / momentum : 0.0;

    return profile;
}

// Point-in-polygon test (standard even-odd ray-casting) for a simple CCW
// polygon cell -- works for any simple (non-self-intersecting) polygon,
// convex or concave, since this project's cells have no convexity guarantee.
bool point_in_cell(const UnstructuredMesh& mesh, const Cell& cell, double px, double py) {
    bool inside = false;
    int n = (int)cell.node_ids.size();
    for (int i = 0, j = n - 1; i < n; j = i++) {
        const Node& pi = mesh.nodes[cell.node_ids[i]];
        const Node& pj = mesh.nodes[cell.node_ids[j]];
        bool crosses = (pi.y > py) != (pj.y > py);
        if (crosses) {
            double x_intersect = pi.x + (py - pi.y) * (pj.x - pi.x) / (pj.y - pi.y);
            if (px < x_intersect) inside = !inside;
        }
    }
    return inside;
}

// Brute-force O(cells) point-location: returns the index of the first cell
// containing (px, py), or -1 if none does (e.g. the point lies outside the
// mesh domain). No spatial acceleration structure -- consistent with this
// project's existing "simple first, document scaling limits" stance (see
// WallDistance.h's own brute-force point-to-segment search).
int locate_cell(const UnstructuredMesh& mesh, double px, double py) {
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        if (point_in_cell(mesh, mesh.cells[c], px, py)) return (int)c;
    }
    return -1;
}

} // namespace

// See WallTraction.h for the input/output contract and methodology.
std::vector<WallFaceSample> compute_wall_traction(
    const UnstructuredMesh& mesh, const std::vector<int>& wall_faces,
    const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& p, const std::vector<double>& rho,
    const std::vector<double>& effective_viscosity,
    const std::vector<double>& boundary_u, const std::vector<double>& boundary_v,
    const std::vector<Gradient2>& grad_u, const std::vector<Gradient2>& grad_v)
{
    std::vector<WallFaceSample> samples;
    samples.reserve(wall_faces.size());

    for (int face_idx : wall_faces) {
        const Face& face = mesh.faces[face_idx];
        int cl = face.cell_left;

        FaceGradient fu = face_gradient(mesh, face, u[cl], 0.0, grad_u[cl], Gradient2{}, boundary_u[face_idx]);
        FaceGradient fv = face_gradient(mesh, face, v[cl], 0.0, grad_v[cl], Gradient2{}, boundary_v[face_idx]);

        Gradient2 gu = corrected_face_gradient_vector(fu, face.nx, face.ny);
        Gradient2 gv = corrected_face_gradient_vector(fv, face.nx, face.ny);

        double mu = effective_viscosity[cl];
        double div_v = gu.dphidx + gv.dphidy;
        double tau_xx = mu * (2.0 * gu.dphidx - (2.0 / 3.0) * div_v);
        double tau_yy = mu * (2.0 * gv.dphidy - (2.0 / 3.0) * div_v);
        double tau_xy = mu * (gu.dphidy + gv.dphidx);

        double traction_x = tau_xx * face.nx + tau_xy * face.ny;
        double traction_y = tau_xy * face.nx + tau_yy * face.ny;

        auto [tx, ty] = rotate90(face.nx, face.ny);

        WallFaceSample sample;
        sample.face_index = face_idx;
        sample.x_mid = face.x_mid;
        sample.y_mid = face.y_mid;
        sample.p = p[cl];
        sample.rho = rho[cl];
        sample.effective_viscosity = mu;
        sample.tau_wall = -(traction_x * tx + traction_y * ty);
        sample.tx = tx;
        sample.ty = ty;
        sample.y_wall_normal = face_normal_distance(mesh, face);
        samples.push_back(sample);
    }

    return samples;
}

// See WallTraction.h for the input/output contract.
double skin_friction_coefficient(const WallFaceSample& sample, const WallReferenceQuantities& ref) {
    double v2 = ref.velocity_x_ref * ref.velocity_x_ref + ref.velocity_y_ref * ref.velocity_y_ref;
    return sample.tau_wall / (0.5 * ref.rho_ref * v2);
}

// See WallTraction.h for the input/output contract.
double pressure_coefficient(const WallFaceSample& sample, const WallReferenceQuantities& ref) {
    double v2 = ref.velocity_x_ref * ref.velocity_x_ref + ref.velocity_y_ref * ref.velocity_y_ref;
    return (sample.p - ref.p_ref) / (0.5 * ref.rho_ref * v2);
}

// See WallTraction.h for the input/output contract.
double wall_y_plus(const WallFaceSample& sample) {
    double u_tau = std::sqrt(std::fabs(sample.tau_wall) / sample.rho);
    return sample.rho * u_tau * sample.y_wall_normal / sample.effective_viscosity;
}

// See WallTraction.h for the input/output contract and methodology.
std::vector<WallForceReport> compute_wall_forces(const UnstructuredMesh& mesh,
                                                   const std::vector<WallFaceSample>& samples,
                                                   const WallReferenceQuantities& ref)
{
    if (samples.empty()) return {};

    double v_mag = std::hypot(ref.velocity_x_ref, ref.velocity_y_ref);
    double dir_x = ref.velocity_x_ref / v_mag, dir_y = ref.velocity_y_ref / v_mag;
    double lift_dir_x = -dir_y, lift_dir_y = dir_x; // drag direction rotated +90 deg
    double dynamic_pressure_length = 0.5 * ref.rho_ref * v_mag * v_mag * ref.length_ref;
    double dynamic_pressure_length2 = dynamic_pressure_length * ref.length_ref;

    // Accumulate per patch_id, plus a -1 "domain total" bucket, in one pass.
    std::map<int, WallForceReport> by_patch;

    for (const WallFaceSample& s : samples) {
        const Face& face = mesh.faces[s.face_index];
        double friction_x = face.area * s.tau_wall * s.tx;
        double friction_y = face.area * s.tau_wall * s.ty;
        double pressure_x = face.area * s.p * face.nx;
        double pressure_y = face.area * s.p * face.ny;

        double rx = s.x_mid - ref.moment_reference_x;
        double ry = s.y_mid - ref.moment_reference_y;
        double dfx = friction_x + pressure_x, dfy = friction_y + pressure_y;
        double moment_contribution = rx * dfy - ry * dfx;

        for (int key : {face.patch_id, -1}) {
            WallForceReport& r = by_patch[key];
            r.patch_id = key;
            double friction_drag_i = friction_x * dir_x + friction_y * dir_y;
            double pressure_drag_i = pressure_x * dir_x + pressure_y * dir_y;
            r.friction_drag += friction_drag_i;
            r.pressure_drag += pressure_drag_i;
            r.total_drag += friction_drag_i + pressure_drag_i;
            r.lift += (friction_x + pressure_x) * lift_dir_x + (friction_y + pressure_y) * lift_dir_y;
            r.moment += moment_contribution;
        }
    }

    std::vector<WallForceReport> result;
    result.reserve(by_patch.size());
    for (auto& [patch_id, r] : by_patch) {
        r.cd_friction = r.friction_drag / dynamic_pressure_length;
        r.cd_pressure = r.pressure_drag / dynamic_pressure_length;
        r.cd_total = r.total_drag / dynamic_pressure_length;
        r.cl = r.lift / dynamic_pressure_length;
        r.cm = r.moment / dynamic_pressure_length2;
        result.push_back(r);
    }
    // Move the domain total (-1) to the end, matching this module's
    // documented "per patch, then a final domain total" order.
    std::stable_partition(result.begin(), result.end(), [](const WallForceReport& r) { return r.patch_id != -1; });
    return result;
}

// See WallTraction.h for the input/output contract and methodology.
std::vector<BoundaryLayerProfile> compute_boundary_layer_profiles(
    const UnstructuredMesh& mesh, const std::vector<int>& wall_faces,
    const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& boundary_u, const std::vector<double>& boundary_v,
    double u_edge, int max_cells_per_march)
{
    std::vector<BoundaryLayerProfile> profiles;
    profiles.reserve(wall_faces.size());

    for (int wall_face_idx : wall_faces) {
        const Face& wall_face = mesh.faces[wall_face_idx];
        // wall_face.nx/ny is outward FROM the fluid cell (cell_left), i.e. it
        // points OUT of the fluid domain into the wall -- marching into the
        // fluid (away from the wall) needs the opposite direction.
        double march_x = -wall_face.nx, march_y = -wall_face.ny; // fixed marching direction
        auto [tan_x, tan_y] = rotate90(wall_face.nx, wall_face.ny); // fixed tangent, for sampling only

        // Distance-0 sample: the wall's own prescribed tangential velocity.
        // Sampled as a SPEED (fabs), not signed relative to this face's own
        // (arbitrary, per-wall) tangent direction -- e.g. a "top" wall's
        // outward-normal-derived tangent points opposite a "bottom" wall's
        // even when the actual flow moves the same physical direction past
        // both, and delta_99/the thickness integrals below assume u/u_edge
        // is a well-behaved ratio approaching +1, not -1.
        std::vector<double> distance = {0.0};
        std::vector<double> u_tangential = {
            std::fabs(boundary_u[wall_face_idx] * tan_x + boundary_v[wall_face_idx] * tan_y)};

        int current_cell = wall_face.cell_left;
        int entry_face = wall_face_idx;
        double cumulative_distance = 0.0;
        int n_marched = 0;

        while (n_marched < max_cells_per_march) {
            int best_face = -1, best_neighbor = -1;
            double best_alignment = -1.0;

            for (int face_idx : mesh.cells[current_cell].faces) {
                if (face_idx == entry_face) continue;
                const Face& f = mesh.faces[face_idx];
                if (f.cell_right == -1) continue; // never cross a boundary face

                int neighbor = (f.cell_left == current_cell) ? f.cell_right : f.cell_left;
                double outward_nx = (f.cell_left == current_cell) ? f.nx : -f.nx;
                double outward_ny = (f.cell_left == current_cell) ? f.ny : -f.ny;
                double alignment = outward_nx * march_x + outward_ny * march_y;

                if (alignment > best_alignment) {
                    best_alignment = alignment;
                    best_face = face_idx;
                    best_neighbor = neighbor;
                }
            }

            // Require a genuinely forward-pointing face (alignment > 0), not
            // just the "least bad" of a lateral/backward set -- near the top
            // of a column, the only remaining interior faces can be lateral
            // (e.g. the last row's side neighbor, once its own top face is a
            // boundary), and picking one of those anyway would silently walk
            // the march sideways instead of correctly stopping.
            if (best_face == -1 || best_alignment <= 0.0) break;

            cumulative_distance += face_normal_distance(mesh, mesh.faces[best_face]);
            current_cell = best_neighbor;
            entry_face = best_face;
            ++n_marched;

            double u_tan = std::fabs(u[current_cell] * tan_x + v[current_cell] * tan_y);
            distance.push_back(cumulative_distance);
            u_tangential.push_back(u_tan);

            if (u_tan >= 0.99 * std::fabs(u_edge)) break;
        }

        profiles.push_back(derive_boundary_layer_profile(wall_face_idx, wall_face.x_mid, wall_face.y_mid, distance,
                                                           u_tangential, u_edge, n_marched));
    }

    return profiles;
}

// See WallTraction.h for the input/output contract and methodology.
std::vector<BoundaryLayerProfile> compute_boundary_layer_profiles_point_location(
    const UnstructuredMesh& mesh, const std::vector<int>& wall_faces,
    const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& boundary_u, const std::vector<double>& boundary_v,
    double u_edge, double max_distance, int n_samples)
{
    std::vector<BoundaryLayerProfile> profiles;
    profiles.reserve(wall_faces.size());

    double step_size = max_distance / n_samples;

    for (int wall_face_idx : wall_faces) {
        const Face& wall_face = mesh.faces[wall_face_idx];
        // wall_face.nx/ny is outward FROM the fluid cell, i.e. it points OUT
        // of the fluid domain into the wall -- sampling into the fluid needs
        // the opposite (inward) direction.
        double inward_x = -wall_face.nx, inward_y = -wall_face.ny;
        auto [tan_x, tan_y] = rotate90(wall_face.nx, wall_face.ny); // fixed tangent, for sampling only

        // Distance-0 sample: the wall's own prescribed tangential velocity,
        // as an unsigned speed -- see compute_boundary_layer_profiles()'s
        // identical note on why (a "top"-type wall's tangent convention can
        // point opposite a "bottom"-type wall's).
        std::vector<double> distance = {0.0};
        std::vector<double> u_tangential = {
            std::fabs(boundary_u[wall_face_idx] * tan_x + boundary_v[wall_face_idx] * tan_y)};

        int n_taken = 0;
        for (int step = 1; step <= n_samples; ++step) {
            double d = step * step_size;
            double px = wall_face.x_mid + d * inward_x, py = wall_face.y_mid + d * inward_y;
            int cell = locate_cell(mesh, px, py);
            if (cell == -1) break; // walked past the domain boundary

            ++n_taken;
            double u_tan = std::fabs(u[cell] * tan_x + v[cell] * tan_y);
            distance.push_back(d);
            u_tangential.push_back(u_tan);

            if (u_tan >= 0.99 * std::fabs(u_edge)) break;
        }

        profiles.push_back(derive_boundary_layer_profile(wall_face_idx, wall_face.x_mid, wall_face.y_mid, distance,
                                                           u_tangential, u_edge, n_taken));
    }

    return profiles;
}

// See WallTraction.h for the input/output contract.
std::vector<WallNodeSample> average_wall_samples_to_nodes(const UnstructuredMesh& mesh,
                                                            const std::vector<WallFaceSample>& traction,
                                                            const WallReferenceQuantities& ref,
                                                            const std::vector<BoundaryLayerProfile>& profiles)
{
    // Map face_index -> boundary-layer profile, if any were supplied.
    std::map<int, const BoundaryLayerProfile*> profile_by_face;
    for (const auto& bl : profiles) profile_by_face[bl.wall_face_index] = &bl;

    // Group contributions by (node_index, patch_id) -- a corner node shared
    // by two different patches is reported once per patch, never merged.
    struct Accum {
        int patch_id;
        double x, y;
        double sum_tau_wall = 0.0, sum_cf = 0.0, sum_cp = 0.0, sum_yplus = 0.0;
        double sum_delta99 = 0.0, sum_disp = 0.0, sum_mom = 0.0, sum_H = 0.0;
        int count = 0;
    };
    std::map<std::pair<int, int>, Accum> by_node_patch;

    for (const WallFaceSample& s : traction) {
        const Face& face = mesh.faces[s.face_index];
        double cf = skin_friction_coefficient(s, ref);
        double cp = pressure_coefficient(s, ref);
        double yplus = wall_y_plus(s);

        double delta99 = 0.0, disp = 0.0, mom = 0.0, H = 0.0;
        auto it = profile_by_face.find(s.face_index);
        if (it != profile_by_face.end()) {
            delta99 = it->second->delta_99;
            disp = it->second->displacement_thickness;
            mom = it->second->momentum_thickness;
            H = it->second->shape_factor;
        }

        for (int node_idx : {face.node1, face.node2}) {
            std::pair<int, int> key{node_idx, face.patch_id};
            Accum& a = by_node_patch[key];
            a.patch_id = face.patch_id;
            a.x = mesh.nodes[node_idx].x;
            a.y = mesh.nodes[node_idx].y;
            a.sum_tau_wall += s.tau_wall;
            a.sum_cf += cf;
            a.sum_cp += cp;
            a.sum_yplus += yplus;
            a.sum_delta99 += delta99;
            a.sum_disp += disp;
            a.sum_mom += mom;
            a.sum_H += H;
            ++a.count;
        }
    }

    std::vector<WallNodeSample> result;
    result.reserve(by_node_patch.size());
    for (const auto& [key, a] : by_node_patch) {
        WallNodeSample n;
        n.node_index = key.first;
        n.patch_id = a.patch_id;
        n.x = a.x;
        n.y = a.y;
        n.tau_wall = a.sum_tau_wall / a.count;
        n.Cf = a.sum_cf / a.count;
        n.Cp = a.sum_cp / a.count;
        n.y_plus = a.sum_yplus / a.count;
        n.delta_99 = a.sum_delta99 / a.count;
        n.displacement_thickness = a.sum_disp / a.count;
        n.momentum_thickness = a.sum_mom / a.count;
        n.shape_factor = a.sum_H / a.count;
        result.push_back(n);
    }
    return result;
}
