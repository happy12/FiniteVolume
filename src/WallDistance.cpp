// SPDX-License-Identifier: GPL-3.0-only
#include "WallDistance.h"

#include <algorithm>
#include <limits>

namespace {

// Squared distance from point (px, py) to the line segment (x1, y1)-(x2, y2).
// Kept squared (no sqrt) since only the minimum over many segments is
// needed per cell; the final sqrt is applied once, after the minimum over
// all wall faces is found.
double point_segment_distance_sq(double px, double py, double x1, double y1, double x2, double y2) {
    double dx = x2 - x1, dy = y2 - y1;
    double len_sq = dx * dx + dy * dy;
    double t = (len_sq > 0.0) ? ((px - x1) * dx + (py - y1) * dy) / len_sq : 0.0;
    t = std::clamp(t, 0.0, 1.0);
    double ex = px - (x1 + t * dx);
    double ey = py - (y1 + t * dy);
    return ex * ex + ey * ey;
}

} // namespace

// See WallDistance.h for the input/output contract.
std::vector<double> compute_wall_distance(const UnstructuredMesh& mesh, const std::vector<int>& wall_face_indices) {
    std::vector<double> distance(mesh.cells.size(), std::numeric_limits<double>::infinity());
    if (wall_face_indices.empty()) return distance;

    #pragma omp parallel for
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        const Cell& cell = mesh.cells[c];
        double min_dist_sq = std::numeric_limits<double>::infinity();
        for (int face_idx : wall_face_indices) {
            const Face& face = mesh.faces[face_idx];
            const Node& n1 = mesh.nodes[face.node1];
            const Node& n2 = mesh.nodes[face.node2];
            double d2 = point_segment_distance_sq(cell.x_centroid, cell.y_centroid, n1.x, n1.y, n2.x, n2.y);
            min_dist_sq = std::min(min_dist_sq, d2);
        }
        distance[c] = std::sqrt(min_dist_sq);
    }
    return distance;
}
