// SPDX-License-Identifier: GPL-3.0-only
#include "GradientReconstruction.h"

#include <cmath>
#include <utility>

namespace {

// Returns the outward-from-'cell_c' unit normal components of 'face',
// accounting for Face::nx/ny being defined relative to cell_left: negated
// when 'cell_c' is the face's cell_right side.
std::pair<double, double> outward_normal(const Face& face, int cell_c) {
    if (face.cell_left == cell_c) return {face.nx, face.ny};
    return {-face.nx, -face.ny};
}

// Shared geometry behind both face_gradient()'s over-relaxed decomposition
// and face_normal_distance(): the distance and direction from
// face.cell_left's centroid to the "opposite" point (face.cell_right's
// centroid, or the face midpoint for a boundary face), and that direction's
// projection onto the face normal.
struct FaceNormalGeometry {
    double dist;             // centroid-to-opposite-point distance
    double cos_theta;        // dhat . face normal
    double dhat_x, dhat_y;   // unit direction from cell_left's centroid to the opposite point
};

FaceNormalGeometry face_normal_geometry(const UnstructuredMesh& mesh, const Face& face) {
    const Cell& cell_L = mesh.cells[face.cell_left];
    bool is_boundary = (face.cell_right == -1);
    double x_R = is_boundary ? face.x_mid : mesh.cells[face.cell_right].x_centroid;
    double y_R = is_boundary ? face.y_mid : mesh.cells[face.cell_right].y_centroid;

    double dx = x_R - cell_L.x_centroid;
    double dy = y_R - cell_L.y_centroid;
    double dist = std::hypot(dx, dy);
    double dhat_x = dx / dist, dhat_y = dy / dist;
    double cos_theta = face.nx * dhat_x + face.ny * dhat_y;
    return {dist, cos_theta, dhat_x, dhat_y};
}

} // namespace

// See GradientReconstruction.h for the input/output contract.
GradientCalculator::GradientCalculator(const UnstructuredMesh& mesh, GradientScheme scheme)
    : scheme(scheme)
{
    if (scheme != GradientScheme::LeastSquares) return;

    // Precompute each cell's weighted normal-equation matrix
    // M = sum_i w_i * [dx_i^2, dx_i*dy_i; dx_i*dy_i, dy_i^2] (w_i = 1/dist_i^2)
    // over its neighbor points (opposite-cell centroids for internal faces,
    // boundary face midpoints for boundary faces), then invert it once so
    // compute() only needs a 2x2 matrix-vector product per cell.
    lsq_inv.resize(mesh.cells.size());
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        const Cell& cell = mesh.cells[c];
        double m00 = 0.0, m01 = 0.0, m11 = 0.0;
        for (int face_idx : cell.faces) {
            const Face& face = mesh.faces[face_idx];
            double nx, ny;
            if (face.cell_right == -1) {
                nx = face.x_mid;
                ny = face.y_mid;
            } else {
                int neighbor = (face.cell_left == (int)c) ? face.cell_right : face.cell_left;
                nx = mesh.cells[neighbor].x_centroid;
                ny = mesh.cells[neighbor].y_centroid;
            }
            double dx = nx - cell.x_centroid;
            double dy = ny - cell.y_centroid;
            double w = 1.0 / (dx * dx + dy * dy);
            m00 += w * dx * dx;
            m01 += w * dx * dy;
            m11 += w * dy * dy;
        }
        double det = m00 * m11 - m01 * m01;
        lsq_inv[c] = {m11 / det, -m01 / det, m00 / det};
    }
}

// See GradientReconstruction.h for the input/output contract.
std::vector<Gradient2> GradientCalculator::compute(const UnstructuredMesh& mesh, const std::vector<double>& phi,
                                                     const std::vector<double>& boundary_phi) const
{
    std::vector<Gradient2> grad(mesh.cells.size());

    if (scheme == GradientScheme::GreenGauss) {
        #pragma omp parallel for
        for (size_t c = 0; c < mesh.cells.size(); ++c) {
            const Cell& cell = mesh.cells[c];
            double sx = 0.0, sy = 0.0;
            for (int face_idx : cell.faces) {
                const Face& face = mesh.faces[face_idx];
                double phi_face;
                if (face.cell_right == -1) {
                    phi_face = boundary_phi[face_idx];
                } else {
                    int other = (face.cell_left == (int)c) ? face.cell_right : face.cell_left;
                    // Inverse-distance-weighted average between the two cell
                    // centroids' values, weighted toward whichever centroid
                    // is closer to the face midpoint.
                    double d_c = std::hypot(cell.x_centroid - face.x_mid, cell.y_centroid - face.y_mid);
                    double d_o = std::hypot(mesh.cells[other].x_centroid - face.x_mid,
                                             mesh.cells[other].y_centroid - face.y_mid);
                    double w_c = d_o / (d_c + d_o);
                    phi_face = w_c * phi[c] + (1.0 - w_c) * phi[other];
                }
                auto [nx, ny] = outward_normal(face, (int)c);
                sx += phi_face * nx * face.area;
                sy += phi_face * ny * face.area;
            }
            grad[c].dphidx = sx / cell.volume;
            grad[c].dphidy = sy / cell.volume;
        }
        return grad;
    }

    // Least-Squares: b = sum_i w_i * dphi_i * [dx_i, dy_i], then grad = Minv * b,
    // using the matrix precomputed (and already inverted) at construction.
    #pragma omp parallel for
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        const Cell& cell = mesh.cells[c];
        double bx = 0.0, by = 0.0;
        for (int face_idx : cell.faces) {
            const Face& face = mesh.faces[face_idx];
            double nx, ny, phi_n;
            if (face.cell_right == -1) {
                nx = face.x_mid;
                ny = face.y_mid;
                phi_n = boundary_phi[face_idx];
            } else {
                int neighbor = (face.cell_left == (int)c) ? face.cell_right : face.cell_left;
                nx = mesh.cells[neighbor].x_centroid;
                ny = mesh.cells[neighbor].y_centroid;
                phi_n = phi[neighbor];
            }
            double dx = nx - cell.x_centroid;
            double dy = ny - cell.y_centroid;
            double w = 1.0 / (dx * dx + dy * dy);
            double dphi = phi_n - phi[c];
            bx += w * dx * dphi;
            by += w * dy * dphi;
        }
        const std::array<double, 3>& inv = lsq_inv[c]; // {inv00, inv01, inv11}
        grad[c].dphidx = inv[0] * bx + inv[1] * by;
        grad[c].dphidy = inv[1] * bx + inv[2] * by;
    }
    return grad;
}

// See GradientReconstruction.h for the input/output contract and methodology.
FaceGradient face_gradient(const UnstructuredMesh& mesh, const Face& face, double phi_L, double phi_R,
                            const Gradient2& grad_L, const Gradient2& grad_R, double boundary_phi_value)
{
    const Cell& cell_L = mesh.cells[face.cell_left];
    bool is_boundary = (face.cell_right == -1);

    double x_R = is_boundary ? face.x_mid : mesh.cells[face.cell_right].x_centroid;
    double y_R = is_boundary ? face.y_mid : mesh.cells[face.cell_right].y_centroid;
    double phi_R_value = is_boundary ? boundary_phi_value : phi_R;

    FaceNormalGeometry geom = face_normal_geometry(mesh, face);

    double direct_term = (phi_R_value - phi_L) / (geom.dist * geom.cos_theta);

    double kx = face.nx - geom.dhat_x / geom.cos_theta;
    double ky = face.ny - geom.dhat_y / geom.cos_theta;

    Gradient2 grad_f;
    if (is_boundary) {
        grad_f = grad_L; // No second cell gradient to blend with.
    } else {
        double d_L = std::hypot(cell_L.x_centroid - face.x_mid, cell_L.y_centroid - face.y_mid);
        double d_R = std::hypot(x_R - face.x_mid, y_R - face.y_mid);
        double w_L = d_R / (d_L + d_R);
        grad_f.dphidx = w_L * grad_L.dphidx + (1.0 - w_L) * grad_R.dphidx;
        grad_f.dphidy = w_L * grad_L.dphidy + (1.0 - w_L) * grad_R.dphidy;
    }

    FaceGradient result;
    result.dphidn = direct_term + grad_f.dphidx * kx + grad_f.dphidy * ky;
    result.grad_f = grad_f;
    return result;
}

// See GradientReconstruction.h for the input/output contract and methodology.
double face_normal_distance(const UnstructuredMesh& mesh, const Face& face)
{
    FaceNormalGeometry geom = face_normal_geometry(mesh, face);
    return geom.dist * geom.cos_theta;
}

// See GradientReconstruction.h for the input/output contract and methodology.
Gradient2 corrected_face_gradient_vector(const FaceGradient& fg, double nx, double ny)
{
    double avg_normal_component = fg.grad_f.dphidx * nx + fg.grad_f.dphidy * ny;
    double delta = fg.dphidn - avg_normal_component;
    return {fg.grad_f.dphidx + delta * nx, fg.grad_f.dphidy + delta * ny};
}
