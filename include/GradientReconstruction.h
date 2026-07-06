// SPDX-License-Identifier: GPL-3.0-only
#ifndef GRADIENTRECONSTRUCTION_H_INCLUDED
#define GRADIENTRECONSTRUCTION_H_INCLUDED

#include <array>
#include <vector>

#include "UnstructuredMesh.h"

// Selects how GradientCalculator reconstructs a per-cell gradient from a
// cell-centered scalar field.
enum class GradientScheme {
    GreenGauss,   // Divergence-theorem sum of face-interpolated values over each cell's bounding faces.
    LeastSquares  // Weighted least-squares fit against neighboring cell/boundary-face values.
};

// A single cell's reconstructed 2D gradient (dphi/dx, dphi/dy), in
// (field units) / (mesh length units).
struct Gradient2 {
    double dphidx = 0.0;
    double dphidy = 0.0;
};

// Reconstructs cell-centered gradients of an arbitrary scalar field on an
// arbitrary polygon mesh, using either of two standard finite-volume
// schemes (see GradientScheme above). Geometry-dependent setup
// (Least-Squares' per-cell normal-equation matrix) is precomputed once at
// construction, since it depends only on mesh geometry, not on the field
// being differentiated -- compute() can then be called repeatedly (e.g. once
// per time step) without repeating that setup.
//
// Methodology:
//   - Green-Gauss: grad(phi)_C = (1/V_C) * sum_faces (phi_face * n_outward * A_face),
//     where phi_face is the inverse-distance-weighted average of the two
//     neighboring cell values (or the caller-supplied boundary value for a
//     boundary face). Exact for a linear field only when every face midpoint
//     lies on the line joining its two cell centroids (e.g. Cartesian/orthogonal
//     meshes); has an O(skewness) error on general unstructured meshes.
//   - Least-Squares: fits grad(phi)_C by minimizing, over every neighboring
//     cell centroid and boundary face midpoint, the weighted squared error
//     between (phi_neighbor - phi_C) and grad(phi)_C . (x_neighbor - x_C),
//     with weight 1/distance^2. Exact for a linear field on any
//     non-degenerate mesh (any neighbor point set spanning both
//     dimensions), since the normal equations are constructed to reproduce
//     affine functions exactly.
class GradientCalculator {
public:
    // Input:
    //   mesh   - mesh to reconstruct gradients on; must already have
    //            compute_geometry() called (cell volumes/centroids and face
    //            midpoints/areas/normals populated). Not copied/stored --
    //            the same mesh must be passed to every compute() call on
    //            this object.
    //   scheme - which reconstruction scheme compute() uses
    // Output: for LeastSquares, precomputes and caches each cell's inverse
    //         normal-equation matrix from mesh geometry; for GreenGauss, does
    //         no precomputation (each face's weights are cheap enough to
    //         recompute on the fly during compute(), since they're only
    //         ever used once each).
    GradientCalculator(const UnstructuredMesh& mesh, GradientScheme scheme);

    // Reconstructs the gradient of 'phi' in every cell.
    // Input:
    //   mesh          - the same mesh passed to the constructor (same
    //                   topology/geometry; only field values change between calls)
    //   phi           - cell-centered field, one value per mesh.cells entry
    //   boundary_phi  - mesh.faces-sized array; only entries at boundary
    //                   face indices (cell_right == -1) are read, giving the
    //                   field's value AT that face's midpoint. Entries for
    //                   internal faces are ignored.
    // Returns: one Gradient2 per mesh cell (same order as mesh.cells)
    std::vector<Gradient2> compute(const UnstructuredMesh& mesh, const std::vector<double>& phi,
                                    const std::vector<double>& boundary_phi) const;

private:
    GradientScheme scheme;

    // Per-cell inverse 2x2 weighted normal-equation matrix (symmetric,
    // stored as {inv00, inv01, inv11}), used only when scheme ==
    // LeastSquares. Populated once at construction; independent of the
    // field being differentiated.
    std::vector<std::array<double, 3>> lsq_inv;
};

// A face's non-orthogonality-corrected normal derivative, plus the plain
// interpolated full gradient vector at the face (see face_gradient() below).
struct FaceGradient {
    double dphidn = 0.0; // Corrected normal derivative, outward from face.cell_left, field units / mesh length units
    Gradient2 grad_f;    // Interpolated full gradient vector at the face (distance-weighted average of the two cell gradients)
};

// Computes the non-orthogonality-corrected face-normal derivative of a
// scalar field at 'face', using the standard over-relaxed decomposition
// (Jasak 1996). Averaging the two cells' gradients and dotting with the
// face normal ("grad_f . n") would decouple the flux from (phi_R - phi_L)
// entirely and risks odd-even/checkerboard oscillations; a naive two-point
// difference ((phi_R - phi_L)/dist, e.g. UnstructuredFVMSolver's current
// diffusion flux) is only correct when the face normal is parallel to the
// line joining the two cell centroids (an orthogonal mesh) and has a real
// error otherwise. This function keeps the direct two-point term (so the
// flux stays sensitive to (phi_R - phi_L), avoiding the checkerboard
// failure mode) but rescales it by 1/cosTheta to correct for non-alignment,
// then adds a correction term using the interpolated cell-gradient vector
// for the remainder the direct term doesn't capture:
//
//   n           = face.nx, face.ny (outward from face.cell_left)
//   d_hat       = unit vector from face.cell_left's centroid to the "R"
//                 point (face.cell_right's centroid, or the face midpoint
//                 for a boundary face)
//   cosTheta    = n . d_hat
//   dphidn      = (phi_R - phi_L) / (dist * cosTheta)          -- direct term
//               + grad_f . (n - d_hat / cosTheta)               -- correction term
//
// On an orthogonal mesh, n == d_hat, cosTheta == 1, and the correction
// direction (n - d_hat/cosTheta) is exactly zero, so this reduces to the
// naive two-point formula.
//
// Input:
//   mesh                - mesh 'face' belongs to (must have
//                          compute_geometry() called)
//   face                - the face to evaluate (internal or boundary)
//   phi_L, phi_R        - field values in face.cell_left/cell_right
//                          (phi_R ignored for a boundary face)
//   grad_L, grad_R      - reconstructed cell gradients (see
//                          GradientCalculator) in face.cell_left/cell_right
//                          (grad_R ignored for a boundary face)
//   boundary_phi_value  - prescribed field value at the face midpoint; only
//                          used when face.cell_right == -1
// Returns: the FaceGradient described above; grad_f is grad_L itself for a
//          boundary face (no second cell gradient to blend with)
FaceGradient face_gradient(const UnstructuredMesh& mesh, const Face& face, double phi_L, double phi_R,
                            const Gradient2& grad_L, const Gradient2& grad_R, double boundary_phi_value);

// Reconstructs a face's full 2D gradient vector by taking its plain
// interpolated average gradient (fg.grad_f) and replacing just its normal
// component with face_gradient()'s non-orthogonality-corrected one
// (fg.dphidn): grad_corrected = grad_f + (dphidn - grad_f . n) * n.
//
// face_gradient() only corrects the NORMAL component of a face's gradient --
// that's the component whose accuracy actually matters for avoiding
// checkerboard decoupling (see docs/navier-stokes-tracker.md Phase 1). A
// viscous stress tensor's off-diagonal/tangential terms (NavierStokesFVMSolver,
// RANSFVMSolver, and WallTraction.h) need the full 2D gradient vector, though,
// hence this correction.
//
// Input:  fg     - a face's FaceGradient (see face_gradient() above)
//         nx, ny - the same face normal passed to face_gradient()
// Returns: the corrected full 2D gradient vector at the face
Gradient2 corrected_face_gradient_vector(const FaceGradient& fg, double nx, double ny);

// Face-normal-projected distance between face.cell_left's centroid and the
// opposite point (face.cell_right's centroid, or the face midpoint for a
// boundary face) -- i.e. dist * cosTheta from face_gradient()'s over-relaxed
// decomposition above. This is the correct physical length scale for a
// diffusive/viscous stability estimate across the face, since it's the same
// length scale face_gradient() itself divides by to compute the flux -- a
// stability estimate built on a different length (e.g. cell volume/face
// area) can systematically understate stiffness on anisotropic cells.
//
// Input:  mesh, face - as in face_gradient() above
// Returns: dist * cosTheta (see face_gradient()'s methodology comment)
double face_normal_distance(const UnstructuredMesh& mesh, const Face& face);

#endif // GRADIENTRECONSTRUCTION_H_INCLUDED
