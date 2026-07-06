// SPDX-License-Identifier: GPL-3.0-only
#ifndef WALLDISTANCE_H_INCLUDED
#define WALLDISTANCE_H_INCLUDED

#include <vector>

#include "UnstructuredMesh.h"

// Computes, for every cell in 'mesh', the minimum Euclidean distance from
// that cell's centroid to the nearest point on any face listed in
// 'wall_face_indices' (the nearest point on the face's line SEGMENT, not
// just its midpoint).
//
// This is the wall-distance field the SA turbulence model's destruction
// term needs (see docs/archive/rans-spalart-allmaras-tracker.md Phase 1/2). Kept
// mesh-geometry-only, with no dependency on RANSFVMSolver or NSBoundaryType,
// since it depends only on which faces are walls, not on any flow field --
// same reasoning as GradientCalculator's precomputed Least-Squares matrix
// (see GradientReconstruction.h). A future RANSFVMSolver just needs to hand
// it the indices of whichever faces have NSBoundaryType::NoSlipWall.
//
// Methodology: brute-force O(cells * wall_face_indices.size()) point-to-
// segment distance, no acceleration structure -- consistent with this
// project's existing "simple first, document scaling limits" stance (see
// MeshReader's edge-matching, GradientCalculator's per-cell matrix). Would
// need a spatial acceleration structure (e.g. a tree over wall face
// midpoints) to scale to a mesh with many more wall faces or cells than this
// project currently targets.
//
// Input:
//   mesh              - mesh to compute distances on; must already have
//                        compute_geometry() called (cell centroids and face
//                        node positions populated)
//   wall_face_indices - indices into mesh.faces treated as no-slip walls
// Returns: one distance per mesh cell (same order as mesh.cells), in mesh
//          length units. If wall_face_indices is empty, every entry is
//          +infinity (there is no wall to measure distance to).
std::vector<double> compute_wall_distance(const UnstructuredMesh& mesh, const std::vector<int>& wall_face_indices);

#endif // WALLDISTANCE_H_INCLUDED
