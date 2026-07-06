// SPDX-License-Identifier: GPL-3.0-only
#ifndef UNSTRUCTUREDMESH_H_INCLUDED
#define UNSTRUCTUREDMESH_H_INCLUDED

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

/* notes
Cache Friendliness: The face array loop (for (const auto& face : mesh.faces)) forces linear reads through continuous memory vectors.
However, accessing the mesh.cells[face.cell_left] metrics causes indirect pointer lookup (random memory jumps).
To mitigate this in a large production solver, you should sort your cell and face indices using an algorithm
like Reverse Cuthill-McKee (RCM) to group adjacent nodes together in RAM.
*/

// Represents a geometric vertex position.
// Units: same length units as the input mesh file (e.g. metres if the mesh
// was generated in metres) -- the solver itself is unit-agnostic.
struct Node {
    double x, y; // Coordinates in mesh length units
};

// Represents a 1D line segment interface between two 2D Control Volumes.
// Faces are the flux-integration surfaces of the finite volume method: every
// conserved quantity that leaves cell_left through this face enters cell_right
// (or crosses a domain boundary, for boundary faces).
struct Face {
    int node1, node2;       // Local vertex indices (into UnstructuredMesh::nodes) forming the face
    int cell_left;          // Index of Owner/Left Cell (into UnstructuredMesh::cells)
    int cell_right;         // Index of Neighbor/Right Cell (-1 if it's a boundary face)
    int patch_id = -1;      // Index into UnstructuredMesh::patches (-1 = internal face)
    double nx, ny;          // Outward pointing unit normal vector components (dimensionless, relative to cell_left)
    double area;            // Length of the 2D face line segment, in mesh length units
    double x_mid, y_mid;    // Midpoint coordinates of the face, in mesh length units
};

// Represents a 2D Control Volume cell (arbitrary simple polygon: triangle, quad, ...).
struct Cell {
    double volume;              // Face-bounded area of the control volume, in (mesh length units)^2
    double x_centroid;          // Center of mass coordinates, in mesh length units
    double y_centroid;
    std::vector<int> faces;     // Indices (into UnstructuredMesh::faces) of the faces bounding this cell
    std::vector<int> node_ids;  // Ordered (CCW) polygon vertex indices, used for geometry + output
};

// Boundary condition kind assigned to a patch, specified as an integer code
// in the case input file (see CaseInput.h).
enum class BoundaryType : int {
    Dirichlet = 0, // Fixed field value at the boundary (value = prescribed phi)
    Neumann = 1    // Fixed diffusive flux through the boundary (value = -alpha * dPhi/dn)
};

// A named group of boundary faces (matches a Gmsh Physical Group of dimension 1)
// sharing a common boundary condition, assigned from the case input file.
struct BoundaryPatch {
    std::string name;                             // Patch name, matched against the case file's "boundary" entries
    BoundaryType type = BoundaryType::Dirichlet;   // Kind of boundary condition applied on this patch
    double value = 0.0;                            // Prescribed value or flux (units depend on 'type', see BoundaryType)
};


// Unified database containing the complete grid graph topology: the mesh
// nodes, the faces connecting cells (and their per-face geometry), the cells
// themselves, and the named boundary patches used to apply boundary conditions.
struct UnstructuredMesh {
    std::vector<Node> nodes;
    std::vector<Face> faces;
    std::vector<Cell> cells;
    std::vector<BoundaryPatch> patches;

    // Computes per-face geometric quantities (midpoint, length/area, outward unit
    // normal) from the node coordinates already stored on each face.
    //
    // Methodology: for a 2D face defined by two nodes n1->n2, the face "area"
    // (really a length, since a 2D face is a line segment) is |n2 - n1|, and the
    // outward normal is that edge vector rotated -90 degrees and normalized:
    // (nx, ny) = (dy, -dx) / |edge|. This convention is only correct if node1/node2
    // are ordered consistently with the CCW winding of cell_left's polygon
    // (guaranteed by MeshReader, which builds faces by walking each cell's
    // node_ids in order).
    //
    // Cell volumes/centroids are NOT computed here: they are filled in directly
    // during mesh parsing (see MeshReader::read_gmsh) via the shoelace formula,
    // since that is where the ordered per-cell node list is available.
    //
    // Input:  faces[i].node1/node2 and nodes[] must already be populated.
    // Output: faces[i].x_mid/y_mid/area/nx/ny are (re)computed in place.
    void compute_geometry() {
        for (auto& face : faces) {
            const Node& n1 = nodes[face.node1];
            const Node& n2 = nodes[face.node2];

            face.x_mid = 0.5 * (n1.x + n2.x);
            face.y_mid = 0.5 * (n1.y + n2.y);

            double dx = n2.x - n1.x;
            double dy = n2.y - n1.y;
            face.area = std::sqrt(dx*dx + dy*dy);

            // Outward normal vector components (edge vector rotated -90 degrees, then normalized)
            face.nx = dy / face.area;
            face.ny = -dx / face.area;
        }
    }
};


#endif // UNSTRUCTUREDMESH_H_INCLUDED
