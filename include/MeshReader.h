// SPDX-License-Identifier: GPL-3.0-only
#ifndef MESHREADER_H_INCLUDED
#define MESHREADER_H_INCLUDED

#include <string>
#include "UnstructuredMesh.h"
#include "MeshCheckReport.h"

// Parses Gmsh ASCII (.msh) files into an UnstructuredMesh. Both legacy format
// 2.2 and format 4.x are supported; the version is auto-detected from the
// file's $MeshFormat header (defaulting to 2.2 if that section is absent).
// Binary mesh files are not supported.
//
// $Nodes            -> mesh.nodes
// $Entities         -> (format 4.x only) per-entity physical group tags
// $Elements (2D)    -> mesh.cells (triangles/quads), mesh.faces derived from shared edges
// $Elements (1D)    -> assigns Face::patch_id to the boundary edges they cover
// $PhysicalNames    -> mesh.patches (only dimension-1 groups become boundary patches)
namespace MeshReader {
    // Reads a Gmsh ASCII mesh file (format 2.2 or 4.x) and fills 'mesh' with
    // its topology and geometry (cell volumes/centroids, face midpoints/areas/
    // normals are NOT computed here for faces -- call
    // UnstructuredMesh::compute_geometry() afterwards; cell volumes/centroids
    // ARE computed here since they need the ordered per-cell node list that is
    // only available at parse time).
    //
    // Input:
    //   filename - path to a .msh file (Gmsh ASCII format 2.2 or 4.x)
    // Output:
    //   mesh     - populated with nodes, faces, cells and boundary patches
    // Returns:
    //   true on success; false if the file could not be opened, uses an
    //   unsupported format version/feature (e.g. binary, parametric node
    //   coordinates), contains no nodes/cells, has an untagged boundary edge,
    //   or fails the post-build geometry validation (non-finite node
    //   coordinate, or a degenerate zero/non-finite cell volume or face
    //   length -- see validate_mesh_geometry() in MeshReader.cpp).
    //
    // 'report', if non-null, records --validate-mesh's 3-step check
    // progress/issues as parsing/connectivity-building/geometry-validation
    // complete (see MeshCheckReport.h and docs/mesh-validation-criteria.md).
    // Real solver runs never pass one; behavior is identical to omitting it.
    bool read_gmsh(const std::string& filename, UnstructuredMesh& mesh, MeshCheckReport* report = nullptr);

    // Reads this project's own FVMESH format (see docs/fvmesh-format.md) into
    // an UnstructuredMesh. Unlike read_gmsh, cells may have any number of
    // vertices/faces (no triangle/quad restriction), and boundary edges name
    // their patch directly (no numeric physical-group indirection). The
    // "# FVMESH <major>.<minor>" header is version-checked (currently only
    // 1.0 is supported), the same way read_gmsh dispatches on $MeshFormat.
    //
    // Input:
    //   filename - path to a .fvmesh file
    // Output:
    //   mesh     - populated with nodes, faces, cells and boundary patches
    //              (same post-conditions as read_gmsh: face midpoint/area/
    //              normal still require a separate compute_geometry() call)
    // Returns:
    //   true on success; false if the file could not be opened, its header
    //   is missing/malformed/an unsupported version, a CELLS/BOUNDARY entry
    //   references an out-of-range node index, has an untagged boundary edge,
    //   or fails the post-build geometry validation (same checks as
    //   read_gmsh -- see validate_mesh_geometry() in MeshReader.cpp). A
    //   BOUNDARY edge that matches no cell edge is silently ignored (same
    //   behavior as read_gmsh's boundary line elements).
    //
    // 'report' -- see read_gmsh's parameter comment above.
    bool read_fvmesh(const std::string& filename, UnstructuredMesh& mesh, MeshCheckReport* report = nullptr);

    // Dispatches to read_gmsh or read_fvmesh based on 'filename's extension
    // (".msh" -> read_gmsh, ".fvmesh" -> read_fvmesh).
    //
    // 'report' -- see read_gmsh's parameter comment above. If parsing fails
    // before the shared build_cells_faces_patches() stage is reached (e.g.
    // file not found, malformed header), this falls back to marking
    // report's Parsing step failed itself, since in that case neither
    // read_gmsh nor read_fvmesh got far enough to do so.
    //
    // Returns: false immediately (no file access attempted) if the extension
    // is neither; otherwise the chosen reader's return value.
    bool read(const std::string& filename, UnstructuredMesh& mesh, MeshCheckReport* report = nullptr);
}

#endif // MESHREADER_H_INCLUDED
