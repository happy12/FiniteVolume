// SPDX-License-Identifier: GPL-3.0-only
#ifndef FVMESHWRITER_H_INCLUDED
#define FVMESHWRITER_H_INCLUDED

#include <string>
#include "UnstructuredMesh.h"

// Writes this project's own FVMESH 1.0 mesh format (see docs/fvmesh-format.md),
// readable back via MeshReader::read_fvmesh. Intended for round-trip testing
// against the Gmsh reader, and as a reference for an eventual in-house mesh
// generator targeting this format.
namespace FvMeshWriter {
    // Input:
    //   filename  - path to write the .fvmesh file to
    //   mesh      - provides node coordinates, per-cell ordered node_ids
    //               polygons, and boundary faces (cell_right == -1 with a
    //               patch_id set) with their owning patch name
    //   precision - significant digits for every written node coordinate;
    //               valid range 1-17, default 6
    // Output: none (writes to disk)
    // Returns: true if the file was opened and written; false otherwise.
    bool write(const std::string& filename, const UnstructuredMesh& mesh, int precision = 6);
}

#endif // FVMESHWRITER_H_INCLUDED
