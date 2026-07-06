// SPDX-License-Identifier: GPL-3.0-only
#include "FvMeshWriter.h"

#include <fstream>
#include <iomanip>

// See FvMeshWriter.h for the input/output contract, and docs/fvmesh-format.md
// for the grammar being emitted.
bool FvMeshWriter::write(const std::string& filename, const UnstructuredMesh& mesh, int precision) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        return false;
    }
    out << std::setprecision(precision);

    out << "# FVMESH 1.0\n";

    out << "POINTS " << mesh.nodes.size() << "\n";
    for (const auto& node : mesh.nodes) {
        out << node.x << " " << node.y << "\n";
    }

    out << "CELLS " << mesh.cells.size() << "\n";
    for (const auto& cell : mesh.cells) {
        out << cell.node_ids.size();
        for (int id : cell.node_ids) {
            out << " " << id;
        }
        out << "\n";
    }

    // Boundary edges are exactly the faces with no right-hand cell and a
    // patch assigned; internal faces (cell_right != -1) are re-derived by
    // the reader from shared cell edges and must not be listed here.
    size_t boundary_count = 0;
    for (const auto& face : mesh.faces) {
        if (face.cell_right == -1 && face.patch_id != -1) ++boundary_count;
    }

    out << "BOUNDARY " << boundary_count << "\n";
    for (const auto& face : mesh.faces) {
        if (face.cell_right == -1 && face.patch_id != -1) {
            out << face.node1 << " " << face.node2 << " " << mesh.patches[face.patch_id].name << "\n";
        }
    }

    return true;
}
