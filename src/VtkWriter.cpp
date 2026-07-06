// SPDX-License-Identifier: GPL-3.0-only
#include "VtkWriter.h"

#include <fstream>
#include <iomanip>

#include "Version.h"

namespace {
    // VTK cell type code for an arbitrary simple polygon (legacy VTK file format).
    const int VTK_POLYGON = 7;

    // Emits the legacy VTK "DATASET UNSTRUCTURED_GRID" geometry shared by both
    // write() overloads below:
    //   POINTS      - one line per mesh node (x y z=0, since this is a 2D solver)
    //   CELLS       - one line per cell: "<vertex count> <vertex indices...>",
    //                 using each cell's ordered node_ids polygon directly
    //   CELL_TYPES  - VTK_POLYGON for every cell (works for triangles, quads, ...)
    //
    // Input:  mesh - provides node coordinates (mesh length units; z is
    //                always written as 0 for this 2D solver) and, per cell,
    //                the ordered node_ids polygon
    // Output: out  - the DATASET/POINTS/CELLS/CELL_TYPES blocks are appended
    //                (at the stream's current precision -- callers set that
    //                before calling this); caller is responsible for
    //                following up with a CELL_DATA block (neither overload's
    //                field data is written here)
    void write_grid_header(std::ofstream& out, const UnstructuredMesh& mesh) {
        out << "# vtk DataFile Version 3.0\n";
        out << "FiniteVolume " FV_VERSION_FULL "\n";
        out << "ASCII\n";
        out << "DATASET UNSTRUCTURED_GRID\n";

        // Node coordinates (z is always 0 for this 2D solver).
        out << "POINTS " << mesh.nodes.size() << " double\n";
        for (const auto& node : mesh.nodes) {
            out << node.x << " " << node.y << " 0.0\n";
        }

        // VTK's CELLS block is prefixed by the total entry count: one "vertex
        // count" value plus each cell's vertex indices, summed over all cells.
        size_t total_size = 0;
        for (const auto& cell : mesh.cells) {
            total_size += cell.node_ids.size() + 1;
        }

        out << "CELLS " << mesh.cells.size() << " " << total_size << "\n";
        for (const auto& cell : mesh.cells) {
            out << cell.node_ids.size();
            for (int id : cell.node_ids) {
                out << " " << id;
            }
            out << "\n";
        }

        // Every cell is a generic polygon, regardless of its vertex count.
        out << "CELL_TYPES " << mesh.cells.size() << "\n";
        for (size_t i = 0; i < mesh.cells.size(); ++i) {
            out << VTK_POLYGON << "\n";
        }
    }
}

// See VtkWriter.h for the input/output contract.
bool VtkWriter::write(const std::string& filename, const UnstructuredMesh& mesh, const std::vector<double>& phi,
                       int precision) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        return false;
    }
    out << std::setprecision(precision);

    write_grid_header(out, mesh);

    // Solver state field phi, one scalar per cell (same order as mesh.cells).
    out << "CELL_DATA " << mesh.cells.size() << "\n";
    out << "SCALARS phi double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (double value : phi) {
        out << value << "\n";
    }

    return true;
}

// See VtkWriter.h for the input/output contract.
bool VtkWriter::write(const std::string& filename, const UnstructuredMesh& mesh, const std::vector<NamedField>& fields,
                       int precision) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        return false;
    }
    out << std::setprecision(precision);

    write_grid_header(out, mesh);

    // Legacy VTK allows multiple SCALARS blocks to share one CELL_DATA count.
    out << "CELL_DATA " << mesh.cells.size() << "\n";
    for (const auto& field : fields) {
        out << "SCALARS " << field.name << " double 1\n";
        out << "LOOKUP_TABLE default\n";
        for (double value : field.values) {
            out << value << "\n";
        }
    }

    return true;
}
