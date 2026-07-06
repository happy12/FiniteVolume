// SPDX-License-Identifier: GPL-3.0-only
#ifndef VTKWRITER_H_INCLUDED
#define VTKWRITER_H_INCLUDED

#include <string>
#include <vector>
#include "UnstructuredMesh.h"

// A named per-cell scalar field, one value per mesh.cells entry, for the
// multi-field write() overload below.
struct NamedField {
    std::string name;
    std::vector<double> values;
};

// Writes results as a legacy ASCII VTK UNSTRUCTURED_GRID file, openable
// directly in ParaView: cell polygons with one or more cell-data scalar fields.
namespace VtkWriter {
    // Input:
    //   filename  - path to write the .vtk file to
    //   mesh      - provides node coordinates and, per cell, the ordered
    //               node_ids polygon used to emit VTK_POLYGON cells
    //   phi       - solver state field, one value per mesh.cells entry, written
    //               as CELL_DATA under the name "phi" (units are whatever the
    //               transported quantity is)
    //   precision - significant digits for every written double (node
    //               coordinates and phi values); valid range 1-17, default 6
    // Output: none (writes to disk)
    // Returns: true if the file was opened and written; false otherwise.
    bool write(const std::string& filename, const UnstructuredMesh& mesh, const std::vector<double>& phi,
               int precision = 6);

    // Input:
    //   filename  - path to write the .vtk file to
    //   mesh      - provides node coordinates and, per cell, the ordered
    //               node_ids polygon used to emit VTK_POLYGON cells
    //   fields    - one or more named per-cell scalar fields, each written as
    //               its own CELL_DATA SCALARS block
    //   precision - significant digits for every written double (node
    //               coordinates and field values); valid range 1-17, default 6
    // Output: none (writes to disk)
    // Returns: true if the file was opened and written; false otherwise.
    bool write(const std::string& filename, const UnstructuredMesh& mesh, const std::vector<NamedField>& fields,
               int precision = 6);
}

#endif // VTKWRITER_H_INCLUDED
