// SPDX-License-Identifier: GPL-3.0-only
#include "MeshReader.h"

#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace {

// A raw Gmsh $Elements entry before it is sorted into mesh.cells (2D
// elements) or used to tag boundary faces (1D line elements).
struct TempElement {
    int type;                  // Gmsh element type code (1=line, 2=triangle, 3=quad)
    int physical_id;           // Gmsh physical group id (first element tag), or -1 if untagged
    std::vector<int> node_ids; // 0-based mesh node indices (already remapped from Gmsh's 1-based ids)
};

// Computes the signed area of a simple polygon via the shoelace formula.
//
// Methodology: Area = 1/2 * sum_i (x_i * y_{i+1} - x_{i+1} * y_i). The sign
// is positive for counter-clockwise (CCW) vertex winding and negative for
// clockwise -- callers that need an unsigned area should take fabs().
//
// Input:  mesh - provides node coordinates; ids - ordered polygon vertex indices
// Output: none (mesh unmodified)
// Returns: signed polygon area, in (mesh length units)^2
double polygon_signed_area(const UnstructuredMesh& mesh, const std::vector<int>& ids) {
    double a = 0.0;
    int n = (int)ids.size();
    for (int i = 0; i < n; ++i) {
        const Node& p0 = mesh.nodes[ids[i]];
        const Node& p1 = mesh.nodes[ids[(i + 1) % n]];
        a += p0.x * p1.y - p1.x * p0.y;
    }
    return 0.5 * a;
}

// Populates volume/centroid for a cell via the polygon shoelace formula
// (works for arbitrary convex or concave simple polygons, any winding).
//
// Methodology: the standard shoelace centroid formula
//   Cx = 1/(6A) * sum_i (x_i + x_{i+1}) * cross_i
//   Cy = 1/(6A) * sum_i (y_i + y_{i+1}) * cross_i
// with cross_i = x_i*y_{i+1} - x_{i+1}*y_i and A the SIGNED area (not fabs),
// so the result is correct regardless of the polygon's winding direction.
//
// Input:  mesh.cells[cell_index].node_ids must already be populated (ordered
//         CCW/CW polygon vertex list) and mesh.nodes must hold their coordinates.
// Output: mesh.cells[cell_index].volume/x_centroid/y_centroid are filled in.
//         volume is in (mesh length units)^2, centroid in mesh length units.
void compute_cell_geometry(UnstructuredMesh& mesh, int cell_index) {
    Cell& cell = mesh.cells[cell_index];
    const std::vector<int>& ids = cell.node_ids;
    int n = (int)ids.size();

    double area = polygon_signed_area(mesh, ids);
    double cx = 0.0, cy = 0.0;
    for (int i = 0; i < n; ++i) {
        const Node& p0 = mesh.nodes[ids[i]];
        const Node& p1 = mesh.nodes[ids[(i + 1) % n]];
        double cross = p0.x * p1.y - p1.x * p0.y;
        cx += (p0.x + p1.x) * cross;
        cy += (p0.y + p1.y) * cross;
    }
    cx /= (6.0 * area);
    cy /= (6.0 * area);

    cell.volume = std::fabs(area);
    cell.x_centroid = cx;
    cell.y_centroid = cy;
}

// Validates basic mesh geometry once nodes/cells/faces are built, catching
// the failure mode CLAUDE.md documents as a known gap: a degenerate mesh
// (collinear/coincident cell vertices, a zero-length face) otherwise loads
// without complaint and only surfaces later as a silent NaN/Inf once
// UnstructuredMesh::compute_geometry() divides a zero face length into its
// normal, or a solver divides by a zero cell volume. Only exactly-degenerate
// (zero or non-finite) geometry is a hard error; a merely very small nonzero
// face/cell is reported as a warning, since a legitimately stretched
// boundary-layer mesh can have very thin cells without being invalid (see
// --verify-ns-stretched-cfl).
//
// Input:  mesh - nodes/cells/faces must already be populated (cell volumes
//         via compute_cell_geometry above; face.area/nx/ny are NOT yet
//         populated at this point in the pipeline -- see
//         UnstructuredMesh::compute_geometry -- so face length is computed
//         fresh here from node coordinates instead of relying on that field)
// Output: none (mesh unmodified); diagnostics printed to stderr
// Returns: true if no hard error was found (warnings do not fail the load); false otherwise
bool validate_mesh_geometry(const UnstructuredMesh& mesh) {
    bool ok = true;

    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        const Node& n = mesh.nodes[i];
        if (!std::isfinite(n.x) || !std::isfinite(n.y)) {
            std::cerr << "Error: mesh node " << i << " has a non-finite coordinate ("
                       << n.x << ", " << n.y << ")\n";
            ok = false;
        }
    }
    if (!ok) return false; // cell/face geometry below is meaningless with garbage coordinates

    // Exact-duplicate node positions (distinct indices, identical
    // coordinates) aren't fatal by themselves, but usually indicate a
    // mesh-generation mistake worth flagging.
    std::vector<int> order(mesh.nodes.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = (int)i;
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (mesh.nodes[a].x != mesh.nodes[b].x) return mesh.nodes[a].x < mesh.nodes[b].x;
        return mesh.nodes[a].y < mesh.nodes[b].y;
    });
    for (size_t i = 1; i < order.size(); ++i) {
        const Node& a = mesh.nodes[order[i - 1]];
        const Node& b = mesh.nodes[order[i]];
        if (a.x == b.x && a.y == b.y) {
            std::cerr << "Warning: mesh nodes " << order[i - 1] << " and " << order[i]
                       << " are at the same coordinate (" << a.x << ", " << a.y << ")\n";
        }
    }

    // Cell volumes: exactly zero/non-finite (a fully collinear/collapsed
    // polygon) is a hard error; merely very small relative to the mesh's own
    // mean cell volume is only a warning.
    double volume_sum = 0.0;
    int finite_volume_count = 0;
    for (const Cell& cell : mesh.cells) {
        if (std::isfinite(cell.volume) && cell.volume > 0.0) {
            volume_sum += cell.volume;
            ++finite_volume_count;
        }
    }
    double mean_volume = (finite_volume_count > 0) ? volume_sum / finite_volume_count : 0.0;

    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        double volume = mesh.cells[c].volume;
        if (!std::isfinite(volume) || volume <= 0.0) {
            std::cerr << "Error: cell " << c << " has a degenerate (zero or non-finite) volume ("
                       << volume << ") -- its polygon vertices are collinear or coincident\n";
            ok = false;
        } else if (mean_volume > 0.0 && volume < 1e-8 * mean_volume) {
            std::cerr << "Warning: cell " << c << " has volume " << volume
                       << ", far smaller than the mesh's mean cell volume (" << mean_volume
                       << ") -- check for a near-degenerate element\n";
        }
    }

    // Face lengths: exactly zero/non-finite (coincident endpoints) is a hard
    // error -- it would otherwise divide by zero in compute_geometry()'s
    // normal calculation; merely very small is only a warning.
    std::vector<double> lengths(mesh.faces.size());
    double length_sum = 0.0;
    int finite_length_count = 0;
    for (size_t f = 0; f < mesh.faces.size(); ++f) {
        const Face& face = mesh.faces[f];
        double dx = mesh.nodes[face.node2].x - mesh.nodes[face.node1].x;
        double dy = mesh.nodes[face.node2].y - mesh.nodes[face.node1].y;
        double length = std::sqrt(dx * dx + dy * dy);
        lengths[f] = length;
        if (std::isfinite(length) && length > 0.0) {
            length_sum += length;
            ++finite_length_count;
        }
    }
    double mean_length = (finite_length_count > 0) ? length_sum / finite_length_count : 0.0;

    for (size_t f = 0; f < mesh.faces.size(); ++f) {
        double length = lengths[f];
        if (!std::isfinite(length) || length <= 0.0) {
            std::cerr << "Error: face " << f << " (nodes " << mesh.faces[f].node1 << ", "
                       << mesh.faces[f].node2 << ") has a degenerate (zero or non-finite) length ("
                       << length << ")\n";
            ok = false;
        } else if (mean_length > 0.0 && length < 1e-8 * mean_length) {
            std::cerr << "Warning: face " << f << " (nodes " << mesh.faces[f].node1 << ", "
                       << mesh.faces[f].node2 << ") has length " << length
                       << ", far smaller than the mesh's mean face length (" << mean_length
                       << ") -- check for a near-degenerate element\n";
        }
    }

    return ok;
}

// Builds mesh.cells/mesh.faces/mesh.patches from the raw elements collected
// by either format's section parser. Shared so that the 2.2 and 4.1 readers
// (which differ only in how they populate 'cell_elements'/'line_elements'/
// 'physical_names') are guaranteed to produce identical mesh topology.
//
// Methodology: see MeshReader::read_gmsh's class comment -- cells are built
// directly from 'cell_elements' (preserving node winding), faces are derived
// by walking each cell's polygon edge by edge and matching shared edges via a
// node-pair key, and boundary line elements are matched back to the derived
// boundary faces by the same key to assign patch_id (creating BoundaryPatch
// entries from 'physical_names' lazily, as needed).
//
// Input:  mesh            - nodes must already be populated
//         cell_elements   - 2D (triangle/quad) elements to become mesh.cells
//         line_elements   - 1D (line) elements tagging boundary edges
//         physical_names  - physical group id -> name, for patches
// Output: mesh.cells/faces/patches are filled in
// Returns: false if any boundary face is untagged or validate_mesh_geometry()
//          finds a degenerate (zero/non-finite volume or face length) element
//          (a descriptive message is printed to stderr in either case);
//          true otherwise
bool build_cells_faces_patches(UnstructuredMesh& mesh,
                                const std::vector<TempElement>& cell_elements,
                                const std::vector<TempElement>& line_elements,
                                const std::map<int, std::string>& physical_names) {
    // --- Build cells and derive faces from shared edges ---
    // edge_to_face maps an (unordered) node-index pair to the face it has
    // already produced, so that the second cell to visit a shared edge links
    // up with the first instead of creating a duplicate face.
    mesh.cells.resize(cell_elements.size());
    std::map<std::pair<int, int>, int> edge_to_face;

    for (size_t ci = 0; ci < cell_elements.size(); ++ci) {
        const std::vector<int>& ids = cell_elements[ci].node_ids;
        mesh.cells[ci].node_ids = ids; // Preserve CCW/CW ordering for geometry + VTK output

        // Walk the cell's polygon edge by edge (node[e] -> node[e+1], wrapping around).
        int n = (int)ids.size();
        for (int e = 0; e < n; ++e) {
            int a = ids[e];
            int b = ids[(e + 1) % n];
            auto key = std::make_pair(std::min(a, b), std::max(a, b));

            auto it = edge_to_face.find(key);
            if (it == edge_to_face.end()) {
                // First cell to visit this edge: create the face, owned as cell_left.
                // node1/node2 keep THIS cell's traversal order, which is what
                // UnstructuredMesh::compute_geometry() needs to get the outward
                // normal direction (relative to cell_left) right.
                Face face;
                face.node1 = a;
                face.node2 = b;
                face.cell_left = (int)ci;
                face.cell_right = -1;
                int face_index = (int)mesh.faces.size();
                mesh.faces.push_back(face);
                mesh.cells[ci].faces.push_back(face_index);
                edge_to_face[key] = face_index;
            } else {
                // Second cell to visit this edge: it's the neighbor across an
                // internal face. Adjacent CCW polygons always traverse a shared
                // edge in opposite directions, so node1/node2 are left untouched.
                int face_index = it->second;
                mesh.faces[face_index].cell_right = (int)ci;
                mesh.cells[ci].faces.push_back(face_index);
            }
        }

        compute_cell_geometry(mesh, (int)ci);
    }

    // --- Create boundary patches from physical names, assign patch_id to boundary faces ---
    // Every collected line element represents an explicitly-tagged boundary
    // edge; match it back to the Face built above by the same node-pair key,
    // and stamp that face with the (lazily created) patch index for its
    // physical group.
    std::map<int, int> physical_id_to_patch; // physical id -> mesh.patches index
    for (const auto& elem : line_elements) {
        auto key = std::make_pair(std::min(elem.node_ids[0], elem.node_ids[1]),
                                   std::max(elem.node_ids[0], elem.node_ids[1]));
        auto face_it = edge_to_face.find(key);
        if (face_it == edge_to_face.end()) continue; // edge not part of any cell boundary

        int patch_index;
        auto patch_it = physical_id_to_patch.find(elem.physical_id);
        if (patch_it == physical_id_to_patch.end()) {
            BoundaryPatch patch;
            auto name_it = physical_names.find(elem.physical_id);
            patch.name = (name_it != physical_names.end()) ? name_it->second
                                                             : ("patch_" + std::to_string(elem.physical_id));
            patch_index = (int)mesh.patches.size();
            mesh.patches.push_back(patch);
            physical_id_to_patch[elem.physical_id] = patch_index;
        } else {
            patch_index = patch_it->second;
        }

        mesh.faces[face_it->second].patch_id = patch_index;
    }

    // --- Reject domain-boundary faces that never got a patch tag ---
    // A face with cell_right == -1 is a genuine domain boundary; if its
    // patch_id is still -1 here, its edge never appeared in the tagged
    // boundary list above, and downstream solvers index bcs[]/patches[] by
    // patch_id without a bounds check. Fail loudly now rather than let that
    // reach an out-of-bounds access deep into a run.
    int untagged_count = 0;
    for (const Face& face : mesh.faces) {
        if (face.cell_right == -1 && face.patch_id == -1) ++untagged_count;
    }
    if (untagged_count > 0) {
        std::cerr << "Error: mesh has " << untagged_count << " boundary faces with no assigned "
                     "patch (untagged boundary edges). Every boundary face must belong to a "
                     "tagged patch.\n";
        return false;
    }

    if (!validate_mesh_geometry(mesh)) {
        return false;
    }

    return true;
}

// Reads a Gmsh legacy ASCII format 2.2 mesh from an already-open stream
// (positioned at the start of the file). See MeshReader::read_gmsh's class
// comment for the section-by-section methodology.
//
// Input:  in   - open input stream, positioned at offset 0
// Output: mesh - populated with nodes, faces, cells and boundary patches
// Returns: true on success; false if the file contains no nodes/cells.
bool read_gmsh_v22(std::ifstream& in, UnstructuredMesh& mesh) {
    std::unordered_map<int, int> node_id_to_index; // gmsh node id -> mesh.nodes index
    std::map<int, std::string> physical_names;      // physical id -> name (dimension 1 only)
    std::vector<TempElement> cell_elements;          // triangles/quads
    std::vector<TempElement> line_elements;          // boundary edges

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        // $PhysicalNames block: maps a physical group id to its name for
        // dimension-1 (boundary edge) groups; dimension-2 groups (the fluid
        // domain itself) are read but discarded, since cells don't need a name.
        if (line.rfind("$PhysicalNames", 0) == 0) {
            int count = 0;
            in >> count;
            std::getline(in, line); // consume rest of count line
            for (int i = 0; i < count; ++i) {
                std::getline(in, line);
                std::istringstream iss(line);
                int dim, id;
                std::string name;
                iss >> dim >> id;
                std::getline(iss, name);
                size_t first_quote = name.find('"');
                size_t last_quote = name.rfind('"');
                if (first_quote != std::string::npos && last_quote != std::string::npos && last_quote > first_quote) {
                    name = name.substr(first_quote + 1, last_quote - first_quote - 1);
                }
                if (dim == 1) {
                    physical_names[id] = name;
                }
            }
        // $Nodes block: "id x y z" per line (z is discarded, this is a 2D solver).
        // Node ids in Gmsh are not guaranteed contiguous/sorted, so we keep an
        // explicit id -> 0-based index map rather than assuming id-1 == index.
        } else if (line.rfind("$Nodes", 0) == 0) {
            int count = 0;
            in >> count;
            mesh.nodes.resize(count);
            for (int i = 0; i < count; ++i) {
                int id; double x, y, z;
                in >> id >> x >> y >> z;
                node_id_to_index[id] = i;
                mesh.nodes[i] = {x, y};
            }
            std::getline(in, line); // consume rest of last node line
        // $Elements block: "elm_number elm_type num_tags tag... node...".
        // Only the first tag (the physical group id) is kept; only element
        // types 1 (2-node line), 2 (3-node triangle) and 3 (4-node quad) are
        // understood -- anything else (points, higher-order elements, 3D
        // elements, ...) is skipped since this is a 2D, linear-element solver.
        } else if (line.rfind("$Elements", 0) == 0) {
            int count = 0;
            in >> count;
            std::getline(in, line); // consume rest of count line
            for (int i = 0; i < count; ++i) {
                std::getline(in, line);
                std::istringstream iss(line);
                int elm_number, elm_type, num_tags;
                iss >> elm_number >> elm_type >> num_tags;

                int physical_id = -1;
                for (int t = 0; t < num_tags; ++t) {
                    int tag;
                    iss >> tag;
                    if (t == 0) physical_id = tag;
                }

                int node_count = 0;
                if (elm_type == 1) node_count = 2;      // 2-node line
                else if (elm_type == 2) node_count = 3; // 3-node triangle
                else if (elm_type == 3) node_count = 4; // 4-node quadrangle
                else continue;                          // unsupported element type, skip

                TempElement elem;
                elem.type = elm_type;
                elem.physical_id = physical_id;
                elem.node_ids.resize(node_count);
                for (int k = 0; k < node_count; ++k) {
                    int nid;
                    iss >> nid;
                    elem.node_ids[k] = node_id_to_index.at(nid);
                }

                if (elm_type == 1) line_elements.push_back(std::move(elem));
                else cell_elements.push_back(std::move(elem));
            }
        }
    }

    if (mesh.nodes.empty() || cell_elements.empty()) {
        std::cerr << "Error reading mesh: no nodes or cells were found in the $Nodes/$Elements sections\n";
        return false;
    }

    return build_cells_faces_patches(mesh, cell_elements, line_elements, physical_names);
}

// Reads a Gmsh ASCII format 4.x mesh from an already-open stream (positioned
// at the start of the file).
//
// Format differences from 2.2 (see MeshReader::read_gmsh's class comment):
//   $Entities   -- new section; each curve (dim 1) / surface (dim 2) entity
//                  lists the physical group tags that apply to every element
//                  built on that entity (elements carry no per-element tags
//                  of their own, unlike 2.2). Point (dim 0) and volume (dim 3)
//                  entities are consumed to keep the stream aligned but not
//                  recorded, since this is a 2D, linear-element solver.
//   $Nodes      -- grouped into per-entity blocks; each block lists its node
//                  tags, then separately lists their "x y z" coordinates.
//   $Elements   -- grouped into per-entity blocks; the element type is a
//                  per-block header field (not per-element), and each element
//                  line is just "elementTag nodeTag...".
//
// Input:  in   - open input stream, positioned at offset 0
// Output: mesh - populated with nodes, faces, cells and boundary patches
// Returns: true on success; false if the file contains no nodes/cells, or
//          uses a feature this reader does not support (parametric node
//          coordinates; a stderr message is printed in that case).
bool read_gmsh_v41(std::ifstream& in, UnstructuredMesh& mesh) {
    std::unordered_map<long long, int> node_id_to_index; // gmsh node tag -> mesh.nodes index
    std::map<int, std::string> physical_names;            // physical id -> name (dimension 1 only)
    std::vector<TempElement> cell_elements;               // triangles/quads
    std::vector<TempElement> line_elements;               // boundary edges

    // (entityDim, entityTag) -> physical group ids assigned to that entity in $Entities.
    std::map<std::pair<int, int>, std::vector<int>> entity_physical_tags;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        // $PhysicalNames block: identical format to 2.2.
        if (line.rfind("$PhysicalNames", 0) == 0) {
            int count = 0;
            in >> count;
            std::getline(in, line); // consume rest of count line
            for (int i = 0; i < count; ++i) {
                std::getline(in, line);
                std::istringstream iss(line);
                int dim, id;
                std::string name;
                iss >> dim >> id;
                std::getline(iss, name);
                size_t first_quote = name.find('"');
                size_t last_quote = name.rfind('"');
                if (first_quote != std::string::npos && last_quote != std::string::npos && last_quote > first_quote) {
                    name = name.substr(first_quote + 1, last_quote - first_quote - 1);
                }
                if (dim == 1) {
                    physical_names[id] = name;
                }
            }
        // $Entities block: one line per entity, in order dim 0 (points), 1
        // (curves), 2 (surfaces), 3 (volumes). Only curves/surfaces are
        // recorded, keyed by (dim, tag) -> their physical group tag list;
        // points/volumes are skipped (not relevant to a 2D solver) but their
        // lines must still be consumed to stay aligned with the file.
        } else if (line.rfind("$Entities", 0) == 0) {
            int num_points = 0, num_curves = 0, num_surfaces = 0, num_volumes = 0;
            in >> num_points >> num_curves >> num_surfaces >> num_volumes;
            std::getline(in, line); // consume rest of the counts line

            for (int i = 0; i < num_points; ++i) {
                std::getline(in, line); // "tag x y z numPhysicalTags physicalTag..." -- discarded
            }
            for (int dim : {1, 2}) {
                int count = (dim == 1) ? num_curves : num_surfaces;
                for (int i = 0; i < count; ++i) {
                    std::getline(in, line);
                    std::istringstream iss(line);
                    int tag;
                    double min_x, min_y, min_z, max_x, max_y, max_z;
                    iss >> tag >> min_x >> min_y >> min_z >> max_x >> max_y >> max_z;
                    int num_physical_tags = 0;
                    iss >> num_physical_tags;
                    std::vector<int> tags(num_physical_tags);
                    for (int t = 0; t < num_physical_tags; ++t) {
                        iss >> tags[t];
                    }
                    // (the trailing "numBoundingEntities entityTag..." list is left unparsed)
                    entity_physical_tags[{dim, tag}] = std::move(tags);
                }
            }
            for (int i = 0; i < num_volumes; ++i) {
                std::getline(in, line); // discarded -- 3D, out of scope
            }
        // $Nodes block: per-entity blocks of node tags followed by their
        // coordinates. Parametric coordinates (extra u/v/w columns) are not
        // supported; a block requesting them fails the parse rather than
        // silently desyncing the reader.
        } else if (line.rfind("$Nodes", 0) == 0) {
            long long num_entity_blocks = 0, num_nodes = 0, min_tag = 0, max_tag = 0;
            in >> num_entity_blocks >> num_nodes >> min_tag >> max_tag;
            std::getline(in, line); // consume rest of the header line
            mesh.nodes.reserve((size_t)num_nodes);

            for (long long b = 0; b < num_entity_blocks; ++b) {
                int entity_dim, entity_tag, parametric;
                long long num_in_block;
                in >> entity_dim >> entity_tag >> parametric >> num_in_block;
                std::getline(in, line); // consume rest of the block header line

                if (parametric != 0) {
                    std::cerr << "Error reading mesh: parametric $Nodes coordinates are not supported\n";
                    return false;
                }

                std::vector<long long> tags((size_t)num_in_block);
                for (long long i = 0; i < num_in_block; ++i) {
                    in >> tags[(size_t)i];
                }
                std::getline(in, line); // consume rest of the last tag line

                for (long long i = 0; i < num_in_block; ++i) {
                    double x, y, z;
                    in >> x >> y >> z;
                    int index = (int)mesh.nodes.size();
                    node_id_to_index[tags[(size_t)i]] = index;
                    mesh.nodes.push_back({x, y});
                }
                std::getline(in, line); // consume rest of the last coordinate line
            }
        // $Elements block: per-entity blocks; the element type is fixed for
        // the whole block, and every element in it inherits the block's
        // entity's physical group (first tag, or -1 if none) -- see the
        // function comment above for why this reproduces 2.2's per-element
        // physical_id field.
        } else if (line.rfind("$Elements", 0) == 0) {
            long long num_entity_blocks = 0, num_elements = 0, min_tag = 0, max_tag = 0;
            in >> num_entity_blocks >> num_elements >> min_tag >> max_tag;
            std::getline(in, line); // consume rest of the header line

            for (long long b = 0; b < num_entity_blocks; ++b) {
                int entity_dim, entity_tag, elm_type;
                long long num_in_block;
                in >> entity_dim >> entity_tag >> elm_type >> num_in_block;
                std::getline(in, line); // consume rest of the block header line

                int node_count = 0;
                if (elm_type == 1) node_count = 2;      // 2-node line
                else if (elm_type == 2) node_count = 3; // 3-node triangle
                else if (elm_type == 3) node_count = 4; // 4-node quadrangle
                // else: unsupported element type -- node_count stays 0, lines
                // below are still consumed but no TempElement is built.

                int physical_id = -1;
                auto ent_it = entity_physical_tags.find({entity_dim, entity_tag});
                if (ent_it != entity_physical_tags.end() && !ent_it->second.empty()) {
                    physical_id = ent_it->second.front();
                }

                for (long long i = 0; i < num_in_block; ++i) {
                    std::getline(in, line);
                    if (node_count == 0) continue; // unsupported element type, skip

                    std::istringstream iss(line);
                    long long element_tag;
                    iss >> element_tag;

                    TempElement elem;
                    elem.type = elm_type;
                    elem.physical_id = physical_id;
                    elem.node_ids.resize(node_count);
                    for (int k = 0; k < node_count; ++k) {
                        long long nid;
                        iss >> nid;
                        elem.node_ids[k] = node_id_to_index.at(nid);
                    }

                    if (elm_type == 1) line_elements.push_back(std::move(elem));
                    else cell_elements.push_back(std::move(elem));
                }
            }
        }
    }

    if (mesh.nodes.empty() || cell_elements.empty()) {
        std::cerr << "Error reading mesh: no nodes or cells were found in the $Nodes/$Elements sections\n";
        return false;
    }

    return build_cells_faces_patches(mesh, cell_elements, line_elements, physical_names);
}

} // namespace

// See MeshReader.h for the input/output contract. Sniffs the mesh format
// version from the file's $MeshFormat header (defaulting to 2.2 if that
// section is absent, matching every mesh this reader supported before format
// 4.x was added), rewinds, then dispatches to the matching parser:
//
//   $PhysicalNames  dimension id "name"           -- boundary group names (dim 1 only)
//   $Entities       (format 4.x only)              -- per-entity physical group tags
//   $Nodes          id x y z (2.2) or per-block (4.x) -- vertex coordinates
//   $Elements       (2.2) or per-block (4.x)       -- cells (type 2/3) and boundary edges (type 1)
//
// Cells and faces are then derived in a second pass, identically for both
// formats (see build_cells_faces_patches above): for every 2D element we walk
// its node list edge by edge; the first cell to visit an edge (node-pair)
// creates a new Face, the second cell to visit the same edge (order-independent,
// via a sorted node-pair key) fills in that Face's cell_right, turning it into
// an internal face. Edges visited by only one cell remain boundary faces
// (cell_right == -1). Finally, the explicit boundary line elements are matched
// back to these boundary faces by the same node-pair key, to assign patch_id.
bool MeshReader::read_gmsh(const std::string& filename, UnstructuredMesh& mesh) {
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Error: could not open mesh file '" << filename << "'\n";
        return false;
    }

    double version = 2.2; // files with no $MeshFormat section are legacy 2.2
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("$MeshFormat", 0) == 0) {
            in >> version;
            std::getline(in, line); // consume rest of the version line
            break;
        }
        if (line.rfind("$PhysicalNames", 0) == 0 || line.rfind("$Nodes", 0) == 0) {
            break; // reached mesh data with no $MeshFormat section present
        }
    }
    in.clear();
    in.seekg(0); // both parsers re-scan the file from the top by section tag

    if (version == 2.2) {
        return read_gmsh_v22(in, mesh);
    }
    if (version >= 4.0) {
        return read_gmsh_v41(in, mesh);
    }

    std::cerr << "Error in mesh file '" << filename << "': unsupported Gmsh format version "
               << version << " (supported: 2.2, 4.x ASCII)\n";
    return false;
}

// See MeshReader.h for the input/output contract, and docs/fvmesh-format.md
// for the full grammar. Unlike read_gmsh, cells may have any number of
// vertices (no triangle/quad restriction) and boundary edges name their
// patch directly, so patch ids are synthesized here (interned by first
// appearance of each patch_name) before handing off to the same
// build_cells_faces_patches() the Gmsh readers use -- that function is
// otherwise completely format-agnostic.
bool MeshReader::read_fvmesh(const std::string& filename, UnstructuredMesh& mesh) {
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Error: could not open mesh file '" << filename << "'\n";
        return false;
    }

    std::string header;
    if (!std::getline(in, header)) {
        std::cerr << "Error in mesh file '" << filename << "': file is empty (missing '# FVMESH <major>.<minor>' header)\n";
        return false;
    }
    if (!header.empty() && header.back() == '\r') header.pop_back();

    // Header is "# FVMESH <major>.<minor>" -- parsed (not exact-string-matched)
    // so an unrecognized future version fails with a clear message, the same
    // way read_gmsh dispatches on its $MeshFormat version.
    std::istringstream header_iss(header);
    std::string hash_token, fvmesh_token, version_token;
    header_iss >> hash_token >> fvmesh_token >> version_token;

    int major = -1, minor = -1;
    size_t dot = version_token.find('.');
    if (hash_token == "#" && fvmesh_token == "FVMESH" && dot != std::string::npos) {
        try {
            major = std::stoi(version_token.substr(0, dot));
            minor = std::stoi(version_token.substr(dot + 1));
        } catch (...) {
            major = -1;
        }
    }

    if (major != 1 || minor != 0) {
        std::cerr << "Error in mesh file '" << filename << "': unsupported FVMESH version '"
                   << version_token << "' (supported: 1.0)\n";
        return false;
    }

    std::vector<TempElement> cell_elements;    // CELLS -> mesh.cells (any vertex count)
    std::vector<TempElement> line_elements;    // BOUNDARY -> boundary edges
    std::map<int, std::string> physical_names; // synthesized patch id -> patch_name
    std::unordered_map<std::string, int> patch_name_to_id;

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string tag;
        iss >> tag;

        if (tag == "POINTS") {
            int count = 0;
            iss >> count;
            mesh.nodes.resize(count);
            for (int i = 0; i < count; ++i) {
                double x, y;
                in >> x >> y;
                mesh.nodes[i] = {x, y};
            }
            std::getline(in, line); // consume rest of the last point line
        } else if (tag == "CELLS") {
            int count = 0;
            iss >> count;
            for (int i = 0; i < count; ++i) {
                int n_verts = 0;
                in >> n_verts;

                TempElement elem;
                elem.type = 0;
                elem.physical_id = -1;
                elem.node_ids.resize(n_verts);
                for (int k = 0; k < n_verts; ++k) {
                    int idx;
                    in >> idx;
                    if (idx < 0 || (size_t)idx >= mesh.nodes.size()) {
                        std::cerr << "Error in mesh file '" << filename
                                  << "': CELLS entry references out-of-range node index " << idx << "\n";
                        return false;
                    }
                    elem.node_ids[k] = idx;
                }
                cell_elements.push_back(std::move(elem));
            }
            std::getline(in, line); // consume rest of the last cell line
        } else if (tag == "BOUNDARY") {
            int count = 0;
            iss >> count;
            for (int i = 0; i < count; ++i) {
                int a, b;
                std::string patch_name;
                in >> a >> b >> patch_name;
                if (a < 0 || (size_t)a >= mesh.nodes.size() || b < 0 || (size_t)b >= mesh.nodes.size()) {
                    std::cerr << "Error in mesh file '" << filename
                              << "': BOUNDARY entry references out-of-range node index\n";
                    return false;
                }

                auto id_it = patch_name_to_id.find(patch_name);
                int physical_id;
                if (id_it == patch_name_to_id.end()) {
                    physical_id = (int)patch_name_to_id.size();
                    patch_name_to_id[patch_name] = physical_id;
                    physical_names[physical_id] = patch_name;
                } else {
                    physical_id = id_it->second;
                }

                TempElement elem;
                elem.type = 0;
                elem.physical_id = physical_id;
                elem.node_ids = {a, b};
                line_elements.push_back(std::move(elem));
            }
            std::getline(in, line); // consume rest of the last boundary line
        } else {
            std::cerr << "Error in mesh file '" << filename << "': unrecognized section '" << tag << "'\n";
            return false;
        }
    }

    if (mesh.nodes.empty() || cell_elements.empty()) {
        std::cerr << "Error in mesh file '" << filename << "': no POINTS or CELLS sections found\n";
        return false;
    }

    return build_cells_faces_patches(mesh, cell_elements, line_elements, physical_names);
}

// See MeshReader.h for the input/output contract. Dispatches purely on
// 'filename's extension -- no content sniffing, unlike read_gmsh's internal
// version detection, since the two formats' extensions are unambiguous.
bool MeshReader::read(const std::string& filename, UnstructuredMesh& mesh) {
    auto has_extension = [&](const std::string& ext) {
        return filename.size() >= ext.size() &&
               filename.compare(filename.size() - ext.size(), ext.size(), ext) == 0;
    };

    if (has_extension(".msh")) {
        return read_gmsh(filename, mesh);
    }
    if (has_extension(".fvmesh")) {
        return read_fvmesh(filename, mesh);
    }

    std::cerr << "Error in mesh file '" << filename
              << "': unrecognized extension (expected .msh or .fvmesh)\n";
    return false;
}
