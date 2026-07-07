// SPDX-License-Identifier: GPL-3.0-only
#ifndef MESHCHECKREPORT_H_INCLUDED
#define MESHCHECKREPORT_H_INCLUDED

#include <map>
#include <string>
#include <vector>

// One check-criterion violation found during --validate-mesh's check. Kept
// deliberately terse (name + value, no prose) -- see
// docs/mesh-validation-criteria.md for what each 'criterion' name means, its
// formula, and its pass/fail threshold.
struct MeshCheckIssue {
    std::string severity;  // "error" or "warning"
    std::string entity;    // "node", "cell", or "face"
    int id;                 // index into UnstructuredMesh::nodes/cells/faces
    std::string criterion; // see docs/mesh-validation-criteria.md
    double value;           // the measured value that triggered this issue
};

// Progress/result tracking for --validate-mesh's 3 fixed steps (Parsing,
// Building connectivity, Validating geometry). Optionally threaded as a
// nullable pointer through MeshReader's parsing functions so they can record
// step completion/issues without changing behavior for real solver runs,
// which never construct one and always pass nullptr.
//
// When a JSON report path is set, the full report is rewritten to that file
// after every step completes (or fails), so a parent process monitoring the
// file sees incremental progress; a one-line-per-step stdout summary is
// always printed regardless of whether a JSON path is set.
class MeshCheckReport {
public:
    // Input:  mesh_file - path being checked, recorded in the report
    //         json_path - if non-empty, the report is (re)written to this
    //                      path after construction and after every
    //                      complete_step()/fail_step() call; if empty, no
    //                      file is written (stdout-only)
    // Output: an initial JSON snapshot (all 3 steps "pending") is written
    //         immediately if json_path is non-empty, so a parent process
    //         sees the file exist before the first step finishes
    MeshCheckReport(std::string mesh_file, std::string json_path);

    // Prints the one-line intro naming all 3 steps up front.
    void print_intro() const;

    // Records one criterion violation against step 'step_index' (0=Parsing,
    // 1=Building connectivity, 2=Validating geometry), to be included in
    // that step's "issues" list once complete_step()/fail_step() is called
    // for it.
    void add_issue(int step_index, const MeshCheckIssue& issue);

    // Marks step 'step_index' completed with the given per-step counts.
    // 'passed' should reflect whether any 'error'-severity issue was
    // recorded against this step (warnings alone still count as a pass).
    // Prints the step's one-line stdout summary (plus one line per issue)
    // and rewrites the JSON report file if a path was given. A no-op if
    // this step was already resolved (completed/failed) by an earlier call
    // -- see MeshReader::read's fallback use of fail_step().
    void complete_step(int step_index, const std::map<std::string, long long>& counts, bool passed);

    // Marks step 'step_index' as failed and every later step as "skipped"
    // (they depend on this step's output, so they never ran). Prints/writes
    // like complete_step(). A no-op if this step was already resolved.
    void fail_step(int step_index, const std::map<std::string, long long>& counts);

    // Input:  none
    // Returns: true if every step resolved with status "completed" and no
    //          'error'-severity issue was recorded against any step
    bool passed() const;

private:
    struct Step {
        std::string name;
        std::string status = "pending"; // pending | completed | failed | skipped
        std::string result;             // "" | pass | fail
        std::map<std::string, long long> counts;
        std::vector<MeshCheckIssue> issues;
    };

    std::string mesh_file;
    std::string json_path;
    std::vector<Step> steps; // fixed size 3, see constructor

    void print_step_line(int step_index) const;
    void write_json() const;
};

#endif // MESHCHECKREPORT_H_INCLUDED
