// SPDX-License-Identifier: GPL-3.0-only
#include "MeshCheckReport.h"

#include <cstdio>
#include <fstream>
#include <iostream>

namespace {

// Escapes a string for embedding in a JSON string literal. The only strings
// this ever needs to escape are file paths and mesh criterion/entity names,
// but paths on Windows contain backslashes, so this can't skip escaping.
std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned char>(c));
                    out += buf;
                } else {
                    out += c;
                }
        }
    }
    return out;
}

} // namespace

MeshCheckReport::MeshCheckReport(std::string mesh_file_, std::string json_path_)
    : mesh_file(std::move(mesh_file_)), json_path(std::move(json_path_)) {
    steps.resize(3);
    steps[0].name = "Parsing";
    steps[1].name = "Building connectivity";
    steps[2].name = "Validating geometry";

    if (!json_path.empty()) write_json();
}

void MeshCheckReport::print_intro() const {
    std::cout << "Mesh check: " << steps.size() << " steps planned (";
    for (size_t i = 0; i < steps.size(); ++i) {
        std::cout << steps[i].name;
        if (i + 1 < steps.size()) std::cout << ", ";
    }
    std::cout << ")\n";
}

void MeshCheckReport::add_issue(int step_index, const MeshCheckIssue& issue) {
    steps[step_index].issues.push_back(issue);
}

void MeshCheckReport::complete_step(int step_index, const std::map<std::string, long long>& counts, bool passed) {
    if (steps[step_index].status != "pending") return;
    steps[step_index].status = "completed";
    steps[step_index].counts = counts;
    steps[step_index].result = passed ? "pass" : "fail";
    print_step_line(step_index);
    if (!json_path.empty()) write_json();
}

void MeshCheckReport::fail_step(int step_index, const std::map<std::string, long long>& counts) {
    if (steps[step_index].status != "pending") return;
    steps[step_index].status = "failed";
    steps[step_index].counts = counts;
    steps[step_index].result = "fail";
    for (size_t s = static_cast<size_t>(step_index) + 1; s < steps.size(); ++s) {
        if (steps[s].status == "pending") steps[s].status = "skipped";
    }
    print_step_line(step_index);
    if (!json_path.empty()) write_json();
}

bool MeshCheckReport::passed() const {
    for (const Step& s : steps) {
        if (s.status != "completed") return false;
        for (const MeshCheckIssue& issue : s.issues) {
            if (issue.severity == "error") return false;
        }
    }
    return true;
}

void MeshCheckReport::print_step_line(int step_index) const {
    const Step& s = steps[step_index];
    std::cout << "Step " << (step_index + 1) << "/" << steps.size() << ": " << s.name << " -- ";

    if (s.status == "failed") {
        std::cout << "FAILED\n";
    } else {
        std::cout << "completed (";
        bool first = true;
        for (const auto& kv : s.counts) {
            if (!first) std::cout << ", ";
            std::cout << kv.second << " " << kv.first;
            first = false;
        }
        std::cout << ") -- ";

        int error_count = 0, warning_count = 0;
        for (const MeshCheckIssue& issue : s.issues) {
            if (issue.severity == "error") ++error_count;
            else ++warning_count;
        }
        if (error_count > 0) {
            std::cout << "FAIL (" << error_count << " error(s)";
            if (warning_count > 0) std::cout << ", " << warning_count << " warning(s)";
            std::cout << ")\n";
        } else if (warning_count > 0) {
            std::cout << "PASS WITH " << warning_count << " WARNING(S)\n";
        } else {
            std::cout << "PASS\n";
        }
    }

    for (const MeshCheckIssue& issue : s.issues) {
        std::cout << "  " << (issue.severity == "error" ? "FAIL" : "WARN") << " " << issue.entity << " "
                   << issue.id << ": " << issue.criterion << " = " << issue.value << "\n";
    }
}

void MeshCheckReport::write_json() const {
    std::ofstream out(json_path);
    if (!out.is_open()) return; // best-effort; the stdout summary is authoritative either way

    bool any_failed = false, any_pending = false, any_warning = false;
    for (const Step& s : steps) {
        if (s.status == "failed") any_failed = true;
        if (s.status == "pending") any_pending = true;
        for (const MeshCheckIssue& issue : s.issues) {
            if (issue.severity == "warning") any_warning = true;
        }
    }
    std::string overall = any_pending  ? "running"
                           : any_failed ? "failed"
                           : any_warning ? "passed_with_warnings"
                                         : "passed";

    out << "{\n";
    out << "  \"operation\": \"validate-mesh\",\n";
    out << "  \"mesh_file\": \"" << json_escape(mesh_file) << "\",\n";
    out << "  \"total_steps\": " << steps.size() << ",\n";
    out << "  \"overall_status\": \"" << overall << "\",\n";
    out << "  \"steps\": [\n";
    for (size_t i = 0; i < steps.size(); ++i) {
        const Step& s = steps[i];
        out << "    {\n";
        out << "      \"step\": " << (i + 1) << ",\n";
        out << "      \"name\": \"" << json_escape(s.name) << "\",\n";
        out << "      \"status\": \"" << s.status << "\",\n";
        out << "      \"result\": " << (s.result.empty() ? "null" : ("\"" + s.result + "\"")) << ",\n";

        out << "      \"counts\": {";
        bool first = true;
        for (const auto& kv : s.counts) {
            if (!first) out << ", ";
            out << "\"" << json_escape(kv.first) << "\": " << kv.second;
            first = false;
        }
        out << "},\n";

        out << "      \"issues\": [";
        for (size_t j = 0; j < s.issues.size(); ++j) {
            const MeshCheckIssue& issue = s.issues[j];
            out << "\n        {\"severity\": \"" << issue.severity << "\", \"entity\": \"" << issue.entity
                << "\", \"id\": " << issue.id << ", \"criterion\": \"" << json_escape(issue.criterion)
                << "\", \"value\": " << issue.value << "}";
            if (j + 1 < s.issues.size()) out << ",";
        }
        if (!s.issues.empty()) out << "\n      ";
        out << "]\n";

        out << "    }" << (i + 1 < steps.size() ? "," : "") << "\n";
    }
    out << "  ]\n";
    out << "}\n";
}
