// SPDX-License-Identifier: GPL-3.0-only
#include "Checkpoint.h"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>

namespace {

const char MAGIC[8] = {'F', 'V', 'M', 'C', 'K', 'P', 'T', '\0'};
const uint32_t FORMAT_VERSION = 2;

// Fixed-size checkpoint header; see Checkpoint.h for the byte layout this
// mirrors. Deliberately laid out so no compiler padding is inserted (every
// member is already aligned to its own size at its offset), verified by the
// static_assert below. build_number is a uint64_t (rather than a narrower
// type that would fit its actual range) specifically so it lands on an
// 8-byte boundary and the struct's overall size stays a multiple of its own
// 8-byte alignment, with no trailing padding.
struct Header {
    char magic[8];
    uint32_t version;
    uint32_t equation;
    uint64_t cell_count;
    int64_t step_index;
    uint64_t build_number;
};
static_assert(sizeof(Header) == 40, "Checkpoint::Header must not have compiler-inserted padding");

// Writes the fixed-size header only (no payload) -- shared by every write()
// overload, including the RANS one below whose payload is two back-to-back
// buffers rather than one.
// Input:  out          - open binary output stream
//         equation     - equation tag to embed
//         step_index   - absolute step index reached
//         build_number - solver build number to embed, informational only
//         cell_count   - number of cells the payload covers
// Output: the header is appended to 'out'
void write_header(std::ofstream& out, CheckpointEquation equation, long long step_index,
                   unsigned long long build_number, uint64_t cell_count) {
    Header header;
    std::memcpy(header.magic, MAGIC, sizeof(MAGIC));
    header.version = FORMAT_VERSION;
    header.equation = static_cast<uint32_t>(equation);
    header.cell_count = cell_count;
    header.step_index = static_cast<int64_t>(step_index);
    header.build_number = static_cast<uint64_t>(build_number);

    out.write(reinterpret_cast<const char*>(&header), sizeof(header));
}

// Writes the fixed-size header followed by the raw payload bytes.
// Input:  out          - open binary output stream
//         equation     - equation tag to embed
//         step_index   - absolute step index reached
//         build_number - solver build number to embed, informational only
//         cell_count   - number of cells the payload covers
//         payload      - pointer to the raw field data
//         payload_size - size of the payload, in bytes
// Output: header + payload are appended to 'out'
void write_header_and_payload(std::ofstream& out, CheckpointEquation equation, long long step_index,
                               unsigned long long build_number, uint64_t cell_count, const void* payload,
                               size_t payload_size) {
    write_header(out, equation, step_index, build_number, cell_count);
    out.write(reinterpret_cast<const char*>(payload), payload_size);
}

// Reads and validates the fixed-size header against the caller's expectations.
// Input:  in                  - open binary input stream
//         filename            - path, for diagnostics only
//         expected_equation   - required equation tag
//         expected_cell_count - required cell count
// Output: out_step_index   - the header's step_index, if validation succeeds
//         out_build_number - the header's build_number, if validation succeeds
// Returns: true if the header is well-formed and matches expectations;
//          false otherwise (a descriptive message is printed to stderr)
bool read_and_validate_header(std::ifstream& in, const std::string& filename,
                               CheckpointEquation expected_equation, size_t expected_cell_count,
                               long long& out_step_index, unsigned long long& out_build_number) {
    Header header;
    in.read(reinterpret_cast<char*>(&header), sizeof(header));
    if (!in) {
        std::cerr << "Checkpoint file '" << filename << "' is truncated or corrupt (short header)\n";
        return false;
    }
    if (std::memcmp(header.magic, MAGIC, sizeof(MAGIC)) != 0) {
        std::cerr << "Checkpoint file '" << filename << "' is not a valid FiniteVolume checkpoint (bad magic)\n";
        return false;
    }
    if (header.version != FORMAT_VERSION) {
        std::cerr << "Checkpoint file '" << filename << "' has unsupported version " << header.version
                   << " (expected " << FORMAT_VERSION << ")\n";
        return false;
    }
    if (header.equation != static_cast<uint32_t>(expected_equation)) {
        std::cerr << "Checkpoint file '" << filename
                   << "' was written for a different equation set than the current run\n";
        return false;
    }
    if (header.cell_count != expected_cell_count) {
        std::cerr << "Checkpoint file '" << filename << "' has " << header.cell_count
                   << " cells but the current mesh has " << expected_cell_count << " cells (mesh changed?)\n";
        return false;
    }

    out_step_index = header.step_index;
    out_build_number = header.build_number;
    return true;
}

} // namespace

bool Checkpoint::exists(const std::string& filename) {
    std::ifstream in(filename);
    return in.is_open();
}

bool Checkpoint::write(const std::string& filename, CheckpointEquation equation, long long step_index,
                        unsigned long long build_number, const std::vector<double>& phi) {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        return false;
    }
    write_header_and_payload(out, equation, step_index, build_number, phi.size(), phi.data(),
                              phi.size() * sizeof(double));
    return static_cast<bool>(out);
}

bool Checkpoint::write(const std::string& filename, CheckpointEquation equation, long long step_index,
                        unsigned long long build_number, const std::vector<EulerState>& U) {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        return false;
    }
    write_header_and_payload(out, equation, step_index, build_number, U.size(), U.data(),
                              U.size() * sizeof(EulerState));
    return static_cast<bool>(out);
}

bool Checkpoint::write(const std::string& filename, CheckpointEquation equation, long long step_index,
                        unsigned long long build_number, const std::vector<EulerState>& U,
                        const std::vector<double>& nut) {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        return false;
    }
    write_header(out, equation, step_index, build_number, U.size());
    out.write(reinterpret_cast<const char*>(U.data()), U.size() * sizeof(EulerState));
    out.write(reinterpret_cast<const char*>(nut.data()), nut.size() * sizeof(double));
    return static_cast<bool>(out);
}

bool Checkpoint::write(const std::string& filename, CheckpointEquation equation, long long step_index,
                        unsigned long long build_number, const std::vector<EulerState>& U,
                        const std::vector<double>& k, const std::vector<double>& omega) {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        return false;
    }
    write_header(out, equation, step_index, build_number, U.size());
    out.write(reinterpret_cast<const char*>(U.data()), U.size() * sizeof(EulerState));
    out.write(reinterpret_cast<const char*>(k.data()), k.size() * sizeof(double));
    out.write(reinterpret_cast<const char*>(omega.data()), omega.size() * sizeof(double));
    return static_cast<bool>(out);
}

bool Checkpoint::read(const std::string& filename, CheckpointEquation expected_equation, size_t expected_cell_count,
                       long long& out_step_index, unsigned long long& out_build_number, std::vector<double>& out_phi) {
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        return false;
    }
    if (!read_and_validate_header(in, filename, expected_equation, expected_cell_count, out_step_index, out_build_number)) {
        return false;
    }

    out_phi.resize(expected_cell_count);
    in.read(reinterpret_cast<char*>(out_phi.data()), out_phi.size() * sizeof(double));
    if (!in) {
        std::cerr << "Checkpoint file '" << filename << "' is truncated or corrupt (short payload)\n";
        return false;
    }
    return true;
}

bool Checkpoint::read(const std::string& filename, CheckpointEquation expected_equation, size_t expected_cell_count,
                       long long& out_step_index, unsigned long long& out_build_number, std::vector<EulerState>& out_U) {
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        return false;
    }
    if (!read_and_validate_header(in, filename, expected_equation, expected_cell_count, out_step_index, out_build_number)) {
        return false;
    }

    out_U.resize(expected_cell_count);
    in.read(reinterpret_cast<char*>(out_U.data()), out_U.size() * sizeof(EulerState));
    if (!in) {
        std::cerr << "Checkpoint file '" << filename << "' is truncated or corrupt (short payload)\n";
        return false;
    }
    return true;
}

bool Checkpoint::read(const std::string& filename, CheckpointEquation expected_equation, size_t expected_cell_count,
                       long long& out_step_index, unsigned long long& out_build_number, std::vector<EulerState>& out_U,
                       std::vector<double>& out_nut) {
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        return false;
    }
    if (!read_and_validate_header(in, filename, expected_equation, expected_cell_count, out_step_index, out_build_number)) {
        return false;
    }

    out_U.resize(expected_cell_count);
    in.read(reinterpret_cast<char*>(out_U.data()), out_U.size() * sizeof(EulerState));
    if (!in) {
        std::cerr << "Checkpoint file '" << filename << "' is truncated or corrupt (short payload)\n";
        return false;
    }

    out_nut.resize(expected_cell_count);
    in.read(reinterpret_cast<char*>(out_nut.data()), out_nut.size() * sizeof(double));
    if (!in) {
        std::cerr << "Checkpoint file '" << filename << "' is truncated or corrupt (short payload)\n";
        return false;
    }
    return true;
}

bool Checkpoint::read(const std::string& filename, CheckpointEquation expected_equation, size_t expected_cell_count,
                       long long& out_step_index, unsigned long long& out_build_number, std::vector<EulerState>& out_U,
                       std::vector<double>& out_k, std::vector<double>& out_omega) {
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        return false;
    }
    if (!read_and_validate_header(in, filename, expected_equation, expected_cell_count, out_step_index, out_build_number)) {
        return false;
    }

    out_U.resize(expected_cell_count);
    in.read(reinterpret_cast<char*>(out_U.data()), out_U.size() * sizeof(EulerState));
    if (!in) {
        std::cerr << "Checkpoint file '" << filename << "' is truncated or corrupt (short payload)\n";
        return false;
    }

    out_k.resize(expected_cell_count);
    in.read(reinterpret_cast<char*>(out_k.data()), out_k.size() * sizeof(double));
    if (!in) {
        std::cerr << "Checkpoint file '" << filename << "' is truncated or corrupt (short payload)\n";
        return false;
    }

    out_omega.resize(expected_cell_count);
    in.read(reinterpret_cast<char*>(out_omega.data()), out_omega.size() * sizeof(double));
    if (!in) {
        std::cerr << "Checkpoint file '" << filename << "' is truncated or corrupt (short payload)\n";
        return false;
    }
    return true;
}
