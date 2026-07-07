// SPDX-License-Identifier: GPL-3.0-only
#ifndef CHECKPOINT_H_INCLUDED
#define CHECKPOINT_H_INCLUDED

#include <string>
#include <vector>

#include "EulerState.h"

// Which equation set a checkpoint was written for. A standalone enum (not
// CaseInput's EquationSet) so this module doesn't need to depend on
// CaseInput.h; callers convert between the two at the call site.
enum class CheckpointEquation : int {
    Diffusion = 0,
    Euler = 1,
    AdvectionDiffusion = 2,
    NavierStokes = 3,
    RANS_SA = 4,
    RANS_SST = 5
};

// Saves/restores solver state (absolute step index + field values) to a
// simple binary file, so a run that stops (step limit, residual convergence,
// or a user interrupt) can be resumed later from the same case file.
//
// File layout (fixed-width types, native byte order -- not intended to be
// portable across architectures, only round-tripped on the same machine):
//   offset  size  field
//   0       8     magic bytes "FVMCKPT\0"
//   8       4     format version (uint32_t, currently 2)
//   12      4     equation tag (uint32_t): 0 = Diffusion, 1 = Euler, 2 = AdvectionDiffusion, 3 = NavierStokes, 4 = RANS_SA, 5 = RANS_SST
//   16      8     cell_count (uint64_t)
//   24      8     step_index (int64_t) -- absolute step index reached
//   32      8     build_number (uint64_t) -- FV_BUILD_NUMBER of the solver
//                 build that wrote this checkpoint; informational only, not
//                 used for compatibility gating (the format version is)
//   40      ...   payload: raw doubles, no padding
//                 Diffusion / AdvectionDiffusion: cell_count doubles (phi)
//                 Euler / NavierStokes:           cell_count EulerState structs (rho,rho_u,rho_v,E per cell)
//                 RANS_SA:                         cell_count EulerState structs, immediately followed by
//                                                   cell_count doubles (nut) -- RANSTurbulenceSASolver's mean-flow
//                                                   state plus its one extra transported scalar
//                 RANS_SST:                        cell_count EulerState structs, immediately followed by
//                                                   cell_count doubles (k), immediately followed by cell_count
//                                                   doubles (omega) -- RANSTurbulenceSSTSolver's mean-flow state
//                                                   plus its two extra transported scalars
//
// Format version 2 added the build_number field; a version-1 checkpoint
// (written before this field existed) is rejected by read() as an
// unsupported version, same as any other format mismatch -- there is no
// migration path, since this is meant to round-trip a single in-progress run,
// not archive results long-term. RANS_SA's dual-payload layout above did NOT
// require a version bump: the fixed 40-byte header is unchanged, and the
// payload shape has always been a function of the equation tag alone.
// RANS_SST's triple-payload layout above did not require one either, for the
// identical reason.
namespace Checkpoint {
    // Input:  filename - path to probe
    // Returns: true if a file exists at that path (used to decide whether to
    //          auto-resume); does not validate its contents
    bool exists(const std::string& filename);

    // Input:
    //   filename     - path to write to
    //   equation     - which solver this checkpoint is for
    //   step_index   - absolute step index reached
    //   build_number - the solver build that produced this checkpoint
    //                  (FV_BUILD_NUMBER), stored purely as metadata
    //   phi/U        - the field to save, one entry per mesh cell
    // Output: none (writes to disk)
    // Returns: true if the file was written successfully; false otherwise
    bool write(const std::string& filename, CheckpointEquation equation, long long step_index,
               unsigned long long build_number, const std::vector<double>& phi);
    bool write(const std::string& filename, CheckpointEquation equation, long long step_index,
               unsigned long long build_number, const std::vector<EulerState>& U);

    // RANS_SA overload: writes U immediately followed by nut (see the dual-payload
    // layout above). U.size() and nut.size() must match (both mesh.cells.size()).
    bool write(const std::string& filename, CheckpointEquation equation, long long step_index,
               unsigned long long build_number, const std::vector<EulerState>& U, const std::vector<double>& nut);

    // RANS_SST overload: writes U immediately followed by k, immediately
    // followed by omega (see the triple-payload layout above). U.size(),
    // k.size(), and omega.size() must all match (mesh.cells.size()).
    bool write(const std::string& filename, CheckpointEquation equation, long long step_index,
               unsigned long long build_number, const std::vector<EulerState>& U, const std::vector<double>& k,
               const std::vector<double>& omega);

    // Input:
    //   filename            - path to read from
    //   expected_equation   - the checkpoint must have been written for this equation set
    //   expected_cell_count - the checkpoint must have this many cells (i.e.
    //                         match the currently-loaded mesh)
    // Output:
    //   out_step_index   - the absolute step index the checkpoint reached
    //   out_build_number - the solver build that produced this checkpoint
    //   out_phi/out_U    - the restored field, one entry per mesh cell
    // Returns: true on success; false if the file could not be opened, its
    //          magic/version/equation/cell_count don't match, or it's
    //          truncated (a descriptive message is printed to stderr in
    //          every failure case except a plain open failure)
    bool read(const std::string& filename, CheckpointEquation expected_equation, size_t expected_cell_count,
              long long& out_step_index, unsigned long long& out_build_number, std::vector<double>& out_phi);
    bool read(const std::string& filename, CheckpointEquation expected_equation, size_t expected_cell_count,
              long long& out_step_index, unsigned long long& out_build_number, std::vector<EulerState>& out_U);

    // RANS_SA overload: reads U immediately followed by nut (see the dual-payload
    // layout above). Output: out_U, out_nut - both resized to expected_cell_count.
    bool read(const std::string& filename, CheckpointEquation expected_equation, size_t expected_cell_count,
              long long& out_step_index, unsigned long long& out_build_number, std::vector<EulerState>& out_U,
              std::vector<double>& out_nut);

    // RANS_SST overload: reads U immediately followed by k, immediately
    // followed by omega (see the triple-payload layout above). Output:
    // out_U, out_k, out_omega - all resized to expected_cell_count.
    bool read(const std::string& filename, CheckpointEquation expected_equation, size_t expected_cell_count,
              long long& out_step_index, unsigned long long& out_build_number, std::vector<EulerState>& out_U,
              std::vector<double>& out_k, std::vector<double>& out_omega);
}

#endif // CHECKPOINT_H_INCLUDED
