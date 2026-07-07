// SPDX-License-Identifier: GPL-3.0-only
#ifndef CFLRAMP_H_INCLUDED
#define CFLRAMP_H_INCLUDED

#include <vector>

// Selects between a fixed CFL for the whole run (Fixed, today's behavior) or
// a residual-based ramp that starts at a low CFL and grows toward the case
// file's existing "cfl" ceiling as the run converges (Ramp) -- see
// docs/adaptive-cfl-ramp-plan.md for the design discussion, including the
// SU2 (CSolver::AdaptCFLNumber) precedent this windowed-trend rule adapts
// from an implicit, per-point solver to this project's explicit, single
// global-cfl solvers.
enum class CflMode { Fixed, Ramp };

struct CflRampParams {
    double cfl_min = 0.05;              // starting/floor CFL (case-file "cfl_min")
    double growth_factor = 1.5;         // multiplier applied when the residual trend is decreasing ("cfl_ramp_growth")
    double shrink_factor = 0.5;         // multiplier applied when the residual trend is mildly rising ("cfl_ramp_shrink")
    int window = 15;                    // steps of residual history the trend is computed over ("cfl_ramp_window")
    double divergence_threshold = 2.0;  // log10(residual) rise over the window that forces a hard reset to cfl_min ("cfl_ramp_divergence_threshold")
};

// Stateless: given the current cfl, the ceiling cfl_max (the case file's
// existing "cfl" value), a rolling window of the most recent step residual
// scalars (oldest first, newest = residual_history.back(), already including
// the residual from the step just completed), and the ramp parameters,
// returns the next cfl to use for the following step, clamped to
// [params.cfl_min, cfl_max].
//
// The trend is the change in log10(residual) from the oldest to the newest
// entry in residual_history:
//   - fewer than params.window entries: not enough history yet, hold
//     current_cfl unchanged (this is why a ramp run spends its first
//     params.window steps at cfl_min rather than growing immediately --
//     deliberate, so the very noisy initial-transient steps can't be
//     mistaken for a real decreasing trend);
//   - trend > params.divergence_threshold: hard reset to params.cfl_min
//     (SU2's "totalChange > 2.0" divergence path);
//   - trend < 0: current_cfl *= params.growth_factor, clamped to cfl_max;
//   - trend within [0, params.divergence_threshold]: current_cfl *=
//     params.shrink_factor, clamped to params.cfl_min, UNLESS trend is
//     within a small deadband of 0 (log10(1.1), i.e. an 10% ratio) in which
//     case current_cfl is held unchanged -- avoids reacting to a flat,
//     roughly-converged residual as if it were a real rise.
// Input:   current_cfl      - the cfl used for the step just completed
//          cfl_max          - the ceiling to clamp growth to (case file's "cfl")
//          residual_history - rolling window of raw (non-log) residual scalars,
//                              oldest first; caller owns this buffer (push the
//                              new residual, evict the oldest once it exceeds
//                              params.window entries)
//          params           - see CflRampParams above
// Returns: the cfl to use for the next step
double next_cfl(double current_cfl, double cfl_max, const std::vector<double>& residual_history,
                 const CflRampParams& params);

#endif // CFLRAMP_H_INCLUDED
