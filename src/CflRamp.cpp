// SPDX-License-Identifier: GPL-3.0-only
#include "CflRamp.h"

#include <algorithm>
#include <cmath>

double next_cfl(double current_cfl, double cfl_max, const std::vector<double>& residual_history,
                 const CflRampParams& params) {
    if (static_cast<int>(residual_history.size()) < params.window) {
        return current_cfl;
    }

    double trend = std::log10(residual_history.back()) - std::log10(residual_history.front());

    const double hold_deadband = std::log10(1.1); // ~10% ratio, same tolerance as a flat/roughly-converged residual

    double next = current_cfl;
    if (trend > params.divergence_threshold) {
        next = params.cfl_min;
    } else if (trend < 0.0) {
        next = current_cfl * params.growth_factor;
    } else if (trend > hold_deadband) {
        next = current_cfl * params.shrink_factor;
    }

    return std::clamp(next, params.cfl_min, cfl_max);
}
