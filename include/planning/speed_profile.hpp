/********************************************************************************
 * Copyright (c) 2025 Contributors
 *
 * SPDX-License-Identifier: EPL-2.0
 ********************************************************************************/
#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "dynamics/comfort_settings.hpp"
#include "dynamics/traffic_participant.hpp"
#include "dynamics/vehicle_state.hpp"
#include "planning/drivable_area.hpp"

namespace adore::planner
{

// =============================================================================
// Speed profile output (ST: s(t) with derivatives)
// =============================================================================

struct SpeedProfilePoint
{
  double t = 0.0; // [s] seconds from plan start
  double s = 0.0; // [m] station along DrivableArea::reference_line
  double v = 0.0; // [m/s] ds/dt
  double a = 0.0; // [m/s^2] dv/dt
  double j = 0.0; // [m/s^3] da/dt

  // Optional annotations (filled by planner; useful for debugging/inspection)
  double kappa   = std::numeric_limits<double>::quiet_NaN(); // [1/m] curvature at s
  double v_limit = std::numeric_limits<double>::quiet_NaN(); // [m/s] composed speed limit at s
};

using SpeedProfile = std::vector<SpeedProfilePoint>;

// (Removed redundant Config block)

// Plans a longitudinal speed profile s(t) along DrivableArea::reference_line.
// Method:
//  1) Build ST boundaries from TrafficParticipant trajectories / fallback extrapolation
//  2) DP on a discretized ST graph -> s_dp(t)
//  3) Build a convex tunnel around s_dp excluding obstacles
//  4) QP smoothing in the tunnel -> final SpeedProfile
//
// Convention: returns an empty vector on failure.
SpeedProfile plan_speed_profile( const DrivableArea& area, const dynamics::TrafficParticipantSet& participants,
                                 const dynamics::VehicleStateDynamic& ego_state, const dynamics::PhysicalVehicleParameters& vehicle_params,
                                 const dynamics::ComfortSettings& comfort_settings, const SpeedProfileConfig& config );

// Smooths the new profile by "sticking" to the previous profile where possible (hysteresis).
void apply_temporal_smoothing( SpeedProfile& new_profile, const SpeedProfile& previous_profile, const SpeedProfileConfig& config );

} // namespace adore::planner
