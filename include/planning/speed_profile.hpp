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

// =============================================================================
// Config (single ST -> DP -> Tunnel -> QP method)
// =============================================================================

struct SpeedProfileConfig
{
  // horizon
  double total_time = 8.0;   // [s]
  double dt_dp      = 0.5;   // [s] coarse DP time step
  double dt_qp      = 0.1;   // [s] QP + output time step
  double s_horizon  = 120.0; // [m] forward horizon relative to s0

  // SL projection helpers (used for projecting ego/obstacles to reference line)
  double projection_window = 40.0; // [m] around seed

  // obstacle sampling + inflation (for ST boundary building)
  double obstacle_dt                  = 0.2; // [s]
  double obstacle_longitudinal_buffer = 3.0; // [m]
  double obstacle_lateral_buffer      = 0.0; // [m]

  // DP grid resolution
  double ds_dp = 1.0; // [m]

  // hard kinematic limits
  double v_min = 0.0;  // [m/s]
  double v_max = 25.0; // [m/s]
  double a_min = -6.0; // [m/s^2]
  double a_max = 3.0;  // [m/s^2]
  double j_max = 6.0;  // [m/s^3]

  // tunnel around DP solution (convex corridor)
  double tunnel_back   = 2.0;  // [m]
  double tunnel_front  = 6.0;  // [m]
  double tunnel_expand = 20.0; // [m] extra slack when carving/merging

  // DP costs (starter knobs)
  double w_dp_progress        = 10.0;
  double w_dp_speed           = 1.0;
  double w_dp_accel           = 1.0;
  double w_dp_jerk            = 1.0;
  double w_dp_obstacle        = 10.0;
  double w_dp_limit_violation = 100.0;


  // QP smoothing costs (starter knobs)
  double w_qp_track_dp = 10.0;
  double w_qp_track_s  = 0.0;

  double w_qp_accel           = 100.0;
  double w_qp_jerk            = 1.0;
  double w_qp_v0              = 5.0;
  double w_qp_a0              = 2.0;
  double w_qp_limit_violation = 100.0;

  // speed limit model sampling
  double speed_limit_ds               = 1.0;  // [m]
  double default_legal_speed          = 30.0; // [m/s] fallback if no map limit and no SpeedLimits samples
  double speed_fraction_of_limit      = 1.0;  // e.g. 0.95
  double max_lateral_acc              = 2.0;  // [m/s^2] comfort
  double max_curvature                = 1.0;  // [1/m] clamp
  double curvature_eps                = 1e-4;
  int    curvature_smoothing_half_win = 5;

  // optional: preferred cruise target (if NaN, planner uses v_max / limit)
  double v_cruise = std::numeric_limits<double>::quiet_NaN();

  // keep 1 unless you add re-linearization for speed bound terms
  int qp_relinearize_iters = 1;
};

// Plans a longitudinal speed profile s(t) along DrivableArea::reference_line.
// Method:
//  1) Build ST boundaries from TrafficParticipant trajectories / fallback extrapolation
//  2) DP on a discretized ST graph -> s_dp(t)
//  3) Build a convex tunnel around s_dp excluding obstacles
//  4) QP smoothing in the tunnel -> final SpeedProfile
//
// Convention: returns an empty vector on failure.
SpeedProfile plan_speed_profile( const DrivableArea& area, const dynamics::TrafficParticipantSet& participants,
                                 const dynamics::VehicleStateDynamic& ego_state, const SpeedProfileConfig& config );

} // namespace adore::planner
