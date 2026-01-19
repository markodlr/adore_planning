/********************************************************************************
 * Copyright (c) 2025 Contributors
 *
 * SPDX-License-Identifier: EPL-2.0
 ********************************************************************************/

#pragma once

#include <cstdint>

namespace adore
{
namespace planner
{

// =============================================================================
// Drivable Area Config
// =============================================================================

enum class LaneScope : uint8_t
{
  MyLane = 0,
  SameDirection,
  AllLanes
};

struct DrivableAreaConfig
{
  LaneScope lane_scope             = LaneScope::AllLanes;
  double    lateral_inflation      = 0.2;
  double    longitudinal_inflation = 1.0;
};

// =============================================================================
// Speed Profile Config
// =============================================================================

struct SpeedProfileConfig
{
  double total_time        = 5.0; // [s] Total planning horizon
  double s_horizon         = 100.0;
  double ds_dp             = 0.25; // [m] DP station step
  double dt_dp             = 0.5;  // [s] DP time step
  double dt_qp             = 0.1;  // [s] QP time step (dense)
  double projection_window = 50.0; // [m] Window for s-projection lookup

  // Kinematic limits
  double v_min = 0.0;

  double v_cruise = -1.0; // if > 0, penalize derivation from this

  // Obstacle
  double obstacle_dt                  = 0.1; // [s] sampling for ST-boundaries
  double obstacle_longitudinal_buffer = 3.0;
  double obstacle_lateral_buffer      = 1.5;

  // DP weights
  double w_dp_progress           = 10.0;     // reward progress (negative cost)
  double w_dp_speed              = 1.0;      // penalize dv
  double w_dp_accel              = 1.0;      // penalize a^2
  double w_dp_jerk               = 10.0;     // penalize j^2
  double w_dp_limit_violation    = 100000.0; // penalize exceeding v_limit
  double w_dp_obstacle_proximity = 10.0;     // penalize closeness to ST-obstacles

  // Tunnel
  double tunnel_back   = 5.0;  // [m] range behind DP guess
  double tunnel_front  = 10.0; // [m] range in front of DP guess
  double tunnel_expand = 1.0;  // [m] extra expansion if no obstacles


  // Curvature & Limits
  double max_curvature                = 0.2; // [1/m] (approx 5m radius)
  double curvature_eps                = 1e-4;
  int    curvature_smoothing_half_win = 5;
  double speed_limit_ds               = 1.0;  // [m]
  double default_legal_speed          = 30.0; // [m/s] fallback if no map limit and no SpeedLimits samples
  double speed_tunnel_width           = 2.0;  // [m/s] width of tunnel around DP solution (total width 2x)
};

} // namespace planner
} // namespace adore
