/********************************************************************************
 * Copyright (c) 2025 Contributors to the Eclipse Foundation
 *
 * See the NOTICE file(s) distributed with this work for additional
 * information regarding copyright ownership.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * https://www.eclipse.org/legal/epl-2.0
 *
 * SPDX-License-Identifier: EPL-2.0
 ********************************************************************************/

#include <cmath>
#include <gtest/gtest.h>

#include <utility>

#include "adore_math/angles.h"
#include "adore_math/point.h"

#include "dynamics/traffic_participant.hpp"
#include "dynamics/vehicle_state.hpp"
#include "planning/common/planning_helpers.hpp"

using adore::dynamics::VehicleStateDynamic;

namespace
{

constexpr double kEps = 1e-9;

} // namespace

TEST( PlanningHelpers, ErrorsStraightAheadZeroLateral )
{
  VehicleStateDynamic state;
  state.x         = 0.0;
  state.y         = 0.0;
  state.yaw_angle = 0.0; // pointing along +x

  const double target_x   = 10.0;
  const double target_y   = 0.0;
  const double target_yaw = 0.0;

  const auto [lateral_error, heading_error] = adore::planner::calculate_errors( state, target_x, target_y, target_yaw );

  EXPECT_NEAR( lateral_error, 0.0, kEps );
  EXPECT_NEAR( heading_error, 0.0, kEps );
}

TEST( PlanningHelpers, TargetToLeftPositiveLateralError )
{
  VehicleStateDynamic state;
  state.x         = 0.0;
  state.y         = 0.0;
  state.yaw_angle = 0.0; // pointing along +x

  // 1 m to the left of the vehicle's heading (positive y)
  const double target_x   = 0.0;
  const double target_y   = 1.0;
  const double target_yaw = 0.0;

  const auto [lateral_error, heading_error] = adore::planner::calculate_errors( state, target_x, target_y, target_yaw );

  EXPECT_NEAR( lateral_error, 1.0, 1e-9 );
  EXPECT_NEAR( heading_error, 0.0, kEps );
}

TEST( PlanningHelpers, TargetToRightNegativeLateralErrorInRotatedFrame )
{
  VehicleStateDynamic state;
  state.x = 0.0;
  state.y = 0.0;
  // Vehicle pointing along +y axis (90 degrees)
  const double yaw = M_PI_2;
  state.yaw_angle  = yaw;

  // Place the target at (1, 0): this is to the "right" of +y direction
  const double target_x   = 1.0;
  const double target_y   = 0.0;
  const double target_yaw = yaw;

  const auto [lateral_error, heading_error] = adore::planner::calculate_errors( state, target_x, target_y, target_yaw );

  // For yaw = pi/2, dx = 1, dy = 0:
  // lateral_error = dy*cos(yaw) - dx*sin(yaw) = 0 - 1*1 = -1
  EXPECT_NEAR( lateral_error, -1.0, 1e-9 );
  EXPECT_NEAR( heading_error, 0.0, kEps );
}
