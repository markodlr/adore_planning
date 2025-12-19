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

#pragma once

#include <map>

#include "adore_map/route.hpp"

#include "dynamics/traffic_participant.hpp"

namespace adore
{
namespace planner
{


// Which lanes to include when constructing the drivable area.
enum class LaneScope
{
  MyLane,        // only the current lane of the route
  SameDirection, // all lanes on the same side (same driving direction)
  AllLanes       // all lanes on the road (both directions)
};

struct DrivableAreaConfig
{
  double safety_margin           = 1.4; // [m] margin from boundaries for ref line
  double max_lateral_change_rate = 0.2; // [m/m] max |dl/ds| for ref line smoothing

  double longitudinal_inflation = 1.5; // [m] inflate participants along route
  double lateral_inflation      = 0.2; // [m] inflate participants sideways

  LaneScope lane_scope = LaneScope::AllLanes; // which lanes to include
};

struct DrivableArea
{
  // s -> boundary MapPoint
  using Boundary = std::map<double, adore::map::MapPoint>;

  Boundary left_boundary;  // outer left boundary of drivable corridor
  Boundary right_boundary; // outer right boundary of drivable corridor
  Boundary reference_line; // smoothed path inside corridor

  bool
  empty() const noexcept
  {
    return left_boundary.empty() || right_boundary.empty();
  }
};

DrivableArea create_drivable_area( const adore::map::Route& route, double start_s, double end_s,
                                   const dynamics::TrafficParticipantSet& traffic_participants, const DrivableAreaConfig& config = {} );

} // namespace planner
} // namespace adore
