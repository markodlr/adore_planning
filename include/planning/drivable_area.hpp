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
#include "adore_math/geometry/projection.hpp"

#include "dynamics/traffic_participant.hpp"
#include "planning/common/planning_config.hpp"

namespace adore
{
namespace planner
{

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
                                   const dynamics::TrafficParticipantSet&     traffic_participants,
                                   const dynamics::PhysicalVehicleParameters& vehicle_params, const DrivableAreaConfig& config = {} );

void get_obstacle_envelope( const adore::map::Route& route, double start_s, double end_s, const dynamics::TrafficParticipant& participant,
                            const dynamics::PhysicalVehicleParameters& vehicle_params, const DrivableAreaConfig& config,
                            std::map<double, std::pair<double, double>>& envelope );

std::vector<DrivableArea> generate_candidates( const DrivableArea& initial_area, const adore::map::Route& route, double start_s,
                                               double end_s, const dynamics::TrafficParticipantSet& traffic_participants,
                                               const DrivableAreaConfig& config );

double score_candidate( const DrivableArea& area, const adore::map::Route& route );

/**
 * @brief Check if any moving traffic is detected in oncoming lanes.
 *
 * @param route The route being followed
 * @param start_s Start station of the planning horizon
 * @param end_s End station of the planning horizon
 * @param participants Traffic participants to check
 * @return true if moving traffic is found in lanes with opposite direction to the route lane
 */
bool has_traffic_in_oncoming_lanes( const adore::map::Route& route, double start_s, double end_s,
                                    const dynamics::TrafficParticipantSet& participants );

} // namespace planner
} // namespace adore
