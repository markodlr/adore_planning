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
#include <string>

#include "dynamics/traffic_participant.hpp"
#include "dynamics/trajectory.hpp"
#include "planning/drivable_area.hpp"

namespace adore
{
namespace planner
{

class NominalParticipantPrediction
{
public:

  NominalParticipantPrediction() = default;

  /// Parameters (all optional):
  /// - "dt": integration time step [s]
  /// - "number_of_integration_steps": prediction horizon steps
  /// - "default_speed": used if vx <= 0
  /// - "include_initial_state": 1.0 => push current state as first trajectory sample
  void set_parameters( const std::map<std::string, double>& params );

  /// Fills / overwrites participant.trajectory for all participants that have a valid route.
  /// Participants without a route are left with an empty trajectory.
  void plan_trajectories( dynamics::TrafficParticipantSet& traffic_participant_set ) const;

private:

  double dt                          = 0.1; // [s]
  int    number_of_integration_steps = 50;
  double default_speed               = 5.0; // [m/s]
  bool   include_initial_state       = true;

  static dynamics::VehicleStateDynamic get_current_state( const dynamics::TrafficParticipant& participant );
};

} // namespace planner
} // namespace adore
