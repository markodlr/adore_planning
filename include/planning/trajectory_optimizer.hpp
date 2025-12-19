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

#include <cmath>

#include <chrono>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <vector>

#include "adore_map/map.hpp"
#include "adore_map/route.hpp"
#include "adore_math/angles.h"
#include "adore_math/distance.h"

#include "dynamics/comfort_settings.hpp"
#include "dynamics/traffic_participant.hpp"
#include "dynamics/trajectory.hpp"
#include "multi_agent_solver/multi_agent_solver.hpp"
#include "planning/speed_profile.hpp"

// NEW: driveable-area reference line + ST speed profile
#include "planning/drivable_area.hpp"
#include "planning/speed_profile.hpp"

namespace adore
{
namespace planner
{

class TrajectoryOptimizer
{
public:


  dynamics::Trajectory optimize_trajectory( const dynamics::VehicleStateDynamic& current_state,
                                            const dynamics::Trajectory&          reference_trajectory,
                                            const dynamics::Trajectory&          initial_guess = dynamics::Trajectory() );

  void set_parameters( const std::map<std::string, double>& params );
  void set_vehicle_parameters( const dynamics::PhysicalVehicleParameters& params );
  void set_comfort_settings( const std::shared_ptr<dynamics::ComfortSettings>& settings );

private:

  struct SolverParams
  {
    double max_iterations = 1000;
    double tolerance      = 1e-3;
    double max_ms         = 50;
    double debug          = 1.0;
  } solver_params;

  struct PlannerCostWeights
  {
    double lane_error     = 5.0;
    double long_error     = 0.1;
    double speed_error    = 5.0;
    double heading_error  = 10.0;
    double steering_angle = 1.0;
    double acceleration   = 0.1;
  } weights;

  double dt              = 0.1;
  size_t horizon_steps   = 40;
  double ref_traj_length = 100;

  std::shared_ptr<mas::OCP> problem;
  dynamics::Trajectory      reference_trajectory;
  dynamics::Trajectory      guess_trajectory;

  dynamics::VehicleStateDynamic start_state;

  dynamics::PhysicalVehicleParameters        vehicle_params;
  std::shared_ptr<dynamics::ComfortSettings> comfort_settings;

  void                   setup_problem();
  mas::StageCostFunction make_trajectory_cost( const dynamics::Trajectory& ref_traj );

  mas::MotionModel     get_planning_model( const dynamics::PhysicalVehicleParameters& params );
  dynamics::Trajectory extract_trajectory();
  void                 solve_problem();
};

} // namespace planner
} // namespace adore
