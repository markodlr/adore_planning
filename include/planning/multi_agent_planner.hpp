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
#include <vector>

#include "adore_map/map.hpp"
#include "adore_map/route.hpp"
#include "adore_math/angles.h"
#include "adore_math/distance.h"
#include "adore_math/fast_trig.h"

#include "dynamics/traffic_participant.hpp"
#include "dynamics/trajectory.hpp"
#include "multi_agent_solver/multi_agent_solver.hpp"

namespace adore
{
namespace planner
{

class MultiAgentPlanner
{


public:

  dynamics::TrafficParticipantSet plan_all_participants( const dynamics::TrafficParticipantSet& traffic_participants,
                                                         const std::shared_ptr<map::Map>&       local_map );

  void set_parameters( const std::map<std::string, double>& params );


private:

  struct SolverParams
  {
    double max_iterations = 30;
    double tolerance      = 1e-3;
    double max_ms         = 50;
    double debug          = 0.0;
  } solver_params;

  struct PlannerCostWeights
  {
    double lane_error     = 1.0;
    double speed_error    = 0.5;
    double heading_error  = 0.3;
    double steering_angle = 0.2;
  } weights;

  double dt                       = 0.1;
  size_t horizon_steps            = 30;
  double idm_time_headway         = 5.0;
  double desired_distance         = 6.0;
  double max_lateral_acceleration = 2.0;

  double max_speed = 10.0;

  dynamics::TrafficParticipantSet traffic_participants;
  std::shared_ptr<map::Map>       local_map;

  mas::MultiAgentProblem multi_agent_problem;
  mas::Solution          solution;

  mas::OCP create_single_ocp( size_t index );

  mas::StageCostFunction make_trajectory_cost( const dynamics::Trajectory& ref_traj );
  mas::MotionModel       get_planning_model( const dynamics::PhysicalVehicleParameters& params );
  void                   extract_trajectories();
  void                   solve_problem();
};

} // namespace planner
} // namespace adore
