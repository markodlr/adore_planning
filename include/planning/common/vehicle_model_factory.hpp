/********************************************************************************
 * Copyright (c) 2025 Contributors
 *
 * SPDX-License-Identifier: EPL-2.0
 ********************************************************************************/

#pragma once

#include "dynamics/physical_vehicle_parameters.hpp"
#include "multi_agent_solver/multi_agent_solver.hpp"

namespace adore
{
namespace planner
{

class VehicleModelFactory
{
public:

  static mas::MotionModel get_planning_model( const dynamics::PhysicalVehicleParameters& params );
};

} // namespace planner
} // namespace adore
