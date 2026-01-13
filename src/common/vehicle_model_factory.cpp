/********************************************************************************
 * Copyright (c) 2025 Contributors
 *
 * SPDX-License-Identifier: EPL-2.0
 ********************************************************************************/

#include "planning/common/vehicle_model_factory.hpp"

#include <cmath>

namespace adore
{
namespace planner
{

mas::MotionModel
VehicleModelFactory::get_planning_model( const dynamics::PhysicalVehicleParameters& params )
{
  return [params]( const mas::State& x, const mas::Control& u ) -> mas::StateDerivative {
    mas::StateDerivative dxdt;
    dxdt.setZero( 4 );
    dxdt( 0 ) = x( 3 ) * std::cos( x( 2 ) );
    dxdt( 1 ) = x( 3 ) * std::sin( x( 2 ) );
    dxdt( 2 ) = x( 3 ) * std::tan( u( 0 ) ) / params.wheelbase;
    dxdt( 3 ) = u( 1 );
    return dxdt;
  };
}

} // namespace planner
} // namespace adore
