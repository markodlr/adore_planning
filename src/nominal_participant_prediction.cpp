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

#include "planning/nominal_participant_prediction.hpp"

#include <algorithm>
#include <limits>

#include "adore_map/route.hpp"
#include "adore_math/fast_trig.h"
#include "adore_math/point.h"
#include "adore_math/pose.h" // or whatever header defines math::Pose2d

namespace adore
{
namespace planner
{

void
NominalParticipantPrediction::set_parameters( const std::map<std::string, double>& params )
{
  for( const auto& [name, value] : params )
  {
    if( name == "dt" )
      dt = value;
    else if( name == "number_of_integration_steps" )
      number_of_integration_steps = static_cast<int>( value );
    else if( name == "include_initial_state" )
      include_initial_state = ( value != 0.0 );
  }

  if( dt <= 0.0 )
  {
    dt = 0.1;
  }
  if( number_of_integration_steps <= 0 )
  {
    number_of_integration_steps = 1;
  }
  if( default_speed < 0.0 )
  {
    default_speed = 0.0;
  }
}

dynamics::VehicleStateDynamic
NominalParticipantPrediction::get_current_state( const dynamics::TrafficParticipant& participant )
{
  if( participant.trajectory && !participant.trajectory->states.empty() )
  {
    return participant.trajectory->states.back();
  }
  return participant.state;
}

void
NominalParticipantPrediction::plan_trajectories( dynamics::TrafficParticipantSet& traffic_participant_set ) const
{
  // Reset trajectories and optionally seed with current state.
  for( auto& [id, participant] : traffic_participant_set.participants )
  {
    participant.trajectory = dynamics::Trajectory();
    if( include_initial_state )
    {
      participant.trajectory->states.push_back( participant.state );
    }
  }
  constexpr double min_speed = 0.2;
  // Predict each participant independently.
  for( auto& [id, participant] : traffic_participant_set.participants )
  {
    const bool has_route = participant.route && !participant.route->reference_line.empty();

    // Start prediction from the last state in the trajectory (either initial or already filled).
    dynamics::VehicleStateDynamic current_state;
    if( participant.trajectory && !participant.trajectory->states.empty() )
    {
      current_state = participant.trajectory->states.back();
    }
    else
    {
      current_state = participant.state;
    }

    // Constant velocity for the entire horizon.
    double vx_const = current_state.vx;

    if( vx_const <= min_speed )
    {
      // Nothing meaningful to predict.
      continue;
    }

    // Route-related state (only used if has_route).
    double s_current    = 0.0;
    double route_length = 0.0;

    if( has_route )
    {
      auto& route  = *participant.route;
      route_length = route.get_length();

      if( route_length <= 0.0 )
      {
        // Degenerate route, treat as no-route below.
      }
      else
      {
        // Single get_s call per participant.
        s_current = route.get_s( current_state ).value_or( 0.0 );
      }
    }

    for( int step = 0; step < number_of_integration_steps; ++step )
    {
      dynamics::VehicleStateDynamic next_state = current_state;

      if( has_route && route_length > 0.0 )
      {
        auto& route = *participant.route;

        const double ds = vx_const * dt;
        if( ds <= 0.0 )
        {
          break;
        }

        s_current                    = std::min( s_current + ds, route_length );
        const math::Pose2d pose_next = route.get_pose_at_s( s_current );

        next_state.x         = pose_next.x;
        next_state.y         = pose_next.y;
        next_state.yaw_angle = pose_next.yaw;
        next_state.vx        = vx_const;
        next_state.ax        = 0.0;
      }
      else
      {
        // No route: straight-line motion at constant speed in current yaw direction.
        const double dx = vx_const * dt * math::fast_cos( current_state.yaw_angle );
        const double dy = vx_const * dt * math::fast_sin( current_state.yaw_angle );

        next_state.x  += dx;
        next_state.y  += dy;
        next_state.vx  = vx_const;
        next_state.ax  = 0.0;
      }

      // Validity area check.
      if( traffic_participant_set.validity_area.has_value() && !traffic_participant_set.validity_area.value().point_inside( next_state ) )
      {
        // Stop predicting once we leave the validity area.
        break;
      }

      participant.trajectory->states.push_back( next_state );
      current_state = next_state;

      // If we reached the end of the route, no need to continue for this participant.
      if( has_route && route_length > 0.0 && s_current >= route_length )
      {
        break;
      }
    }
  }
}

} // namespace planner
} // namespace adore
