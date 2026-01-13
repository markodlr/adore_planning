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

#include <deque>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

#include "adore_math/spline.h"

#include "controllers/pure_pursuit.hpp"
#include "dynamics/integration.hpp"
#include "dynamics/physical_vehicle_model.hpp"
#include "dynamics/vehicle_state.hpp"
#include "planning/idm.hpp"
#include "planning/speed_profile.hpp"
#include <eigen3/Eigen/Dense>

namespace adore
{
namespace planner
{
// Checks if the point is to the right of the first segment in the line
template<typename Point, typename Line>
bool
is_point_to_right_of_line( const Point& point, const Line& line )
{
  if( line.size() < 2 )
  {
    throw std::invalid_argument( "Line must have at least two points to define a direction." );
  }

  // Define the start and end points of the first segment
  const auto& p1 = line[0];
  const auto& p2 = line[1];

  // Calculate the direction vector of the first line segment
  double line_dx = p2.x - p1.x;
  double line_dy = p2.y - p1.y;

  // Calculate the vector from p1 to the point in question
  double point_dx = point.x - p1.x;
  double point_dy = point.y - p1.y;

  // Calculate the cross product of the line direction and point vector
  double cross_product = line_dx * point_dy - line_dy * point_dx;

  // If cross_product is negative, the point is to the right of the first line segment
  return cross_product < 0;
}

template<typename Line, typename Pose>
Line
filter_points_in_front( const Line& line, const Pose& reference_pose )
{
  Line filtered_points;

  // Calculate unit vector for the reference orientation
  double ref_dx = std::cos( reference_pose.yaw_angle );
  double ref_dy = std::sin( reference_pose.yaw_angle );

  for( const auto& point : line )
  {
    // Vector from reference position to the current point
    double dx = point.x - reference_pose.x;
    double dy = point.y - reference_pose.y;

    // Dot product to determine if the point is in front of the reference pose
    double dot_product = dx * ref_dx + dy * ref_dy;
    if( dot_product > 0 )
    { // If dot product is positive, the point is in front
      filtered_points.push_back( point );
    }
  }

  return filtered_points;
}

template<typename Line>
Line
shift_points_right( const Line& points, double shift_distance )
{
  // Check if we have enough points to calculate headings
  if( points.size() < 2 )
    return points;

  Line shifted_points;

  for( size_t i = 0; i < points.size(); ++i )
  {
    double dx, dy;

    if( i == points.size() - 1 )
    {
      // For the last point, use the direction to the previous point
      dx = points[i].x - points[i - 1].x;
      dy = points[i].y - points[i - 1].y;
    }
    else
    {
      // For others, direction to next point
      dx = points[i + 1].x - points[i].x;
      dy = points[i + 1].y - points[i].y;
    }

    // Normalize direction vector
    double length  = std::sqrt( dx * dx + dy * dy );
    dx            /= length;
    dy            /= length;

    // Calculate perpendicular direction to the right
    double shift_dx = dy * shift_distance;
    double shift_dy = -dx * shift_distance;

    // Apply shift to the current point
    auto shifted_point  = points[i];
    shifted_point.x    += shift_dx;
    shifted_point.y    += shift_dy;

    // Add shifted point to the result deque
    shifted_points.emplace_back( shifted_point );
  }

  return shifted_points;
}

// Helper for waypoints_to_trajectory - Calculate the distance error and direction error
static std::pair<double, double>
calculate_errors( const dynamics::VehicleStateDynamic& state, double target_x, double target_y, double target_yaw )
{
  double dx            = target_x - state.x;
  double dy            = target_y - state.y;
  double lateral_error = dy * std::cos( state.yaw_angle ) - dx * std::sin( state.yaw_angle );

  double heading_error = target_yaw - state.yaw_angle;
  heading_error        = adore::math::normalize_angle( heading_error );
  return { lateral_error, heading_error };
}

// compute the distance to the nearest object, if within a certain radius from the center of the waypoint lane, considering obstacle fixed
static double
get_distance_to_nearest_obstacle( const tk::spline& waypoint_spline_x, const tk::spline& waypoint_spline_y, const double waypoints_length,
                                  const dynamics::TrafficParticipantSet& traffic_participants )
{
  double ds                     = 0.5; // check step size
  int    number_of_steps        = static_cast<int>( waypoints_length / ds );
  double min_distance_to_object = std::numeric_limits<double>::max();
  double treshold_within_lane   = 1.0;

  for( const auto& [id, participant] : traffic_participants.participants )
  {
    math::Point2d object_position;
    object_position.x = participant.state.x;
    object_position.y = participant.state.y;
    double s          = 0.0;

    for( int i = 0; i < number_of_steps; i++ )
    {
      math::Point2d waypoint_position;
      waypoint_position.x = waypoint_spline_x( s );
      waypoint_position.y = waypoint_spline_y( s );

      double object_distance = adore::math::distance_2d( object_position, waypoint_position );

      if( object_distance < treshold_within_lane && s < min_distance_to_object )
      {
        min_distance_to_object = s;
        break;
      }

      s = ( i + 1 ) * ds;
    }
  }
  return min_distance_to_object;
}

template<typename Line>
inline dynamics::Trajectory
waypoints_to_trajectory( const dynamics::VehicleStateDynamic& start_state, const Line& waypoints,
                         const dynamics::TrafficParticipantSet& traffic_participants, const dynamics::PhysicalVehicleModel& model,
                         double target_speed = 2.0, double dt = 0.1, double k_speed = 0.5, double k_lateral = 0.5, double k_heading = 2.0,
                         double cg_ratio = 0.5 )
{
  dynamics::Trajectory trajectory;
  trajectory.states.push_back( start_state );

  // Initialize splines for waypoints
  tk::spline          spline_x, spline_y;
  std::vector<double> s_vec, x_vec, y_vec;
  double              cumulative_dist = 0.0;

  for( size_t i = 0; i < waypoints.size(); ++i )
  {
    if( i > 0 )
    {
      double dist = adore::math::distance_2d( waypoints[i], waypoints[i - 1] );
      if( dist < 0.01 || std::isnan( dist ) )
        continue;

      cumulative_dist += dist;
    }
    s_vec.push_back( cumulative_dist );
    x_vec.push_back( waypoints[i].x );
    y_vec.push_back( waypoints[i].y );
  }

  spline_x.set_points( s_vec, x_vec );
  spline_y.set_points( s_vec, y_vec );

  dynamics::VehicleStateDynamic current_state = start_state;
  double                        s             = 0.0;

  for( double time = 0; time <= cumulative_dist / target_speed; time += dt )
  {
    // calculate the distance to closest object
    double closest_obstacle_distance = get_distance_to_nearest_obstacle( spline_x, spline_y, s_vec.back(), traffic_participants );

    // Calculate acceleration based on speed error
    double acceleration = idm::calculate_idm_acc( 100, closest_obstacle_distance, target_speed, 3.0, 7.0, current_state.vx, 2.0, 0.0 );

    dynamics::VehicleCommand control;
    control.acceleration = acceleration;

    // Update distance `s` based on current speed and calculated acceleration
    s += current_state.vx * dt + 0.5 * acceleration * dt * dt;

    // Get target position and heading from the splines using `s`
    double target_x   = spline_x( s );
    double target_y   = spline_y( s );
    double target_yaw = std::atan2( spline_y.deriv( 1, s ), spline_x.deriv( 1, s ) );

    // Calculate errors
    auto [lateral_error, heading_error] = calculate_errors( current_state, target_x, target_y, target_yaw );

    // Set steering angle based on direction and lateral error
    control.steering_angle = heading_error * k_heading + lateral_error * k_lateral;

    control.clamp_within_limits( model.params );

    // Update vehicle state using Euler integration
    current_state                = dynamics::integrate_euler( current_state, control, dt, model.motion_model );
    current_state.ax             = control.acceleration;
    current_state.steering_angle = control.steering_angle;

    // Add the updated state to the trajectory
    trajectory.states.push_back( current_state );
  }

  return trajectory;
}

inline dynamics::Trajectory
generate_reference_trajectory( const SpeedProfile& speed_profile, const std::map<double, adore::map::MapPoint>& reference_line, double dt,
                               size_t horizon )
{
  dynamics::Trajectory ref_trajectory;

  if( horizon == 0 || dt <= 0.0 || speed_profile.empty() || reference_line.size() < 2 )
  {
    std::cerr << "reference_line.size(): " << reference_line.size() << ", speed_profile.size(): " << speed_profile.size() << ", dt: " << dt
              << ", horizon: " << horizon << std::endl;
    return ref_trajectory;
  }
  ref_trajectory.states.reserve( horizon );


  auto clamp01 = []( double u ) { return std::max( 0.0, std::min( 1.0, u ) ); };

  auto sample_speed_at_t = [&]( double t_query, SpeedProfilePoint& out ) -> bool {
    if( speed_profile.empty() )
      return false;

    if( t_query <= speed_profile.front().t )
    {
      out   = speed_profile.front();
      out.t = t_query;
      return true;
    }

    if( t_query >= speed_profile.back().t )
    {
      out   = speed_profile.back();
      out.t = t_query;
      return true;
    }

    auto it_hi = std::lower_bound( speed_profile.begin(), speed_profile.end(), t_query,
                                   []( const SpeedProfilePoint& p, double t ) { return p.t < t; } );

    if( it_hi == speed_profile.begin() )
    {
      out   = *it_hi;
      out.t = t_query;
      return true;
    }

    const auto it_lo = std::prev( it_hi );

    const double t0    = it_lo->t;
    const double t1    = it_hi->t;
    const double denom = t1 - t0;
    const double u     = ( std::abs( denom ) < 1e-9 ) ? 0.0 : clamp01( ( t_query - t0 ) / denom );

    out   = *it_lo;
    out.t = t_query;
    out.s = it_lo->s + u * ( it_hi->s - it_lo->s );
    out.v = it_lo->v + u * ( it_hi->v - it_lo->v );
    out.a = it_lo->a + u * ( it_hi->a - it_lo->a );
    return true;
  };

  auto interpolate_pose_at_s = [&]( double s_query, double& x, double& y, double& yaw ) -> bool {
    if( reference_line.size() < 2 )
      return false;

    auto it_hi = reference_line.lower_bound( s_query );
    auto it_lo = it_hi;

    if( it_hi == reference_line.begin() )
    {
      it_lo = it_hi;
      it_hi = std::next( it_hi );
    }
    else if( it_hi == reference_line.end() )
    {
      it_hi = std::prev( reference_line.end() );
      it_lo = std::prev( it_hi );
    }
    else
    {
      it_lo = std::prev( it_hi );
    }

    const double s0 = it_lo->first;
    const double s1 = it_hi->first;

    const auto& p0 = it_lo->second;
    const auto& p1 = it_hi->second;

    const double denom = s1 - s0;
    const double u     = ( std::abs( denom ) < 1e-9 ) ? 0.0 : clamp01( ( s_query - s0 ) / denom );

    x = p0.x + u * ( p1.x - p0.x );
    y = p0.y + u * ( p1.y - p0.y );

    const double dx = p1.x - p0.x;
    const double dy = p1.y - p0.y;
    yaw             = ( std::abs( dx ) < 1e-9 && std::abs( dy ) < 1e-9 ) ? 0.0 : std::atan2( dy, dx );

    return true;
  };

  for( size_t k = 0; k < horizon; ++k )
  {
    const double t = static_cast<double>( k ) * dt;

    SpeedProfilePoint sp{};
    if( !sample_speed_at_t( t, sp ) )
      break;

    double x   = 0.0;
    double y   = 0.0;
    double yaw = 0.0;
    if( !interpolate_pose_at_s( sp.s, x, y, yaw ) )
      break;

    dynamics::VehicleStateDynamic st;
    st.time      = t; // relative time for optimizer
    st.x         = x;
    st.y         = y;
    st.yaw_angle = yaw;
    st.vx        = sp.v;
    st.ax        = sp.a;

    ref_trajectory.states.push_back( st );
  }

  return ref_trajectory;
}

inline dynamics::Trajectory
initial_guess_pure_pursuit( const dynamics::Trajectory& ref_traj, const dynamics::VehicleStateDynamic& start_state,
                            dynamics::PhysicalVehicleModel& model )
{
  dynamics::Trajectory guess_traj;

  if( ref_traj.states.empty() )
  {
    return guess_traj;
  }

  // Infer dt/horizon from the reference (ref is expected to be relative-time, uniform dt).
  const size_t horizon_steps = ref_traj.states.size();
  double       dt            = 0.1;
  if( horizon_steps >= 2 )
  {
    const double dt_ref = ref_traj.states[1].time - ref_traj.states[0].time;
    if( std::isfinite( dt_ref ) && dt_ref > 1e-6 )
      dt = dt_ref;
  }

  guess_traj.states.reserve( horizon_steps );

  controllers::PurePursuit pp;
  pp.model = model;

  // Initialize guess state from the first reference sample (safe default if caller doesn't pass ego).
  dynamics::VehicleStateDynamic st = start_state;
  st.time                          = ref_traj.states.front().time;

  for( size_t k = 0; k < horizon_steps; ++k )
  {
    const auto cmd = pp.get_next_vehicle_command( ref_traj, st );

    st.ax             = cmd.acceleration;
    st.steering_angle = cmd.steering_angle;
    guess_traj.states.push_back( st );
    st = dynamics::integrate_rk4( st, cmd, dt, pp.model.motion_model );
  }

  return guess_traj;
}


} // namespace planner

} // namespace adore
