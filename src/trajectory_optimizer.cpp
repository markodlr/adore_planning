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

#include "planning/trajectory_optimizer.hpp"

#include <cmath>

#include <algorithm>
#include <limits>

#include "adore_math/fast_trig.h"

#include "controllers/pure_pursuit.hpp"
#include "dynamics/motion_models.hpp"
#include "planning/common/planning_helpers.hpp"
#include "planning/common/vehicle_model_factory.hpp"

namespace adore
{
namespace planner
{

void
TrajectoryOptimizer::set_parameters( const std::map<std::string, double>& params )
{
  for( const auto& [name, value] : params )
  {
    if( name == "dt" && value > 0 )
      dt = value;
    if( name == "horizon_steps" && value > 0 )
      horizon_steps = static_cast<size_t>( value );
    if( name == "lane_error" )
      weights.lane_error = value;
    if( name == "speed_error" )
      weights.speed_error = value;
    if( name == "heading_error" )
      weights.heading_error = value;
    if( name == "steering_angle" )
      weights.steering_angle = value;
    if( name == "acceleration" )
      weights.acceleration = value;
    if( name == "progress_error" )
      weights.progress_error = value;
    if( name == "max_iterations" )
      solver_params.max_iterations = value;
    if( name == "tolerance" )
      solver_params.tolerance = value;
    if( name == "max_ms" )
      solver_params.max_ms = value;
    if( name == "debug" )
      solver_params.debug = value;
  }
}

void
TrajectoryOptimizer::set_comfort_settings( const std::shared_ptr<dynamics::ComfortSettings>& settings )
{
  comfort_settings = settings;
  comfort_settings->clamp( vehicle_params );
}

void
TrajectoryOptimizer::set_vehicle_parameters( const dynamics::PhysicalVehicleParameters& params )
{
  vehicle_params = params;
}

// -----------------------------------------------------------------------------
// Stage Cost: Penalize errors relative to reference trajectory
// -----------------------------------------------------------------------------
mas::StageCostFunction
TrajectoryOptimizer::make_trajectory_cost()
{
  return [this]( const mas::State& state_vec, const mas::Control& u, std::size_t k ) -> double {
    double cost = 0.0;

    // State: [x, y, yaw, v]
    const double x   = state_vec( 0 );
    const double y   = state_vec( 1 );
    const double yaw = state_vec( 2 );
    const double v   = state_vec( 3 );

    // Get reference state at step k
    // If k exceeds reference size, clamp to last
    size_t      ref_idx   = std::min( k, reference_trajectory.states.size() - 1 );
    const auto& ref_state = reference_trajectory.states[ref_idx];

    double rx   = ref_state.x;
    double ry   = ref_state.y;
    double ryaw = ref_state.yaw_angle;
    double rv   = ref_state.vx;

    // Global Errors
    double dx = x - rx;
    double dy = y - ry;

    // Lateral error: Perpendicular distance to reference heading
    // rotate (dx, dy) by -ryaw
    // y_local = -dx * sin(ryaw) + dy * cos(ryaw)
    double lat_err = -dx * std::sin( ryaw ) + dy * std::cos( ryaw );

    // Heading error
    double hdg_err = math::normalize_angle( yaw - ryaw );

    // Speed Error
    const double spd_err = v - rv;

    // Longitudinal Error (optional "progress" error in local frame)
    // x_local = dx * cos(ryaw) + dy * sin(ryaw)
    // We can penalize this to ensure we stick to the time-parameterized points
    double long_err = dx * std::cos( ryaw ) + dy * std::sin( ryaw );

    cost += weights.lane_error * lat_err * lat_err;
    cost += weights.heading_error * hdg_err * hdg_err;
    cost += weights.speed_error * spd_err * spd_err;
    cost += weights.long_error * long_err * long_err; // Using long_error as progress tracking

    // Control costs: [steer, accel]
    // We can also penalize deviation from reference controls if we trust them?
    // For now, penalize absolute effort
    cost += weights.steering_angle * u( 0 ) * u( 0 );
    cost += weights.acceleration * u( 1 ) * u( 1 );

    return cost;
  };
}

// -----------------------------------------------------------------------------
// Main Optimization
// -----------------------------------------------------------------------------
dynamics::Trajectory
TrajectoryOptimizer::optimize_trajectory( const dynamics::VehicleStateDynamic& current_state, const dynamics::Trajectory& ref_traj )
{
  start_state          = current_state;
  reference_trajectory = ref_traj;

  if( reference_trajectory.states.empty() )
  {
    // Fallback or empty result
    return {};
  }

  // Cap horizon to reference trajectory size - don't plan beyond what we have reference for
  const size_t effective_horizon = std::min( horizon_steps, reference_trajectory.states.size() );
  if( effective_horizon < 2 )
  {
    return {};
  }

  // Define start_state_vec [x, y, yaw, v]
  start_state_vec = Eigen::VectorXd( 4 );
  start_state_vec << start_state.x, start_state.y, start_state.yaw_angle, start_state.vx;

  // 1. Setup Problem (with effective horizon)
  setup_problem( effective_horizon );

  // 2. Generate Initial Guess (populates problem->initial_*)
  generate_initial_guess( effective_horizon );

  // 3. Solve
  solve_problem();

  // 4. Extract Result
  auto result = extract_trajectory();


  return result;
}

void
TrajectoryOptimizer::setup_problem( size_t effective_horizon )
{
  problem = std::make_shared<mas::OCP>();

  problem->state_dim     = 4; // [x, y, yaw, v]
  problem->control_dim   = 2; // [steer, accel]
  problem->horizon_steps = effective_horizon;
  problem->dt            = dt;
  problem->initial_state = start_state_vec;

  // Dynamics: Use VehicleModelFactory
  problem->dynamics = VehicleModelFactory::get_planning_model( vehicle_params );

  // Bounds
  // Controls: [steer, accel]
  Eigen::VectorXd lower_bounds( 2 ), upper_bounds( 2 );
  lower_bounds << -vehicle_params.steering_angle_max, vehicle_params.acceleration_min;
  upper_bounds << vehicle_params.steering_angle_max, vehicle_params.acceleration_max;
  problem->input_lower_bounds = lower_bounds;
  problem->input_upper_bounds = upper_bounds;

  problem->stage_cost = make_trajectory_cost();

  Eigen::VectorXd  state_lower( 4 );
  constexpr double k_inf = std::numeric_limits<double>::infinity();
  // x, y, yaw, v
  state_lower << -k_inf, -k_inf, -k_inf, 0.0;
  problem->state_lower_bounds = state_lower;

  Eigen::VectorXd state_upper( 4 );
  state_upper << k_inf, k_inf, k_inf, 50.0; // Higher max speed limit
  problem->state_upper_bounds = state_upper;

  // Initialize standard fields used by solver
  problem->initial_states   = Eigen::MatrixXd::Zero( 4, effective_horizon );
  problem->initial_controls = Eigen::MatrixXd::Zero( 2, effective_horizon );
}

// -----------------------------------------------------------------------------
// Initial Guess Strategies
// -----------------------------------------------------------------------------
void
TrajectoryOptimizer::generate_initial_guess( size_t effective_horizon )
{
  // 1. Construct PhysicalVehicleModel
  dynamics::PhysicalVehicleModel model;
  model.params       = vehicle_params;
  model.motion_model = [params = vehicle_params]( const dynamics::VehicleStateDynamic& state, const dynamics::VehicleCommand& command ) {
    return dynamics::kinematic_bicycle_model( state, params, command );
  };

  // 2. Call helper
  dynamics::Trajectory guess = initial_guess_pure_pursuit( reference_trajectory, start_state, model );

  // 3. Fill problem
  size_t len = std::min( effective_horizon, guess.states.size() );
  for( size_t k = 0; k < len; ++k )
  {
    const auto& s                   = guess.states[k];
    problem->initial_states( 0, k ) = s.x;
    problem->initial_states( 1, k ) = s.y;
    problem->initial_states( 2, k ) = s.yaw_angle;
    problem->initial_states( 3, k ) = s.vx;

    problem->initial_controls( 0, k ) = s.steering_angle;
    problem->initial_controls( 1, k ) = s.ax;
  }

  // Fill remaining if reference is shorter (extrapolate)
  for( size_t k = len; k < effective_horizon; ++k )
  {
    if( k > 0 )
    {
      problem->initial_states.col( k )   = problem->initial_states.col( k - 1 );
      problem->initial_controls.col( k ) = problem->initial_controls.col( k - 1 );

      // Propagate dynamics simple
      mas::State           x_prev      = problem->initial_states.col( k - 1 );
      mas::Control         u_prev      = problem->initial_controls.col( k - 1 );
      mas::StateDerivative x_dot       = problem->dynamics( x_prev, u_prev );
      problem->initial_states.col( k ) = x_prev + x_dot * dt;
    }
  }

  // Force start state match
  problem->initial_states.col( 0 ) = start_state_vec;

  // Ensure the guess is set in the problem
  problem->best_states   = problem->initial_states;
  problem->best_controls = problem->initial_controls;
}

void
TrajectoryOptimizer::solve_problem()
{
  mas::SolverParams params;
  params["max_iterations"] = solver_params.max_iterations;
  params["tolerance"]      = solver_params.tolerance;
  params["max_ms"]         = solver_params.max_ms;
  params["debug"]          = solver_params.debug;

  auto solve_with = [&]( auto&& solver, double max_ms ) {
    params["max_ms"] = max_ms;
    solver.set_params( params );
    solver.solve( *problem );
    problem->update_initial_with_best();
  };

  solve_with( mas::OSQPCollocation{}, 60 );
}

dynamics::Trajectory
TrajectoryOptimizer::extract_trajectory()
{
  dynamics::Trajectory trajectory;
  trajectory.states.reserve( problem->horizon_steps );

  double current_time = start_state.time;

  for( size_t i = 0; i < problem->horizon_steps; ++i )
  {
    dynamics::VehicleStateDynamic state;
    auto                          x_vec = problem->best_states.col( i );   // [x, y, yaw, v]
    auto                          u_vec = problem->best_controls.col( i ); // [steer, accel]

    state.x         = x_vec( 0 );
    state.y         = x_vec( 1 );
    state.yaw_angle = math::normalize_angle( x_vec( 2 ) );
    state.vx        = x_vec( 3 );

    state.time           = current_time + i * dt;
    state.steering_angle = std::clamp( u_vec( 0 ), -vehicle_params.steering_angle_max, vehicle_params.steering_angle_max );
    state.ax             = u_vec( 1 );

    trajectory.states.push_back( state );
  }

  return trajectory;
}

} // namespace planner
} // namespace adore
