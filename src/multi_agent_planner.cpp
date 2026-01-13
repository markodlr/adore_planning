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

#include "planning/multi_agent_planner.hpp"
#include "planning/common/vehicle_model_factory.hpp"

#include "planning/common/planning_helpers.hpp"

namespace adore
{
namespace planner
{

void
MultiAgentPlanner::set_parameters( const std::map<std::string, double>& params )
{
  for( const auto& [name, value] : params )
  {
    if( name == "dt" && value > 0 ) // Ensure dt > 0
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
    if( name == "max_iterations" )
      solver_params.max_iterations = value;
    if( name == "tolerance" )
      solver_params.tolerance = value;
    if( name == "max_ms" )
      solver_params.max_ms = value;
    if( name == "debug" )
      solver_params.debug = value;
    if( name == "max_lateral_acceleration" )
      max_lateral_acceleration = value;
    if( name == "idm_time_headway" )
      idm_time_headway = value;
    if( name == "desired_distance" )
      desired_distance = value;
  }
}

dynamics::TrafficParticipantSet
MultiAgentPlanner::plan_all_participants( const dynamics::TrafficParticipantSet& in_participants,
                                          const std::shared_ptr<map::Map>&       in_local_map )
{
  traffic_participants = in_participants;
  local_map            = in_local_map;

  solve_problem();
  extract_trajectories();
  return traffic_participants;
}

#include "planning/common/vehicle_model_factory.hpp"

mas::MotionModel
MultiAgentPlanner::get_planning_model( const dynamics::PhysicalVehicleParameters& params )
{
  return VehicleModelFactory::get_planning_model( params );
}

mas::OCP
MultiAgentPlanner::create_single_ocp( size_t id )
{

  auto& participant = traffic_participants.participants.at( id );

  mas::OCP problem;
  problem.state_dim     = 5;
  problem.control_dim   = 2;
  problem.horizon_steps = horizon_steps;
  problem.dt            = dt;
  problem.initial_state = Eigen::VectorXd( problem.state_dim );
  problem.dynamics      = get_planning_model( participant.physical_parameters );

  // set bounds
  Eigen::VectorXd lower_bounds( problem.control_dim ), upper_bounds( problem.control_dim );
  lower_bounds << -participant.physical_parameters.steering_angle_max, participant.physical_parameters.acceleration_min;
  upper_bounds << participant.physical_parameters.steering_angle_max, participant.physical_parameters.acceleration_max;
  problem.input_lower_bounds = lower_bounds;
  problem.input_upper_bounds = upper_bounds;

  Eigen::VectorXd state_lower_bounds( problem.state_dim ), state_upper_bounds( problem.state_dim );
  state_lower_bounds << -100000000, -100000000, -2 * M_PI, -0.2, -0.5; // x, y, yaw_angle, speed, s
  state_upper_bounds << 100000000, 100000000, 2 * M_PI, 15.0, 1000;
  problem.state_lower_bounds = state_lower_bounds;
  problem.state_upper_bounds = state_upper_bounds;

  problem.stage_cost = [&, weights = weights, participant = participant, id = id]( const mas::State& x, const mas::State& u,
                                                                                   size_t time_idx ) -> double {
    double cost        = 0;
    double guide_speed = 8.0;
    // lane deviation cost
    if( participant.route )
    {
      double s               = x( 4 );
      auto   ref             = participant.route->get_pose_at_s( s );
      double curvature       = participant.route->get_curvature_at_s( s );
      auto   map_point       = participant.route->get_map_point_at_s( s );
      double max_curve_speed = max_speed;
      if( curvature > 0 )
        max_curve_speed = std::sqrt( max_lateral_acceleration / curvature );

      guide_speed = std::min( max_speed, max_curve_speed );
      if( map_point.max_speed )
        guide_speed = std::min( guide_speed, map_point.max_speed.value() );

      const double dx = x( 0 ) - ref.x;
      const double dy = x( 1 ) - ref.y;

      const double cref_yaw = math::fast_cos( ref.yaw );
      const double sref_yaw = math::fast_sin( ref.yaw );

      const double lon_err = dx * cref_yaw + dy * sref_yaw;
      const double lat_err = -dx * sref_yaw + dy * cref_yaw;

      cost += weights.lane_error * lat_err * lat_err;
      cost += pow( math::normalize_angle( ref.yaw - x( 2 ) ), 2 ) * weights.heading_error;
    }

    double       heading = x( 2 );
    const double cyaw    = math::fast_cos( heading );
    const double syaw    = math::fast_sin( heading );

    double max_dist   = 1000.0;
    double speed_obj  = 0.0;
    double other_size = 0.0;

    for( const auto& agent : multi_agent_problem.agents )
    {

      auto other_id      = agent->id;
      auto other_ocp_ptr = agent->ocp;

      if( other_id == id )
        continue;

      const auto& other_state = other_ocp_ptr->best_states.col( time_idx );

      const double rel_x = other_state( 0 ) - x( 0 );
      const double rel_y = other_state( 1 ) - x( 1 );

      const double long_dist = rel_x * cyaw + rel_y * syaw;
      const double lat_dist  = -rel_x * syaw + rel_y * cyaw;

      if( long_dist > 2 && std::fabs( lat_dist ) < 3 && long_dist < max_dist )
      {
        max_dist                = long_dist;
        speed_obj               = other_state( 3 );
        auto& other_participant = traffic_participants.participants.at( id );
        other_size              = other_participant.physical_parameters.get_total_length();
      }
    }
    constexpr double DEFAULT_ROUTE_DIST = 1000;

    double remaining_route_length = participant.route ? participant.route->get_length() - x( 4 ) : DEFAULT_ROUTE_DIST;

    const double to_front = participant.physical_parameters.wheelbase + participant.physical_parameters.front_axle_to_front_border;

    const double spacing = desired_distance + to_front + other_size / 2;

    // double idm_acc = idm::calculate_idm_acc( remaining_route_length, max_dist, guide_speed, idm_time_headway, spacing, x( 3 ),
    //                                          participant.physical_parameters.acceleration_max, speed_obj );
    double idm_acc = idm::calculate_idm_acc( remaining_route_length, max_dist, guide_speed, idm_time_headway, spacing, x( 3 ),
                                             participant.physical_parameters.acceleration_max, speed_obj );

    cost += weights.speed_error * std::pow( idm_acc - u( 1 ), 2 );
    cost += weights.steering_angle * u( 0 ) * u( 0 );
    return cost;
  };

  auto& state = participant.state;

  double s = participant.route ? participant.route->get_s( state, 10 ).value_or( 0.0 ) : 0;

  problem.initial_state << state.x, state.y, state.yaw_angle, state.vx, s;

  // populate initial controls from participant trajectory if available
  if( participant.trajectory && !participant.trajectory->states.empty() )
  {
    problem.initial_controls = mas::ControlTrajectory::Zero( problem.control_dim, problem.horizon_steps );
    for( size_t t = 0; t < problem.horizon_steps; ++t )
    {
      if( t < participant.trajectory->states.size() )
      {
        const auto& traj_state           = participant.trajectory->states[t];
        problem.initial_controls( 0, t ) = traj_state.steering_angle;
        problem.initial_controls( 1, t ) = ( traj_state.ax );
      }
    }
  }
  // -------------------------
  // Fast single-agent trajectory to warm-start controls
  // -------------------------
  dynamics::Trajectory fast_traj;

  const bool have_existing_traj = ( participant.trajectory && !participant.trajectory->states.empty() );

  if( !have_existing_traj && participant.route && !participant.route->reference_line.empty() )
  {
    // 1) Build waypoint list along the route ahead of the vehicle
    std::vector<map::MapPoint> points;
    const double               s0    = s;
    const double               s_max = s0 + static_cast<double>( horizon_steps ) * dt * max_speed;

    for( const auto& [s_ref, mp] : participant.route->reference_line )
    {
      if( s_ref + 0.5 < s0 )
        continue;
      if( s_ref > s_max )
        break;
      points.push_back( mp );
    }

    if( points.size() >= 2 )
    {
      // 2) Build simple physical model (kinematic bicycle)
      dynamics::PhysicalVehicleModel model;
      model.params       = participant.physical_parameters;
      model.motion_model = [params = model.params]( const dynamics::VehicleStateDynamic& x, const dynamics::VehicleCommand& u ) {
        return dynamics::kinematic_bicycle_model( x, params, u );
      };

      // 3) Other traffic as obstacles (remove self)
      dynamics::TrafficParticipantSet others = traffic_participants;
      others.participants.erase( id );

      // 4) Target speed for the warm-start
      const double guess_speed = std::max( 1.0, std::min( max_speed, state.vx > 0.0 ? state.vx : 8.0 ) );

      fast_traj = waypoints_to_trajectory( state, points, others, model, guess_speed,
                                           dt // same dt as planner
      );
    }
  }

  // -------------------------
  // Initialize controls
  // -------------------------
  if( !fast_traj.states.empty() )
  {
    // Prefer fast single-agent warm-start if available
    problem.initial_controls = mas::ControlTrajectory::Zero( problem.control_dim, problem.horizon_steps );

    const std::size_t steps = std::min<std::size_t>( problem.horizon_steps, fast_traj.states.size() );

    for( std::size_t t = 0; t < steps; ++t )
    {
      const auto& st                   = fast_traj.states[t];
      problem.initial_controls( 0, t ) = st.steering_angle;
      problem.initial_controls( 1, t ) = st.ax;
    }
  }
  else if( have_existing_traj )
  {
    // Fall back to existing trajectory on the participant, if provided
    problem.initial_controls = mas::ControlTrajectory::Zero( problem.control_dim, problem.horizon_steps );

    for( std::size_t t = 0; t < problem.horizon_steps; ++t )
    {
      if( t < participant.trajectory->states.size() )
      {
        const auto& traj_state           = participant.trajectory->states[t];
        problem.initial_controls( 0, t ) = traj_state.steering_angle;
        problem.initial_controls( 1, t ) = traj_state.ax;
      }
    }
  }


  problem.initialize_problem();
  return problem;
}

void
MultiAgentPlanner::solve_problem()
{
  multi_agent_problem.agents.clear();

  // build OCPs and add as agents
  for( auto& [id, participant] : traffic_participants.participants )
  {
    auto ocp = std::make_shared<mas::OCP>( create_single_ocp( id ) );
    multi_agent_problem.add_agent( std::make_shared<mas::Agent>( id, ocp ) );
  }

  // Nash line‑search with OSQP collocation solver
  mas::SolverParams inner_params{
    { "max_iterations",  100 },
    {      "tolerance", 1e-2 },
    {         "max_ms",   20 },
    {          "debug",  0.0 }
  };
  size_t        max_outer_iterations = 4;
  mas::Solver   solver{ std::in_place_type<mas::iLQR> };
  mas::Strategy strat = mas::TrustRegionNashStrategy{ max_outer_iterations, std::move( solver ), inner_params };
  solution            = mas::solve( strat, multi_agent_problem );
}

void
MultiAgentPlanner::extract_trajectories()
{
  // problem.blocks are sorted by agent id
  for( std::size_t i = 0; i < multi_agent_problem.blocks.size(); ++i )
  {
    std::size_t id          = multi_agent_problem.blocks[i].agent_id;
    auto&       participant = traffic_participants.participants.at( id );
    const auto& X           = solution.states[i];
    const auto& U           = solution.controls[i];

    dynamics::Trajectory traj;
    traj.states.reserve( X.cols() );
    for( int k = 0; k < X.cols(); ++k )
    {
      dynamics::VehicleStateDynamic st;
      st.x              = X( 0, k );
      st.y              = X( 1, k );
      st.yaw_angle      = X( 2, k );
      st.vx             = X( 3, k );
      st.time           = participant.state.time + k * dt;
      st.steering_angle = U( 0, k );
      st.ax             = U( 1, k );
      traj.states.push_back( st );
    }
    participant.trajectory = std::move( traj );
  }
}

} // namespace planner
} // namespace adore
