#include "planning/motion_planner.hpp"
#include "planning/common/planning_helpers.hpp"
#include <iostream>

namespace adore
{
namespace planner
{

MotionPlanner::MotionPlanner()
{  
  // Drivable Area defaults
  config.da_config.lane_scope             = LaneScope::AllLanes;
  config.da_config.lateral_inflation      = 0.3;
  config.da_config.longitudinal_inflation = 2.5;

  config.drivable_area_length = 100.0;
  config.drivable_area_before = 10.0;

  // Speed Profile defaults
  config.sp_config.total_time        = 5.0;
  config.sp_config.s_horizon         = 100.0;
  config.sp_config.ds_dp             = 0.1; // Finer grid to allow acceleration from rest with low a_max
  config.sp_config.projection_window = 50.0;

  config.sp_config.dt_qp = 0.1;
  config.sp_config.dt_dp = 0.5;

  config.sp_config.w_dp_progress           = 1000.0;
  config.sp_config.w_dp_speed              = 1.0;
  config.sp_config.w_dp_accel              = 10.0;
  config.sp_config.w_dp_jerk               = 10.0;
  config.sp_config.w_dp_obstacle_proximity = 10.0;

  config.sp_config.max_curvature = 0.1;
  config.sp_config.obstacle_lateral_buffer = 0.5; // [m] reduced from 1.5 default

  config.sp_config.w_qp_track_dp        = 0.01;
  config.sp_config.w_qp_accel           = 1.0;
  config.sp_config.w_qp_jerk            = 10.0;
  config.sp_config.w_qp_v0              = 10.0;
  config.sp_config.w_qp_a0              = 2.0;
  config.sp_config.w_qp_limit_violation = 100.0;
}

void
MotionPlanner::init( const MotionPlannerConfig& config_in, const dynamics::PhysicalVehicleParameters& vehicle_params_in,
                     const std::shared_ptr<dynamics::ComfortSettings>& comfort_settings_in )
{
  config           = config_in;
  vehicle_params   = vehicle_params_in;
  comfort_settings = comfort_settings_in;

  // Derive geometric defaults if not overridden
  const double safety_buffer = 0.0; // [m] minimum cushion on each side of the vehicle

  config.sp_config.obstacle_lateral_buffer = safety_buffer;

  optimizer.set_vehicle_parameters( vehicle_params );
  optimizer.set_comfort_settings( comfort_settings );
}

void
MotionPlanner::update_parameters( const std::map<std::string, double>& params )
{
  // Pass relevant parameters to optimizer
  optimizer.set_parameters( params );
  
  // Potential TODO: parsing map to update MotionPlannerConfig fields
}

PlannerResult
MotionPlanner::plan( const adore::map::Route& route, const dynamics::VehicleStateDynamic& ego_state,
                     const dynamics::TrafficParticipantSet& participants )
{
  PlannerResult result;

  const double state_s    = route.get_s( ego_state, 10.0 ).value_or( 0.0 );
  const double da_start_s = state_s - config.drivable_area_before;
  const double da_end_s   = state_s + config.drivable_area_length;

  // 1. Predict participants
  // We make a local copy because prediction modifies the trajectories in place
  auto prediction_participants = participants;
  participant_predictor.plan_trajectories( prediction_participants );

  // 2. Create Drivable Area
  auto drivable_area = create_drivable_area( route, da_start_s, da_end_s, prediction_participants, vehicle_params, config.da_config );
  
  // Store for debugging
  result.drivable_area = drivable_area;

  if( drivable_area.empty() )
  {
    result.success = false;
    result.message = "Empty Drivable Area";
    return result;
  }

  // 3. Speed Profile
  std::cout << "Planning speed profile with max_speed=" << comfort_settings->max_speed 
            << ", body_width=" << vehicle_params.body_width << std::endl;
  auto speed_profile = plan_speed_profile( drivable_area, prediction_participants, ego_state, vehicle_params, *comfort_settings, config.sp_config );

  if( speed_profile.empty() ) {
      std::cout << "Speed profile EMPTY!" << std::endl;
  } else {
      std::cout << "Speed profile size: " << speed_profile.size() 
                << ", v_start: " << speed_profile.front().v 
                << ", v_max: " << speed_profile.back().v << std::endl;
  }

  if( speed_profile.empty() ) {
      // It's possible to fail finding a valid speed profile
      // In that case, we might want to return failure, or maybe the existing code just continued with empty profile?
      // Checking old code: "std::cerr << "Speed profile points: " << speed_profile.size() << std::endl;"
      // Then "generate_reference_trajectory" handles empty profile by doing nothing.
      // So we can proceed.
  }

  // 4. Reference Trajectory
  const double dt            = 0.1;
  const size_t horizon_steps = 40;

  const dynamics::Trajectory ref_traj = generate_reference_trajectory( speed_profile, drivable_area.reference_line, dt, horizon_steps );

  if (ref_traj.states.empty()) {
      result.success = false;
      result.message = "Failed to generate reference trajectory";
      return result;
  }

  // 5. Initial Guess (Pure Pursuit)
  dynamics::PhysicalVehicleModel temp_model; // Need a model wrapper for pure pursuit helper
  temp_model.params = vehicle_params;
  // (Assuming default motion_model in temp_model is acceptable for initial guess, or we should copy it if available)

  const dynamics::Trajectory guess_traj = initial_guess_pure_pursuit( ref_traj, ego_state, temp_model );

  // 6. Optimization
  dynamics::Trajectory traj = ref_traj; // fallback if optimization skipped (e.g. short ref)

  if( ref_traj.states.size() > 4 )
  {
    traj = optimizer.optimize_trajectory( ego_state, ref_traj, guess_traj );
  }

  traj.adjust_start_time( ego_state.time );
  
  result.trajectory = std::move( traj );
  result.success    = true;
  result.message    = "Planning Successful";

  if (result.trajectory.has_value() && !result.trajectory->states.empty()) {
    std::cerr << "solve_dp SUCCESS. Path points: " << result.trajectory->states.size() << " v_end=" << result.trajectory->states.back().vx << std::endl;
  }
  return result;
}

} // namespace planner
} // namespace adore
