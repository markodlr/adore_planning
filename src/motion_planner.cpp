#include "planning/motion_planner.hpp"

#include <iostream>
#include <vector>

#include "adore_math/spline.h" // Needed for tk::spline

#include "planning/common/planning_helpers.hpp"

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
  config.drivable_area_before = 0.0;

  // Speed Profile defaults
  config.sp_config.total_time        = 5.0;
  config.sp_config.s_horizon         = 100.0;
  config.sp_config.ds_dp             = 0.1; // Finer grid to allow acceleration from rest with low a_max
  config.sp_config.projection_window = 50.0;

  config.sp_config.dt_qp = 0.1;
  config.sp_config.dt_dp = 0.5;

  config.sp_config.w_dp_progress           = 10.0;
  config.sp_config.w_dp_speed              = 1.0;
  config.sp_config.w_dp_accel              = 10.0;
  config.sp_config.w_dp_jerk               = 10.0;
  config.sp_config.w_dp_obstacle_proximity = 10.0;

  config.sp_config.max_curvature           = 0.1;
  config.sp_config.obstacle_lateral_buffer = 0.5; // [m] reduced from 1.5 default
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
  // std::cout << "Planning speed profile with max_speed=" << comfort_settings->max_speed << ", body_width=" << vehicle_params.body_width
  //           << std::endl;
  auto speed_profile = plan_speed_profile( drivable_area, prediction_participants, ego_state, vehicle_params, *comfort_settings,
                                           config.sp_config );

  if( speed_profile.empty() )
  {
    std::cerr << "Speed profile EMPTY!" << std::endl;
    // We can't proceed without a speed profile for the optimizer
    result.success = false;
    result.message = "Speed Profile EMPTY";
    return result;
  }

  // 3b. Smoothing / Tunneling
  if( !previous_speed_profile_.empty() )
  {
    apply_temporal_smoothing( speed_profile, previous_speed_profile_, config.sp_config );
  }

  // Store for next cycle
  previous_speed_profile_ = speed_profile;

  // 4. Optimization
  if( !drivable_area.reference_line.empty() )
  {
    double dt      = config.sp_config.dt_qp;
    size_t horizon = static_cast<size_t>( config.sp_config.total_time / dt );

    // Generate reference trajectory from speed profile
    dynamics::Trajectory ref_traj = generate_reference_trajectory( speed_profile, drivable_area.reference_line, dt, horizon );

    dynamics::Trajectory traj = optimizer.optimize_trajectory( ego_state, ref_traj );

    traj.adjust_start_time( ego_state.time );
    result.trajectory = std::move( traj );
  }


  result.success = true;
  result.message = "Planning Successful";

  if( result.trajectory.has_value() && !result.trajectory->states.empty() )
  {
    // std::cerr << "optimize_trajectory SUCCESS. Path points: " << result.trajectory->states.size()
    //           << " v_end=" << result.trajectory->states.back().vx << std::endl;
  }
  // fallback if optimization did nothing? optimize_trajectory returns result.

  return result;
}

} // namespace planner
} // namespace adore
